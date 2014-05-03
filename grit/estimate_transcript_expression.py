import sys, os
import time
import traceback

import numpy
import scipy

from pysam import Fastafile, Samfile

from itertools import izip, chain
from collections import defaultdict, namedtuple
import Queue

import random

import multiprocessing
from multiprocessing.sharedctypes import RawArray, RawValue
from lib.multiprocessing_utils import Pool, ThreadSafeFile

from files.gtf import load_gtf, Transcript, Gene
from files.reads import fix_chrm_name_for_ucsc

from f_matrix import DesignMatrix
import frequency_estimation
from frag_len import load_fl_dists, FlDist, build_normal_density

import config

import cPickle as pickle

class SharedData(object):
    """Share data across processes.

    """
    def init_processing_data( self, pickled_gene_fnames, 
                              rnaseq_reads, promoter_reads, polya_reads ):
        pickled_gene_fnames.sort(key=lambda x:x[1], reverse=True)
        for gene_id, n_transcripts, fname in pickled_gene_fnames:
            self.gene_fname_mapping[gene_id] = fname
            self.gene_ntranscripts_mapping[gene_id] = n_transcripts
            self.gene_ids.append(gene_id)
            self.gene_cb_queues[gene_id] = self._manager.list()
            self.gene_cb_queues_locks[gene_id] = self._manager.Lock()
        
        self.rnaseq_reads = rnaseq_reads
        self.promoter_reads = promoter_reads
        self.polya_reads = polya_reads
        
        self.num_rnaseq_reads = multiprocessing.Value('i', 0)
        self.num_cage_reads = multiprocessing.Value('i', 0)
        self.num_polya_reads = multiprocessing.Value('i', 0)
        
        return
    
    def get_gene(self, gene_id):
        if self._cached_gene_id == gene_id and self._cached_gene != None:
            return self._cached_gene
        
        # we don't need to lock this because genes can only be
        # set one time
        fname = self.gene_fname_mapping[gene_id]
        with open(fname) as fp:
            gene = pickle.load(fp)

        self._cached_gene_id = gene_id
        self._cached_gene = gene
        
        return gene
        
    def get_design_matrix(self, gene_id):
        if self._cached_fmat_gene_id == gene_id:
            return self._cached_fmat

        with self.output_dict_lock: 
            fname = self.output_dict[(gene_id, 'design_matrices')]
        with open(fname) as fp:
            f_mat = pickle.load(fp)
        self._cached_fmat_gene_id = gene_id
        self._cached_fmat = f_mat
        return f_mat
    
    def set_design_matrix(self, gene_id, f_mat):
        ofname = os.path.join(config.tmp_dir, gene_id + ".fmat")
        # because there's no cache invalidation mechanism, we're only
        # allowed to set the f_mat object once. This also allows us to
        # move the load outside of the lock
        try: assert (gene_id, 'design_matrices') not in self.output_dict
        except:
            assert False, "%s has already had its design matrix set" % gene_id
        
        with open(ofname, "w") as ofp:
            pickle.dump(f_mat, ofp)
        
        with self.output_dict_lock: 
            self.output_dict[(gene_id, 'design_matrices')] = ofname

        if f_mat.num_rnaseq_reads != None:
            with self.num_rnaseq_reads.get_lock():
                self.num_rnaseq_reads.value += f_mat.num_rnaseq_reads
        if f_mat.num_fp_reads != None:
            with self.num_cage_reads.get_lock():
                self.num_cage_reads.value += f_mat.num_fp_reads
        if f_mat.num_tp_reads != None:
            with self.num_polya_reads.get_lock():
                self.num_polya_reads.value += f_mat.num_tp_reads
        
        return
    
    def get_num_reads_in_bams(self):
        return (self.num_cage_reads.value, 
                self.num_rnaseq_reads.value, 
                self.num_polya_reads.value)
    
    def get_mle(self, gene_id):
        with self.output_dict_lock: 
            return numpy.frombuffer(self.mle_estimates[gene_id])
    
    def set_mle(self, gene, mle):
        assert len(mle) == len(gene.transcripts) + 1
        with self.output_dict_lock: 
            self.mle_estimates[gene.id][:] = mle
            
    def get_cbs(self, gene_id, cb_type):
        if cb_type == 'ub':
            return numpy.frombuffer(self.ubs[gene_id])
        elif cb_type == 'lb':
            return numpy.frombuffer(self.lbs[gene_id])
        else: 
            assert False, "Unrecognized confidence bound type '%s'" % cb_type
    
    def set_cbs(self, gene_id, bnd_type_indices_and_values):
        with self.output_dict_lock: 
            for cb_type, index, value in bnd_type_indices_and_values:
                if cb_type == 'ub':
                    data = self.ubs[gene_id][index] = value
                elif cb_type == 'lb':
                    data = self.lbs[gene_id][index] = value
                else: 
                    assert False, "Unrecognized confidence bound type '%s'" % cb_type
        
        return
    
    def __init__(self):
        self.gene_ids = []
        self.gene_fname_mapping = {}
        self.gene_ntranscripts_mapping = {}
        self.gene_cb_queues = {}
        self.gene_cb_queues_locks = {}
        
        self._manager = multiprocessing.Manager()

        # initialize queues to store the expression commands. Since we dont
        # know the total number of reads until the design matrices have been
        # built, we do that first and then start every commands that's been 
        # stored into expression queue        
        self.input_queue = self._manager.dict()
        self.input_queue_lock = multiprocessing.Lock()

        self.cb_genes = self._manager.list()
        self.cb_genes_being_processed = self._manager.list()
        self.cb_genes_lock = multiprocessing.Lock()
        
        # store the expression estimates
        self.mle_estimates = {}
        self.lbs = {}
        self.ubs = {}
        
        # store data that all children need to be able to access        
        self.output_dict = self._manager.dict()
        self.output_dict_lock = multiprocessing.Lock()    
                
        # create objects to cache gene objects, so that we dont have to do a 
        # fresh load from the shared manager
        self._cached_gene_id = None
        self._cached_gene = None
        self._cached_fmat_gene_id = None
        self._cached_fmat = None
    
    def populate_expression_queue(self):
        for gene_id in self.gene_fname_mapping:
            n_trans = self.gene_ntranscripts_mapping[gene_id]
            self.mle_estimates[gene_id] = RawArray(
                'd', [-1]*(n_trans+1))
            self.ubs[gene_id] = RawArray(
                'd', [-1]*n_trans)
            self.lbs[gene_id] = RawArray(
                'd', [-1]*n_trans)
        
    def populate_estimate_cbs_queue(self, gene_id, bnd_type):
        assert bnd_type in ('ub', 'lb')
        f_mat = self.get_design_matrix(gene_id)
        with self.gene_cb_queues_locks[gene_id]:
            queue = self.gene_cb_queues[gene_id]
            # add the out of gene bin
            #queue.append((None, 0, 'ub'))
            #queue.append((None, 0, 'lb'))
            # add 1 to the row num to account for the out of gene bin
            for row_num, t_index in enumerate(f_mat.transcript_indices()):
                queue.append((t_index, row_num+1, bnd_type))
        
        with self.cb_genes_lock:
            self.cb_genes.append(gene_id)
        
        return
    
    def get_cb_gene(self):
        with self.cb_genes_lock:
            try: 
                gene_id = self.cb_genes.pop()
                self.cb_genes_being_processed.append(gene_id)
                return gene_id
            except IndexError:
                # if there are no untouched genes left, then 
                # provide a gene that is being processed by something else
                while len(self.cb_genes_being_processed) > 0:
                    i = random.randint(0, len(self.cb_genes_being_processed)-1)
                    gene_id = self.cb_genes_being_processed[i]
                    with self.gene_cb_queues_locks[gene_id]:
                        if len(self.gene_cb_queues[gene_id]) == 0:
                            del self.cb_genes_being_processed[i]
                        else:
                            return gene_id
                return None
    
    def get_queue_item(self, gene_id=None, 
                       raise_error_if_no_matching_gene=False):
        """Get an item from the work queue.
        
        If gene id is set, then try to get work from this gene to avoid
        having to reload the genes. 
        """
        with self.input_queue_lock:
            if len(self.input_queue) == 0:
                raise Queue.Empty, "Work queue is empty"
            
            if gene_id == None:
                gene_id, gene_work = self.input_queue.popitem()
            else:
                try:
                    gene_work = self.input_queue.pop(gene_id)
                except KeyError:
                    if raise_error_if_no_matching_gene:
                        raise Queue.Empty, "Work queue is empty for gene_id "
                    else:
                        gene_id, gene_work = self.input_queue.popitem()
            
            work_type, work_data = gene_work.pop()
            if len(gene_work) > 0:
                self.input_queue[gene_id] = gene_work

        assert work_type != 'ERROR'
        assert work_type != 'FINISHED'
        return work_type, gene_id, work_data

def calc_fpkm( gene, fl_dists, freqs, num_reads_in_bam):
    assert len(gene.transcripts) == len(freqs)

    num_reads_in_bam = num_reads_in_bam[1]    
    fl_dist = fl_dists['mean']
    
    def calc_effective_length_and_scale_factor(t):
        length = sum( e[1] - e[0] + 1 for e in t.exons ) 
        # subtract for mappability problems at junctions
        length -= 0*len(t.introns)
        if length < fl_dist.fl_min: 
            return 0, 0
        fl_min, fl_max = fl_dist.fl_min, min(length, fl_dist.fl_max)
        allowed_fl_lens = numpy.arange(fl_min, fl_max+1)
        weights = fl_dist.fl_density[
            fl_min-fl_dist.fl_min:fl_max-fl_dist.fl_min+1]
        sc_factor = min(20, 1./weights.sum())
        mean_fl_len = float((allowed_fl_lens*weights).sum())
        return length - mean_fl_len, sc_factor
    
    fpkms = []
    for t, freq in izip( gene.transcripts, freqs ):
        if freq < 0:
            fpkm = None
        else:
            effective_t_len, t_scale = calc_effective_length_and_scale_factor(t)
            if effective_t_len <= 0: 
                fpkms.append( 0 )
                continue
            assert effective_t_len > 0
            assert num_reads_in_bam > 0
            fpk = freq*num_reads_in_bam/(effective_t_len/1000.)
            fpkm = fpk/(num_reads_in_bam/1000000.)
        fpkms.append( fpkm )
    return fpkms

def find_confidence_bounds_in_gene( gene, num_reads_in_bams,
                                    f_mat, mle_estimate, 
                                    trans_indices, trans_indices_lock,
                                    cb_alpha):
    # update the mle_estimate array to only store observable transcripts
    # add 1 to skip the out of gene bin
    observable_trans_indices = (
        numpy.array([-1,] + f_mat.transcript_indices().tolist())+1 )
    mle_estimate = mle_estimate[observable_trans_indices]

    if config.VERBOSE:
        config.log_statement( 
            "Estimating confidence bounds for gene %s" % gene.id )
    
    #n_skipped = sum( 1 for x in sorted(f_mat.filtered_transcripts)
    #                 if x < trans_indices[0])
    # XXX Make sure that this is beign counted correctly
    #n_skipped_tmp = len(set(xrange(trans_indices[0])) - \
    #    set(x-1 for x in observable_trans_indices[1:] if x-1 < trans_indices[0]))
    #config.log_statement( str([n_skipped_tmp, n_skipped, f_mat.filtered_transcripts, \
    #    observable_trans_indices, trans_indices]), log=True)
    #assert n_skipped == n_skipped_tmp
    
    res = []
    while True:
        try:
            with trans_indices_lock:
                trans_index, exp_mat_row, bnd_type = trans_indices.pop()
        except IndexError:
            config.log_statement("")
            return res
        
        if config.DEBUG_VERBOSE: config.log_statement( 
            "Estimating %s confidence bound for gene %s (%i/%i remain)" % ( 
            bnd_type, gene.id, len(trans_indices), 2*len(gene.transcripts)))
        try:
            p_value, bnd = frequency_estimation.estimate_confidence_bound( 
                f_mat, num_reads_in_bams,
                exp_mat_row, mle_estimate, bnd_type, cb_alpha )
        except Exception, inst:
            p_value = 1.
            bnd = 0.0 if bnd_type == 'lb' else 1.0
            error_msg = "%i: Skipping %s (%s:%s:%i-%i): %s" % (
                os.getpid(), gene.id, 
                gene.chrm, gene.strand, gene.start, gene.stop, inst)
            config.log_statement( error_msg, log=True )
            if config.DEBUG_VERBOSE:
                config.log_statement( traceback.format_exc(), log=True )
        
        if config.DEBUG_VERBOSE: config.log_statement( 
            "FINISHED %s BOUND %s\t%s\t%i/%i\t%.2e\t%.2e" % (
            bnd_type, gene.id, None, 
            trans_index, len(gene.transcripts), 
            bnd, p_value ) )
        res.append((bnd_type, trans_index, bnd))

    if config.VERBOSE:
        config.log_statement( 
            "FINISHED Estimating %s confidence bound for gene %s transcript %i-%i/%i" % ( 
            bnd_type,gene.id,trans_indices[0]+1, trans_indices[-1]+1, 
            len(gene.transcripts)))
    
    return res

def find_confidence_bounds_worker( data ):
    while True:
        # get a gene to process
        gene_id = data.get_cb_gene()
        if gene_id == None:
            config.log_statement("")
            return

        if config.VERBOSE:
            config.log_statement(
                "Estimating confidence bounds for '%s'" % gene_id)
        
        gene = data.get_gene(gene_id)
        f_mat = data.get_design_matrix(gene_id)
        mle_estimate = data.get_mle(gene_id)
        trans_indices = data.gene_cb_queues[gene_id]
        trans_indices_lock = data.gene_cb_queues_locks[gene_id]
        num_reads_in_bams = data.get_num_reads_in_bams()
        
        cbs = find_confidence_bounds_in_gene( 
            gene, num_reads_in_bams,
            f_mat, mle_estimate, 
            trans_indices, trans_indices_lock,
            cb_alpha=config.CB_SIG_LEVEL)
        
        data.set_cbs(gene.id, cbs)
        
        if config.VERBOSE:
            config.log_statement("Finished processing '%s'" % gene_id)
    
    return

def estimate_confidence_bounds( data, bnd_type ):
    # sort so that the biggest genes are processed first
    config.log_statement(
        "Populating estimate confidence bounds queue.")

    for gene_id in sorted(data.gene_ids, 
                          key=lambda x:data.gene_ntranscripts_mapping[x]):
        try: data.populate_estimate_cbs_queue(gene_id, bnd_type)
        except KeyError: continue
    
    if config.NTHREADS == 1:
        find_confidence_bounds_worker( data )
    else:
        ps = []
        for i in xrange(config.NTHREADS):
            p = multiprocessing.Process(
                target=find_confidence_bounds_worker, args=[data,])
            p.daemon=True
            p.start()
            ps.append(p)

        while True:
            config.log_statement(
                "Waiting for estimate confidence bound children processes to finish (%i/%i genes remain)"%(
                    len(data.cb_genes), len(data.gene_ids)))
                    
            time.sleep(1)
            if all( not p.is_alive() for p in ps): break
    
    return


def estimate_mle_worker( gene_ids, gene_ids_lock, data ):
    while True:
        config.log_statement("Retrieving gene from quue")
        try:
            with gene_ids_lock:
                gene_id = gene_ids.pop()
        except IndexError:
            config.log_statement("")
            return        
        
        try:
            config.log_statement(
                "Loading gene %s" % gene_id )
            gene = data.get_gene(gene_id)
              
            config.log_statement(
                "Finding MLE for Gene %s(%s:%s:%i-%i) - %i transcripts" \
                    % (gene.id, gene.chrm, gene.strand, 
                       gene.start, gene.stop, len(gene.transcripts) ) )
            
            f_mat = data.get_design_matrix(gene_id)
            num_reads_in_bams = data.get_num_reads_in_bams()

            expected_array, observed_array = f_mat.expected_and_observed(
                num_reads_in_bams)
            mle = frequency_estimation.estimate_transcript_frequencies( 
                observed_array, expected_array)
        except Exception, inst:
            error_msg = "%i: Skipping %s (%s:%s:%i-%i): %s" % (
                os.getpid(), gene.id, 
                gene.chrm, gene.strand, gene.start, gene.stop, inst)
            config.log_statement( error_msg, log=True )
            if config.DEBUG_VERBOSE:
                config.log_statement( traceback.format_exc(), log=True )
            return None

        log_lhd = frequency_estimation.calc_lhd( 
            mle, observed_array, expected_array)

        # add back in the missing trasncripts
        full_mle = -1*numpy.ones(len(gene.transcripts)+1, dtype=float)
        full_mle[numpy.array([-1,]+f_mat.transcript_indices().tolist())+1] = mle
        
        data.set_mle(gene, full_mle)
        config.log_statement( "FINISHED MLE %s\t%.2f - updating queues" % ( 
                gene.id, log_lhd ) )

def estimate_mles( data ):
    config.log_statement("Initializing MLE queue")
    
    manager = multiprocessing.Manager()
    gene_ids = manager.list()
    gene_ids_lock = manager.Lock()
    # sort so that the biggest genes are processed first
    for gene_id in sorted(data.gene_ids, 
                          key=lambda x:data.gene_ntranscripts_mapping[x]):
        gene_ids.append(gene_id)
    
    args = [ gene_ids, gene_ids_lock, data ]
    if config.NTHREADS == 1:
        estimate_mle_worker(*args)
    else:
        ps = []
        for i in xrange(config.NTHREADS):
            p = multiprocessing.Process(
                target=estimate_mle_worker, args=args)
            p.daemon=True
            p.start()
            ps.append(p)

        while True:
            config.log_statement(
                "Waiting for MLE children processes to finish (%i/%i genes remain)"%(
                    len(gene_ids), len(data.gene_ids)))
            time.sleep(1)
            if all( not p.is_alive() for p in ps): break
    
    manager.shutdown()
    return

def build_design_matrices_worker( gene_ids, gene_ids_lock,
                                  data, fl_dists,
                                  (rnaseq_reads, promoter_reads, polya_reads)):
    if rnaseq_reads != None: rnaseq_reads.reload()
    if promoter_reads != None: promoter_reads.reload()
    if polya_reads != None: polya_reads.reload()
    
    while True:
        with gene_ids_lock:
            try: gene_id = gene_ids.pop()
            except IndexError: 
                config.log_statement("")
                return
        
        try:
            gene = data.get_gene(gene_id)
            config.log_statement( 
                "Finding design matrix for Gene %s(%s:%s:%i-%i) - %i transcripts"%(
                    gene.id, gene.chrm, gene.strand, 
                    gene.start, gene.stop, len(gene.transcripts) ) )
            
            f_mat = DesignMatrix(gene, fl_dists, 
                                 rnaseq_reads, promoter_reads, polya_reads,
                                 config.MAX_NUM_TRANSCRIPTS)
        except Exception, inst:
            error_msg = "%i: Skipping %s (%s:%s:%i-%i): %s" % (
                os.getpid(), gene.id, 
                gene.chrm, gene.strand, gene.start, gene.stop, inst)
            if config.DEBUG_VERBOSE:
                config.log_statement( 
                    error_msg + "\n" + traceback.format_exc(), log=True )
        else:
            config.log_statement( "WRITING DESIGN MATRIX TO DISK %s" % gene.id )
            data.set_design_matrix(gene.id, f_mat)
            config.log_statement( "FINISHED DESIGN MATRICES %s" % gene.id )
    
def build_design_matrices( data, fl_dists,
                           (rnaseq_reads, promoter_reads, polya_reads)):    
    manager = multiprocessing.Manager()
    gene_ids = manager.list()
    gene_ids_lock = manager.Lock()
    if config.VERBOSE:
        config.log_statement( "Populating estimate cofidence bounds queue" )
    # sort so that the biggest genes are processed first
    for gene_id in sorted(data.gene_ids, 
                          key=lambda x:data.gene_ntranscripts_mapping[x]):
        gene_ids.append(gene_id)
    if config.VERBOSE:
        config.log_statement( "FINISHED Populating estimate confidence bounds queue" )
    
    args = [ gene_ids, gene_ids_lock, data, fl_dists, 
             (rnaseq_reads, promoter_reads, polya_reads)]
    if config.NTHREADS == 1:
        build_design_matrices_worker(*args)
    else:
        ps = []
        for i in xrange(config.NTHREADS):
            p = multiprocessing.Process(
                target=build_design_matrices_worker, args=args)
            p.daemon=True
            p.start()
            ps.append(p)

        while True:
            config.log_statement(
                "Waiting for design matrix children processes to finish (%i/%i genes remain)"%(
                    len(gene_ids), len(data.gene_ids)))
            time.sleep(1)
            if all( not p.is_alive() for p in ps): break
    
    manager.shutdown()
    return

def write_data_to_tracking_file(data, fl_dists, ofp):
    num_reads_in_bams = data.get_num_reads_in_bams()
    ofp.write("\t".join(
            ["tracking_id", "gene_id ",
             "coverage", "FPKM    ",
             "FPKM_lo ", "FPKM_hi ", "status"] 
            ) + "\n")

    try: 
        sorted_gene_ids = sorted(
            data.gene_ids, key=lambda x: int(x.split("_")[-1]))
    except:
        sorted_gene_ids = data.gene_ids
    for gene_id in sorted_gene_ids:
        gene = data.get_gene(gene_id)
        mles = data.get_mle(gene_id)
        mle_fpkms = calc_fpkm( gene, fl_dists, mles[1:], num_reads_in_bams)
        ubs = data.get_cbs(gene_id, 'ub')
        if ubs != None:
            ub_fpkms = calc_fpkm( gene, fl_dists, ubs, num_reads_in_bams)
        lbs = data.get_cbs(gene_id, 'lb')
        if lbs != None:
            lb_fpkms = calc_fpkm( gene, fl_dists, lbs, num_reads_in_bams)
        try: sorted_transcripts = sorted(gene.transcripts,
                                         key=lambda x: int(x.id.split("_")[-1]))
        except: sorted_transcripts = gene.transcripts
        for i, t in enumerate(sorted_transcripts):
            line = []
            line.append(t.id.ljust(11))
            line.append(t.gene_id.ljust(11))
            line.append('-'.ljust(8))
            if mles == None or mle_fpkms[i] == None: line.append('-       ')
            else: line.append(('%.2e' % mle_fpkms[i]).ljust(8))
            if lbs == None or lb_fpkms[i] == None: line.append('-       ')
            else: line.append(('%.2e' % lb_fpkms[i]).ljust(8))
            if ubs == None or ub_fpkms[i] == None: line.append('-       ')
            else: line.append(('%.2e' % ub_fpkms[i]).ljust(8))
            line.append( "OK" )

            ofp.write("\t".join(line) + "\n" )

def quantify_transcript_expression(
    promoter_reads, rnaseq_reads, polya_reads,
    pickled_gene_fnames, fl_dists,
    ofname, 
    do_estimate_confidence_bounds=True ):
    """Build transcripts
    """
    write_design_matrices=False
        
    data = SharedData()
    if config.VERBOSE: config.log_statement( 
        "Initializing processing data" )
    data.init_processing_data(pickled_gene_fnames, 
                              rnaseq_reads, promoter_reads, polya_reads)
    if config.VERBOSE: config.log_statement( 
        "Building design matrices" )
    build_design_matrices( data, fl_dists,
                           (rnaseq_reads, promoter_reads, polya_reads))
    
    if config.VERBOSE: config.log_statement( 
        "Populating input queue from expression queue" )
    data.populate_expression_queue()
    if config.VERBOSE: config.log_statement( 
        "Estimating MLEs" )
    estimate_mles( data )

    if config.VERBOSE: config.log_statement( 
        "Calculating FPKMS and Writing mle's to output mle" )
    
    if do_estimate_confidence_bounds:
        if config.VERBOSE: config.log_statement( 
            "Estimating lower confidence bounds" )
        estimate_confidence_bounds(data, 'lb')
        if config.VERBOSE: config.log_statement( 
            "FINISHED Estimating lower confidence bounds" )

        if config.VERBOSE: config.log_statement( 
            "Estimating upper confidence bounds" )
        estimate_confidence_bounds(data, 'ub')
        if config.VERBOSE: config.log_statement( 
            "FINISHED Estimating upper confidence bounds" )
    
    if config.VERBOSE: config.log_statement( 
        "Writing output data to tracking file" )
    
    expression_ofp = ThreadSafeFile(ofname, "w")
    write_data_to_tracking_file(data, fl_dists, expression_ofp)    
    expression_ofp.close()
    
    return
