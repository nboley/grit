import sys, os
import time
import traceback

import numpy
import scipy

import copy

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

import f_matrix
import frequency_estimation
from frag_len import load_fl_dists, FlDist, build_normal_density

import config

import cPickle as pickle

SAMPLE_ID = None
REP_ID = None

class NoDesignMatrixError(Exception):
    pass

class SharedData(object):
    """Share data across processes.

    """    
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

        with self.design_mat_lock: 
            fname = self.design_mat_filenames[gene_id].value
        if fname == '': 
            raise NoDesignMatrixError, "No design matrix for '%s'" % gene_id
        with open(fname) as fp:
            f_mat = pickle.load(fp)
        self._cached_fmat_gene_id = gene_id
        self._cached_fmat = f_mat
        return f_mat
    
    def set_design_matrix(self, gene_id, f_mat):
        ofname = config.get_fmat_tmp_fname(gene_id, SAMPLE_ID, REP_ID)
        
        # because there's no cache invalidation mechanism, we're only
        # allowed to set the f_mat object once. This also allows us to
        # move the load outside of the lock
        try: assert self.design_mat_filenames[gene_id].value == ''
        except:
            config.log_statement(
                "%s has already had its design matrix set (%s)" % (
                    gene_id, self.design_mat_filenames[gene_id].value ), 
                log=True)
            return
        
        with open(ofname, "w") as ofp:
            pickle.dump(f_mat, ofp)
        
        with self.design_mat_lock: 
            self.design_mat_filenames[gene_id].value = ofname
        
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
        return numpy.frombuffer(self.mle_estimates[gene_id])
    
    def set_mle(self, gene, mle):
        assert len(mle) == len(gene.transcripts) + 1
        with self.mle_lock: 
            self.mle_estimates[gene.id][:] = mle
            
    def get_cbs(self, gene_id, cb_type):
        if cb_type == 'ub':
            return numpy.frombuffer(self.ubs[gene_id])
        elif cb_type == 'lb':
            return numpy.frombuffer(self.lbs[gene_id])
        else: 
            assert False, "Unrecognized confidence bound type '%s'" % cb_type
    
    def set_cbs(self, gene_id, bnd_type_indices_and_values):
        with self.cbs_lock: 
            for cb_type, index, value in bnd_type_indices_and_values:
                if cb_type == 'ub':
                    data = self.ubs[gene_id][index] = value
                elif cb_type == 'lb':
                    data = self.lbs[gene_id][index] = value
                else: 
                    assert False, "Unrecognized confidence bound type '%s'" % cb_type
        
        return
    
    def __init__(self, pickled_gene_fnames):        
        self._manager = multiprocessing.Manager()
        
        #self.cb_genes = self._manager.list()
        #self.cb_genes_being_processed = self._manager.list()
        self.cb_genes_lock = multiprocessing.Lock()
        
        # store the expression estimates
        self.mle_estimates = {}
        self.lbs = {}
        self.ubs = {}
        
        # store data that all children need to be able to access        
        self.design_mat_filenames = {}
        self.design_mat_lock = multiprocessing.Lock()    
        
        self.mle_lock = multiprocessing.Lock()    
        self.cbs_lock = multiprocessing.Lock()    
        
        # initialize the gene data
        self.gene_ids = []
        self.gene_fname_mapping = {}
        self.gene_ntranscripts_mapping = {}

        pickled_gene_fnames.sort(key=lambda x:x[1], reverse=True)
        for gene_id, n_transcripts, fname in pickled_gene_fnames:
            self.gene_fname_mapping[gene_id] = fname
            self.gene_ntranscripts_mapping[gene_id] = n_transcripts
            self.gene_ids.append(gene_id)
            
            self.design_mat_filenames[gene_id] = multiprocessing.Array(
                'c', 1000)
            self.design_mat_filenames[gene_id].value = ''
        
        self.num_rnaseq_reads = multiprocessing.Value('i', 0)
        self.num_cage_reads = multiprocessing.Value('i', 0)
        self.num_polya_reads = multiprocessing.Value('i', 0)

                
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
                                    trans_indices, cntr,
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
        with cntr.get_lock():
            index = cntr.value
            if index == -1: 
                config.log_statement('')
                break
            cntr.value -= 1
                
        trans_index, exp_mat_row, bnd_type = trans_indices[index]
        
        config.log_statement( 
            "Estimating %s confidence bound for gene %s (%i/%i remain)" % ( 
                bnd_type, gene.id, cntr.value+1, len(gene.transcripts)))
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
            config.log_statement( traceback.format_exc(), log=True )
        
        if config.DEBUG_VERBOSE: config.log_statement( 
            "FINISHED %s BOUND %s\t%s\t%i/%i\t%.2e\t%.2e" % (
            bnd_type, gene.id, None, 
            trans_index, len(gene.transcripts), 
            bnd, p_value ) )
        res.append((bnd_type, trans_index, bnd))

    if config.VERBOSE:
        config.log_statement( 
            "FINISHED Estimating confidence bound for gene %s" % gene.id )
    
    return res

def find_confidence_bounds_worker( 
        data, gene_ids, trans_index_cntrs, bnd_type ):
    def get_new_gene():
        
        # get a gene to process
        try: gene_id = gene_ids.get(timeout=0.1)
        except Queue.Empty: 
            assert gene_ids.qsize() == 0
            config.log_statement("")
            raise IndexError, "No genes left"
        
        config.log_statement(
            "Loading design matrix for gene '%s'" % gene_id)

        gene = data.get_gene(gene_id)
        try: 
            f_mat = data.get_design_matrix(gene_id)
        except NoDesignMatrixError:
            if config.DEBUG_VERBOSE:
                config.log_statement("No design matrix for '%s'" % gene_id, 
                                     log=True)
            raise

        mle_estimate = data.get_mle(gene_id)
        
        trans_indices = []
        for row_num, t_index in enumerate(f_mat.transcript_indices()):
            trans_indices.append((t_index, row_num+1, bnd_type))

        cntr = trans_index_cntrs[gene_id]
        with cntr.get_lock():
            if cntr.value == -1000: 
                cntr.value = len(trans_indices)-1
        
        return gene, f_mat, mle_estimate, trans_indices, cntr
    
    def get_gene_being_processed():
        longest_gene_id = None
        gene_len = 0
        for gene_id, cntr in trans_index_cntrs.iteritems():
            value = cntr.value
            if value > gene_len:
                longest_gene_id = gene_id
                gene_len = value
        
        if longest_gene_id == None: 
            return None
        
        gene = data.get_gene(longest_gene_id)
        f_mat = data.get_design_matrix(longest_gene_id)
        mle_estimate = data.get_mle(longest_gene_id)

        trans_indices = []
        for row_num, t_index in enumerate(f_mat.transcript_indices()):
            trans_indices.append((t_index, row_num+1, bnd_type))
        
        return ( gene, f_mat, mle_estimate, 
                 trans_indices, trans_index_cntrs[longest_gene_id] )
    
    no_new_genes = False    
    num_reads_in_bams = data.get_num_reads_in_bams()
    while True:
        try:
            try: 
                gene, f_mat, mle_estimate, trans_indices, cntr = get_new_gene()
            except NoDesignMatrixError:
                continue
            except IndexError: 
                res = get_gene_being_processed()
                if res == None: 
                    break
                gene, f_mat, mle_estimate, trans_indices, cntr = res
            
            cbs = find_confidence_bounds_in_gene( 
                gene, num_reads_in_bams,
                f_mat, mle_estimate, 
                trans_indices, cntr,
                cb_alpha=config.CB_SIG_LEVEL)
            data.set_cbs(gene.id, cbs)
            
            if config.VERBOSE:
                config.log_statement("Finished processing '%s'" % gene.id)
        except Exception, inst:
            config.log_statement( traceback.format_exc(), log=True )
    
    config.log_statement("")
    return

def estimate_confidence_bounds( data, bnd_type ):
    config.log_statement(
        "Populating estimate confidence bounds queue.")

    ## populate the queue
    # sort so that the biggest genes are processed first
    gene_ids = multiprocessing.Queue()
    trans_index_cntrs = {}
    sorted_gene_ids = sorted(data.gene_ids, 
                             key=lambda x:data.gene_ntranscripts_mapping[x],
                             reverse=True)
    for i, gene_id in enumerate(sorted_gene_ids):
        gene_ids.put(gene_id)
        trans_index_cntrs[gene_id] = multiprocessing.Value( 'i', -1000)
    
    config.log_statement("Waiting on gene bounds children")

    if False and config.NTHREADS == 1:
        find_confidence_bounds_worker( 
            data, gene_ids, trans_indices_queues, bnd_type )
    else:
        pids = []
        for i in xrange(config.NTHREADS):
            pid = os.fork()
            if pid == 0:
                try: 
                    find_confidence_bounds_worker(
                        data, gene_ids, 
                        trans_index_cntrs, bnd_type)
                except Exception, inst:
                    config.log_statement( traceback.format_exc(), log=True )
                finally:
                    os._exit(0)
            pids.append(pid)
        
        for pid in pids:
            config.log_statement("Waiting on pid '%i'" % pid)
            os.waitpid(pid, 0) 
    
    return


def estimate_mle_worker( gene_ids, data ):
    while True:
        config.log_statement("Retrieving gene from queue")
        gene_id = gene_ids.get()
        if gene_id == 'FINISHED': 
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
            
            try: 
                f_mat = data.get_design_matrix(gene_id)
            except NoDesignMatrixError:
                if config.DEBUG_VERBOSE:
                    config.log_statement("No design matrix for '%s'" % gene_id, 
                                         log=True)
                continue
            num_reads_in_bams = data.get_num_reads_in_bams()

            expected_array, observed_array = f_mat.expected_and_observed(
                num_reads_in_bams)
            if (expected_array, observed_array) == (None, None): 
                continue
            mle = frequency_estimation.estimate_transcript_frequencies( 
                observed_array, expected_array)
        except Exception, inst:
            error_msg = "%i: Skipping %s (%s:%s:%i-%i): %s" % (
                os.getpid(), gene.id, 
                gene.chrm, gene.strand, gene.start, gene.stop, inst)
            config.log_statement( error_msg, log=True )
            config.log_statement( traceback.format_exc(), log=True )
            continue

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

    gene_ids = multiprocessing.Queue()
    sorted_gene_ids = sorted(data.gene_ids, 
                             key=lambda x:data.gene_ntranscripts_mapping[x],
                             reverse=True)
    # sort so that the biggest genes are processed first
    
    args = [ gene_ids, data ]
    if False and config.NTHREADS == 1:
        estimate_mle_worker(*args)
    else:
        ps = []
        for i in xrange(config.NTHREADS):
            pid = os.fork()
            if pid == 0:
                try:
                    estimate_mle_worker(*args)
                except Exception, inst:
                    config.log_statement( str(error_msg), log=True )
                    config.log_statement( traceback.format_exc(), log=True )
                finally:
                    os._exit(0)
            ps.append(pid)
        
        # populate the queue
        config.log_statement("Populating MLE queue")
        for i, gene_id in enumerate(sorted_gene_ids):
            gene_ids.put(gene_id)
        for i in xrange(config.NTHREADS):
            gene_ids.put("FINISHED")
        config.log_statement("Waiting on MLE children")
        
        for pid in ps:
            os.waitpid(pid, 0)
    
    return

def build_design_matrices_worker( gene_ids, 
                                  data, fl_dists,
                                  (rnaseq_reads, promoter_reads, polya_reads)):
    config.log_statement("Reloading read data in subprocess")
    if rnaseq_reads != None: rnaseq_reads.reload()
    if promoter_reads != None: promoter_reads.reload()
    if polya_reads != None: polya_reads.reload()
    
    while True:
        config.log_statement("Acquiring gene to process")        
        gene_id = gene_ids.get()
        if gene_id == 'FINISHED': 
            config.log_statement("")
            return
        try:
            config.log_statement("Loading gene '%s'" % gene_id)
            gene = data.get_gene(gene_id)
            config.log_statement( 
                "Finding design matrix for Gene %s(%s:%s:%i-%i) - %i transcripts"%(
                    gene.id, gene.chrm, gene.strand, 
                    gene.start, gene.stop, len(gene.transcripts) ) )
            
            f_mat = f_matrix.DesignMatrix(
                gene, fl_dists, 
                rnaseq_reads, promoter_reads, polya_reads,
                config.MAX_NUM_TRANSCRIPTS_TO_QUANTIFY)
            
            config.log_statement( "WRITING DESIGN MATRIX TO DISK %s" % gene.id )
            data.set_design_matrix(gene.id, f_mat)
            config.log_statement( "FINISHED DESIGN MATRICES %s" % gene.id )

        except f_matrix.NoObservableTranscriptsError:
            if config.DEBUG_VERBOSE:
                config.log_statement(
                    "No observable transcripts for '%s'" % gene_id, log=True)
            continue
        except Exception, inst:
            error_msg = "%i: Skipping %s: %s" % (
                os.getpid(), gene_id, inst )
            config.log_statement( 
                error_msg + "\n" + traceback.format_exc(), log=True )

def build_design_matrices( data, fl_dists,
                           (rnaseq_reads, promoter_reads, polya_reads)):    
    gene_ids = multiprocessing.Queue()
    config.log_statement( "Populating build design matrices queue" )
    sorted_gene_ids = sorted(data.gene_ids, 
                             key=lambda x:data.gene_ntranscripts_mapping[x],
                             reverse=True)
    # sort so that the biggest genes are processed first
    config.log_statement("FINISHED Populating build design matrices queue")
    
    args = [ gene_ids, data, fl_dists, 
             (rnaseq_reads, promoter_reads, polya_reads)]
    if False and config.NTHREADS == 1:
        build_design_matrices_worker(*args)
    else:
        ps = []
        for i in xrange(config.NTHREADS):
            #p = multiprocessing.Process(
            #    target=build_design_matrices_worker, args=args)
            #p.start()
            #ps.append(p)
            pid = os.fork()
            if pid == 0:
                try:
                    build_design_matrices_worker(*args)
                except Exception, inst:
                    config.log_statement( str(error_msg), log=True )
                    config.log_statement( traceback.format_exc(), log=True )
                finally:
                    os._exit(0)
            ps.append(pid)
        
        # populate the queue
        config.log_statement("Populating design matrix queue")
        for i, gene_id in enumerate(sorted_gene_ids):
            gene_ids.put(gene_id)
        for i in xrange(config.NTHREADS):
            gene_ids.put('FINISHED')
        config.log_statement("Waiting on design matrix children")
        
        #while any( p != None and p.is_alive() for p in ps ):
        #while len(ps) > 0:
        #    config.log_statement(
        #        "Waiting for design matrix children processes to finish (%i/%i genes remain)"%(
        #            gene_ids.qsize(), len(data.gene_ids)))
        #    time.sleep(1.)
        for pid in ps:
            config.log_statement("Waiting on pid '%i'" % pid)
            os.waitpid(pid, 0)

    config.log_statement("Read counts: %s" % str(data.get_num_reads_in_bams()), 
                         log=True)
    
    return

def build_gene_lines_for_tracking_file(
        gene_id, data, num_reads_in_bams, fl_dists):
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

    lines = []
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
        lines.append("\t".join(line))
    
    return lines

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
        try: 
            lines = build_gene_lines_for_tracking_file(
                gene_id, data, num_reads_in_bams, fl_dists)
        except Exception, inst:
            config.log_statement("Skipping '%s': %s" % (gene_id, str(inst)))
            config.log_statement( traceback.format_exc(), log=True )
        else:
            ofp.write("\n".join(lines) + "\n" )

def quantify_transcript_expression(
    promoter_reads, rnaseq_reads, polya_reads,
    pickled_gene_fnames, fl_dists,
    ofname, sample_type=None, rep_id=None ):
    """Build transcripts
    """
    global SAMPLE_ID
    SAMPLE_ID=sample_type
    global REP_ID
    REP_ID = rep_id
    
    write_design_matrices=False

    if config.VERBOSE: config.log_statement( 
        "Initializing processing data" )        
    data = SharedData(pickled_gene_fnames)
    
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
    
    if config.ESTIMATE_LOWER_CONFIDENCE_BOUNDS:
        if config.VERBOSE: config.log_statement( 
            "Estimating lower confidence bounds" )
        estimate_confidence_bounds(data, 'lb')
        if config.VERBOSE: config.log_statement( 
            "FINISHED Estimating lower confidence bounds" )
    
    if config.ESTIMATE_UPPER_CONFIDENCE_BOUNDS:
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
