import sys, os
import time
import traceback

import numpy
import scipy

from pysam import Fastafile, Samfile

from itertools import izip, chain
from collections import defaultdict
import Queue

import multiprocessing
from lib.multiprocessing_utils import Pool

from files.gtf import load_gtf, Transcript, Gene
from files.reads import fix_chrm_name_for_ucsc
from transcript import cluster_exons, build_transcripts
from proteomics.ORF import find_cds_for_gene

from f_matrix import DesignMatrix
import frequency_estimation
from frag_len import load_fl_dists, FlDist, build_normal_density

import config

import cPickle as pickle
import tempfile

class ThreadSafeFile( file ):
    def __init__( *args ):
        file.__init__( *args )
        args[0].lock = multiprocessing.Lock()

    def write( self, string ):
        with self.lock:
            file.write( self, string )
            self.flush()

class SharedData(object):
    """Share data across processes.

    """
    def add_universal_processing_data(self, (contig, strand), gene_id ):
        """Add stuff we need to provide whether we havea  list of 
           already built genes or not.
        """
        self.output_dict[ (gene_id, 'contig') ] = contig
        self.output_dict[ (gene_id, 'strand') ] = strand
        
        self.output_dict[ (gene_id, 'lbs') ] = None
        self.output_dict[ (gene_id, 'ubs') ] = None
        self.output_dict[ (gene_id, 'mle') ] = None
        self.output_dict[ (gene_id, 'fpkm') ] = None
    
    def add_elements_for_contig_and_strand(
            self,
            (contig, strand), grpd_exons, 
            rnaseq_reads, promoter_reads, polya_reads ):
        gene_id_num = 1
        config.log_statement( 
            "Clustering elements into genes for %s:%s" % ( contig, strand ) )
        for ( tss_es, tes_es, internal_es, 
              se_ts, promoters, polyas ) in cluster_exons( 
                set(map(tuple, grpd_exons['tss_exon'].tolist())), 
                set(map(tuple, grpd_exons['internal_exon'].tolist())), 
                set(map(tuple, grpd_exons['tes_exon'].tolist())), 
                set(map(tuple, grpd_exons['single_exon_gene'].tolist())),
                set(map(tuple, grpd_exons['promoter'].tolist())), 
                set(map(tuple, grpd_exons['polya'].tolist())), 
                set(map(tuple, grpd_exons['intron'].tolist())), 
                strand):
            # skip genes without all of the element types
            if len(se_ts) == 0 and (
                    len(tes_es) == 0 
                    or len( tss_es ) == 0 ):
                continue

            gene_id_num += 1
            gene_id = "%s_%s_%i" % ( 
                contig, 'm' if strand == '-' else 'p', gene_id_num )

            self.add_to_input_queue('gene', gene_id)
            
            with self.output_dict_lock:
                self.output_dict[ (gene_id, 'tss_exons') ] = tss_es
                self.output_dict[ (gene_id, 'internal_exons') ] = internal_es
                self.output_dict[ (gene_id, 'tes_exons') ] = tes_es
                self.output_dict[ (gene_id, 'se_transcripts') ] = se_ts
                self.output_dict[ (gene_id, 'promoters') ] = promoters
                self.output_dict[ (gene_id, 'polyas') ] = polyas
                # XXX - BUG - FIXME
                self.output_dict[ (gene_id, 'introns') ] = grpd_exons['intron']
                
                self.add_universal_processing_data(
                    (contig, strand), gene_id )
        
        config.log_statement("")
        return    

    def init_processing_data( self, elements, genes, fasta, fl_dists,
                              rnaseq_reads, promoter_reads, polya_reads ):
        if genes != None:
            for gene in genes:
                with self.output_dict_lock:
                    self.add_universal_processing_data(
                        (gene.chrm, gene.strand), gene.id )
                self.set_gene(gene)
                assert not config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS
                self.add_build_design_matrices_to_queue(gene.id)

        else:
            args_template = [rnaseq_reads, promoter_reads, polya_reads ]
            all_args = []
            for (contig, strand), grpd_exons in elements.iteritems():
                all_args.append( [
                        self, (contig, strand), grpd_exons] + args_template)

            if config.NTHREADS == 1:
                for args in all_args:
                    SharedData.add_elements_for_contig_and_strand(*args)
            else:
                p = Pool(config.NTHREADS)
                p.apply(SharedData.add_elements_for_contig_and_strand, all_args)

        self.output_dict['fasta_fn'] = ( 
            None if fasta == None else fasta.name )
        
        self.output_dict[ 'fl_dists' ] = fl_dists

        self.output_dict['num_rnaseq_reads'] = 0
        self.output_dict['num_cage_reads'] = 0
        self.output_dict['num_polya_reads'] = 0

        return


    def get_gene(self, gene_id):
        if self._cached_gene_id == gene_id and self._cached_gene != None:
            return self._cached_gene
        
        # we don't need to lock this because genes can only be
        # set one time
        fname = self.output_dict[(gene_id, 'gene')]
        with open(fname) as fp:
            gene = pickle.load(fp)

        self._cached_gene_id = gene_id
        self._cached_gene = gene
        
        return gene
    
    def set_gene(self, gene):
        ofname = os.path.join(tempfile.mkdtemp(), gene.id + ".gene")
        # because there's no cache invalidation mechanism, we're only
        # allowed to set the gene object once. This also allows us to
        # move the load outside of the lock
        assert (gene.id, 'gene') not in self.output_dict
        with open(ofname, "w") as ofp:
            pickle.dump(gene, ofp)
        
        # put the gene into the manager
        with self.output_dict_lock: 
            self.output_dict[(gene.id, 'gene')] = ofname

    def add_build_design_matrices_to_queue(self, gene_id):
        self.add_to_input_queue('design_matrices', gene_id)
    
    def get_fl_dists(self):
        with self.output_dict_lock:
            return self.output_dict['fl_dists']
    
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
        ofname = os.path.join(tempfile.mkdtemp(), gene_id + ".fmat")
        # because there's no cache invalidation mechanism, we're only
        # allowed to set the f_mat object once. This also allows us to
        # move the load outside of the lock
        assert (gene_id, 'design_matrices') not in self.output_dict

        with open(ofname, "w") as ofp:
            pickle.dump(f_mat, ofp)
        
        with self.output_dict_lock: 
            self.output_dict[(gene_id, 'design_matrices')] = ofname
        
        return

    def add_estimate_mle_to_work_queue(self, gene_id, f_mat):
        with self.output_dict_lock: 
            if f_mat.num_rnaseq_reads != None:
                self.output_dict['num_rnaseq_reads'] += f_mat.num_rnaseq_reads
            if f_mat.num_fp_reads != None:
                self.output_dict['num_cage_reads'] += f_mat.num_fp_reads
            if f_mat.num_tp_reads != None:
                self.output_dict['num_polya_reads'] += f_mat.num_tp_reads
        
        self.add_to_expression_queue('mle', gene_id, None)
        
        return
    
    def get_num_reads_in_bams(self, gene_id):
        num_rnaseq_reads = self.output_dict['num_rnaseq_reads']
        num_cage_reads = self.output_dict['num_cage_reads']
        num_polya_reads = self.output_dict['num_polya_reads']
        return (num_cage_reads, num_rnaseq_reads, num_polya_reads)

    def get_mle(self, gene_id):
        with self.output_dict_lock: 
            return self.output_dict[(gene_id, 'mle')]
    
    def set_mle(self, gene, mle):
        assert len(mle) == len(gene.transcripts) + 1
        with self.output_dict_lock: 
            self.output_dict[(gene.id, 'mle')] = mle
    
    def add_estimate_cbs_to_work_queue(self, gene_id):
        gene = self.get_gene(gene_id)
        with self.output_dict_lock:
            self.output_dict[(gene_id, 'ub')] = [None]*len(gene.transcripts)
            self.output_dict[(gene_id, 'lb')] = [None]*len(gene.transcripts)

        grouped_indices = []
        for i in xrange(len(gene.transcripts)):
            if i%config.NUM_TRANS_IN_GRP == 0:
                grouped_indices.append( [] )
            grouped_indices[-1].append( i )

        work = []
        for indices in grouped_indices:
            self.add_to_input_queue('lb', gene_id, indices)
            self.add_to_input_queue('ub', gene_id, indices)
        
        return
    
    def set_gene_to_finished(self, gene_id):
        with self.input_queue_lock:
            self.finished_queue.put( ('gtf', gene_id) )
    
    def get_cbs(self, gene_id, cb_type):
        assert cb_type in ('ub', 'lb')
        with self.output_dict_lock: 
            return self.output_dict[(gene_id, cb_type)]

    def set_cbs(self, gene_id, cb_type, indices_and_values):
        assert cb_type in ('ub', 'lb')
        # put the gene into the manager
        with self.output_dict_lock: 
            data = self.output_dict[(gene_id, cb_type)]
            for i, v in indices_and_values:
                data[i] = v
            self.output_dict[(gene_id, cb_type)] = data
        
            # check to see if we're finished - clean this up
            ubs = self.output_dict[(gene_id, 'ub')]
            lbs = self.output_dict[(gene_id, 'lb')]

            # XXX
            #config.log_statement("CNTS: %i/%i %i/%i" % (
            #        sum(int(x != None) for x in ubs), len(ubs), 
            #        sum(int(x != None) for x in lbs), len(lbs)), log=True)
            if all(x != None for x in ubs) and all(x != None for x in lbs):
                self.set_gene_to_finished(gene_id)
        
        return
    
    def __init__(self):
        self._manager = multiprocessing.Manager()

        # initialize queues to store the expression commands. Since we dont
        # know the total number of reads until the design matrices have been
        # built, we do that first and then start every commands that's been 
        # stored into expression queue        
        self.input_queue = self._manager.dict()
        self.input_queue_lock = multiprocessing.Lock()

        # store data to be written out by the writer process
        self.expression_queue = self._manager.list()
        self.finished_queue = self._manager.Queue()

        # store data that all children need to be able to access        
        self.output_dict = self._manager.dict()
        self.output_dict_lock = multiprocessing.Lock()    
        
        # create a scratch directory to store design matrices
        self.scratch_dir = tempfile.mkdtemp()
        
        # create objects to cache gene objects, so that we dont have to do a 
        # fresh load from the shared manager
        self._cached_gene_id = None
        self._cached_gene = None
        self._cached_fmat_gene_id = None
        self._cached_fmat = None
    
    def add_to_input_queue(self, work_type, gene_id, work_data=None):
        with self.input_queue_lock:
            try:
                gene_work = self.input_queue.pop(gene_id)
                gene_work.append((work_type, work_data))
                self.input_queue[gene_id] = gene_work
            except KeyError:
                self.input_queue[gene_id] = [(work_type, work_data),]
    
    def add_to_expression_queue(self, work_type, gene_id, work_data=None):
        with self.input_queue_lock:
            self.expression_queue.append((work_type, gene_id, work_data))
    
    def migrate_expression_to_input_queue(self):
        if config.DEBUG_VERBOSE: config.log_statement( 
            "Populating input queue from expression queue" )
        for work_type, gene_id, work_data in self.expression_queue:
            try:
                gene_work = self.input_queue.pop(gene_id)
                gene_work.append((work_type, work_data))
                self.input_queue[gene_id] = gene_work
            except KeyError:
                self.input_queue[gene_id] = [(work_type, work_data),]
        del self.expression_queue[:]
    
    def get_queue_item(self, gene_id=None, 
                       raise_error_if_no_matching_gene=False):
        """Get an item from th work queue.
        
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
            fpkm = 0
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

class MaxIterError( ValueError ):
    pass

def write_gene_to_gtf( ofp, gene, mles=None, lbs=None, ubs=None, fpkms=None,
                       unobservable_transcripts=set()):
    #if mles != None:
        # we add one to accoutn for the out-of-gene transcript
    #    assert len(gene.transcripts) == len(mles)+len(unobservable_transcripts)
    n_skipped_ts = 0
    
    for index, transcript in enumerate(gene.transcripts):
        if index in unobservable_transcripts:
            n_skipped_ts += 1
            continue
        meta_data = {}
        if mles != None:
            meta_data["frac"] = ("%.2e" % mles[index-n_skipped_ts])
        if lbs != None:
            meta_data["conf_lo"] = "%.2e" % lbs[index-n_skipped_ts]
        if ubs != None:
            meta_data["conf_hi"] = "%.2e" % ubs[index-n_skipped_ts]
        if fpkms != None:
            meta_data["FPKM"] = "%.2e" % fpkms[index-n_skipped_ts]
        # choose the score to be the 1% confidence bound ratio, ie 
        # the current transcripts upper bound over all transcripts'
        # lower bound
        if lbs != None and ubs != None:
            frac = int((1000.*ubs[index-n_skipped_ts])/(1e-8+max(lbs)))
            transcript.score = max(1,min(1000,frac))
        
        if config.FIX_CHRM_NAMES_FOR_UCSC:
            transcript.chrm = fix_chrm_name_for_ucsc(transcript.chrm)
        
        ofp.write( transcript.build_gtf_lines(
                gene.id, meta_data, source="grit") + "\n" )
        ofp.flush()
    
    return

def write_gene_to_fpkm_tracking( ofp, gene, lbs=None, ubs=None, fpkms=None,
                       unobservable_transcripts=set()):
    n_skipped_ts = 0
    lines = []
    for index, transcript in enumerate(gene.transcripts):
        if index in unobservable_transcripts:
            n_skipped_ts += 1
            continue
        line = ['-']*13
        line[0] = str(transcript.id)
        line[3] = str(gene.id)
        line[6] = '%s:%i-%i' % ( gene.chrm, transcript.start, transcript.stop)
        line[7] = str( transcript.calc_length() )
        line[12] = 'OK'
        
        if fpkms != None:
            line[9] =  "%.2e" % fpkms[index-n_skipped_ts]
        if lbs != None:
            line[10] = "%.2e" % lbs[index-n_skipped_ts]
        if ubs != None:
            line[11] = "%.2e" % ubs[index-n_skipped_ts]
        
        lines.append( "\t".join(line) )
    
    ofp.write( "\n".join(lines)+"\n" )
    return

def find_matching_promoter_for_transcript(transcript, promoters):
    # find the promoter that starts at the same basepair
    # If it extends beyond the first exon, we truncate the
    # promoter at the end of the first exon
    tss_exon = transcript.exons[0] if transcript.strand == '+' \
        else transcript.exons[-1] 
    matching_promoter = None
    for promoter in promoters:
        if transcript.strand == '-' and promoter[1] == tss_exon[1]:
            matching_promoter = (max(promoter[0], tss_exon[0]), promoter[1])
        elif transcript.strand == '+' and promoter[0] == tss_exon[0]:
            matching_promoter = (promoter[0], min(promoter[1], tss_exon[1]))
    
    return matching_promoter

def find_matching_polya_region_for_transcript(transcript, polyas):
    # find the polya that ends at the same basepair
    # If it extends beyond the tes exon, we truncate the
    # polya region
    tes_exon = transcript.exons[-1] if transcript.strand == '+' \
        else transcript.exons[0] 
    matching_polya = None
    for polya in polyas:
        if transcript.strand == '+' and polya[1] == tes_exon[1]:
            matching_polya = (max(polya[0], tes_exon[0]), polya[1])
        elif transcript.strand == '-' and polya[0] == tes_exon[0]:
            matching_polya = (polya[0], min(polya[1], tes_exon[1]))
    
    return matching_polya

def build_genes_worker(gene_id, op_lock, output):
    config.log_statement("Building transcript and ORFs for Gene %s" % gene_id)
    with op_lock:
        contig = output[ (gene_id, 'contig') ]
        strand = output[ (gene_id, 'strand') ]
        tss_exons = output[ (gene_id, 'tss_exons') ]
        internal_exons = output[(gene_id, 'internal_exons')]
        tes_exons = output[ (gene_id, 'tes_exons') ]
        se_transcripts = output[ (gene_id, 'se_transcripts') ]
        introns = output[ (gene_id, 'introns') ]
        promoters = output[ (gene_id, 'promoters') ]
        polyas = output[ (gene_id, 'polyas') ]
        fasta_fn = output['fasta_fn']

    gene_min = min( min(e) for e in chain(
            tss_exons, tes_exons, se_transcripts))
    gene_max = max( max(e) for e in chain(
            tss_exons, tes_exons, se_transcripts))
    
    try:
        transcripts = []
        for i, exons in enumerate( build_transcripts( 
                tss_exons, internal_exons, tes_exons,
                se_transcripts, introns, strand ) ):
            transcript = Transcript(
                "%s_%i" % ( gene_id, i ), contig, strand, 
                exons, cds_region=None, gene_id=gene_id)
            transcript.promoter = find_matching_promoter_for_transcript(
                transcript, promoters)
            transcript.polya_region = \
               find_matching_polya_region_for_transcript(transcript, polyas)
            transcripts.append( transcript )

        gene = Gene(gene_id, contig,strand, gene_min, gene_max, transcripts)

        if fasta_fn != None:
            fasta = Fastafile( fasta_fn )
            gene.transcripts = find_cds_for_gene( 
                gene, fasta, only_longest_orf=True )

        config.log_statement("FINISHED Building transcript and ORFs for Gene %s" % gene_id)
    except Exception, inst:
        config.log_statement(
            "ERROR building transcript in %s(%s:%s:%i-%i): %s" % (
                gene_id, contig, strand, gene_min, gene_max, inst), log=True)
        return None
    
    return gene

def build_design_matrices_worker(gene, fl_dists,
                                 (rnaseq_reads, promoter_reads, polya_reads)):
    config.log_statement( 
        "Finding design matrix for Gene %s(%s:%s:%i-%i) - %i transcripts"%(
            gene.id, gene.chrm, gene.strand, 
            gene.start, gene.stop, len(gene.transcripts) ) )
    
    try:
        f_mat = DesignMatrix(gene, fl_dists, 
                             rnaseq_reads, promoter_reads, polya_reads,
                             config.MAX_NUM_TRANSCRIPTS)
    except Exception, inst:
        error_msg = "%i: Skipping %s (%s:%s:%i-%i): %s" % (
            os.getpid(), gene.id, 
            gene.chrm, gene.strand, gene.start, gene.stop, inst)
        config.log_statement( 
            error_msg + "\n" + traceback.format_exc(), log=True )
        return None
    
    config.log_statement( "FINISHED DESIGN MATRICES %s" % gene.id )
    return f_mat

def estimate_mle_worker( gene, fl_dists, f_mat, num_reads_in_bams ):
    config.log_statement(
        "Finding MLE for Gene %s(%s:%s:%i-%i) - %i transcripts" \
            % (gene.id, gene.chrm, gene.strand, 
               gene.start, gene.stop, len(gene.transcripts) ) )
    
    try:
        expected_array, observed_array = f_mat.expected_and_observed(
            num_reads_in_bams)
        mle = frequency_estimation.estimate_transcript_frequencies( 
            observed_array, expected_array)
    except Exception, inst:
        error_msg = "%i: Skipping %s (%s:%s:%i-%i): %s" % (
            os.getpid(), gene.id, 
            gene.chrm, gene.strand, gene.start, gene.stop, inst)
        config.log_statement( error_msg, log=True )
        return None

    log_lhd = frequency_estimation.calc_lhd( 
        mle, observed_array, expected_array)
    
    # add back in the missing trasncripts
    full_mle = -1*numpy.ones(len(gene.transcripts)+1, dtype=float)
    full_mle[numpy.array([-1,]+f_mat.transcript_indices().tolist())+1] = mle

    config.log_statement( "FINISHED MLE %s\t%.2f - updating queues" % ( 
            gene.id, log_lhd ) )
    
    return full_mle


def find_confidence_bounds_worker( gene, fl_dists, num_reads_in_bams,
                                   f_mat, mle_estimate, 
                                   trans_index, bound_type, 
                                   cb_alpha):
    bnd_type = 'LOWER' if bound_type == 'lb' else 'UPPER'

    if type(trans_index) == int:
        trans_indices = [trans_index,]
    else:
        assert isinstance( trans_index, list )
        trans_indices = trans_index
    
    # update the mle_estimate array to only store observable transcripts
    # add 1 to skip the out of gene bin
    try:
        observable_trans_indices = numpy.array([-1,] + f_mat.transcript_indices().tolist())+1
        mle_estimate = mle_estimate[observable_trans_indices]
    except Exception, inst:
        config.log_statement( gene.id + "\n" +
                              str(observable_trans_indices) + "\n" +
                              str(f_mat.filtered_transcripts) + "\n" + 
                              str(mle_estimate) + "\n" + 
                              traceback.format_exc(), log=True )
        return []
        
    res = []
    config.log_statement( 
        "Estimating %s confidence bound for gene %s transcript %i-%i/%i" % ( 
            bnd_type,gene.id,trans_indices[0]+1, trans_indices[-1]+1, 
            mle_estimate.shape[0]))
    
    n_skipped = sum( 1 for x in sorted(f_mat.filtered_transcripts)
                     if x < trans_indices[0])
    # XXX Make sure that this is beign counted correctly
    #n_skipped_tmp = len(set(xrange(trans_indices[0])) - \
    #    set(x-1 for x in observable_trans_indices[1:] if x-1 < trans_indices[0]))
    #config.log_statement( str([n_skipped_tmp, n_skipped, f_mat.filtered_transcripts, \
    #    observable_trans_indices, trans_indices]), log=True)
    #assert n_skipped == n_skipped_tmp
    
    for trans_index in trans_indices:
        # make sure that this isn't a filtered transcript. if it is, return
        # the broadest transcript bounds possible.
        if trans_index+1 not in observable_trans_indices:
            if bound_type == 'lb': res.append((trans_index, 0.0))
            elif bound_type == 'ub': res.append((trans_index, 1.0))
            else: assert False
            n_skipped += 1
            continue

        if config.DEBUG_VERBOSE: config.log_statement( 
            "Estimating %s confidence bound for gene %s transcript %i/%i" % ( 
            bnd_type,gene.id,trans_index+1,len(gene.transcripts)))
        try:
            p_value, bnd = frequency_estimation.estimate_confidence_bound( 
                f_mat, num_reads_in_bams,
                trans_index-n_skipped+1, mle_estimate, bnd_type, cb_alpha )
        except Exception, inst:
            p_value = 1.
            bnd = 0.0 if bound_type == 'lb' else 1.0
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
        res.append((trans_index, bnd))
    
    config.log_statement( 
        "FINISHED Estimating %s confidence bound for gene %s transcript %i-%i/%i" % ( 
            bnd_type,gene.id,trans_indices[0]+1, trans_indices[-1]+1, 
            len(gene.transcripts)))
    
    return res

# output_dict, output_dict_lock, finished_genes_queue, 
def write_finished_data_to_disk( data, 
                                 gtf_ofp, expression_ofp,
                                 compute_confidence_bounds, 
                                 write_design_matrices ):
    config.log_statement("Initializing background writer")
    while True:
        try:
            write_type, key = data.finished_queue.get(timeout=0.1)
            if write_type == 'FINISHED':
                break
        except Queue.Empty:
            config.log_statement( "Waiting for write queue to fill." )
            time.sleep( 1 )
            continue
        
        # write out the design matrix
        try:            
            if write_type == 'gtf':
                config.log_statement( "Writing GENE %s to gtf" % key )

                gene = data.get_gene(key)
                unobservable_transcripts, lbs, mles, ubs, fpkms = (
                    [], None, None, None, None)
                if not config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS:
                    f_mat = data.get_design_matrix(key)
                    unobservable_transcripts = sorted(
                        f_mat.filtered_transcripts)

                    fl_dists = data.get_fl_dists()
                    mles = data.get_mle(key)[1:]
                    lbs, ubs = None, None
                    if compute_confidence_bounds:
                        lbs = data.output_dict[(key, 'lb')]
                        ubs = data.output_dict[(key, 'ub')]
                    num_reads_in_bams = data.get_num_reads_in_bams(key)

                    # XXX
                    #config.log_statement(str(unobservable_transcripts),log=True)
                    #config.log_statement(str(lbs), log=True)
                    #config.log_statement(str(ubs), log=True)
                    #config.log_statement(str(mles), log=True)
                
                    if ubs != None:
                        ubs = calc_fpkm( gene, fl_dists, ubs,
                                         num_reads_in_bams)
                    if lbs != None:
                        lbs = calc_fpkm( gene, fl_dists, lbs,
                                         num_reads_in_bams)
                    fpkms = calc_fpkm( 
                        gene, fl_dists, mles, num_reads_in_bams)
                
                write_gene_to_gtf(gtf_ofp, gene, mles, lbs, ubs, fpkms, 
                                  unobservable_transcripts)

                if expression_ofp != None:
                    write_gene_to_fpkm_tracking( 
                        expression_ofp, gene, lbs, ubs, fpkms, 
                        unobservable_transcripts=unobservable_transcripts)

                # FIX THIS TO DELETE UNUSED DATA
                """
                with data.output_dict_lock:
                    #del data.output_dict[(key, 'gene')]
                    if not config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS:
                        del data.output_dict[(key, 'mle')]
                        del data.output_dict[(key, 'design_matrices')]
                    if compute_confidence_bounds:
                        del data.output_dict[(key, 'lbs')]
                        del data.output_dict[(key, 'ubs')]
                """
                config.log_statement( "" )
        except Exception, inst:
            config.log_statement( "Failed to write %s to GTF" % key, 
                                  log=True )
            config.log_statement( traceback.format_exc(), log=True )
        
    return

def convert_elements_to_arrays(all_elements):
    # convert into array
    all_array_elements = defaultdict( 
        lambda: defaultdict(lambda: numpy.zeros(0)) )
    for key, elements in all_elements.iteritems():
        for element_type, contig_elements in elements.iteritems():
            all_array_elements[key][element_type] \
                = numpy.array( sorted( contig_elements ) )

    return all_array_elements


def load_elements( fp ):
    all_elements = defaultdict( lambda: defaultdict(set) )
    for line in fp:
        if line.startswith( 'track' ): continue
        chrm, start, stop, element_type, score, strand = line.split()[:6]
        # subtract 1 from stop becausee beds are closed open, and we 
        # wnat everything in 0-based closed-closed
        all_elements[(chrm, strand)][element_type].add( 
            (int(start), int(stop)-1) )
    
    return convert_elements_to_arrays(all_elements)

def extract_elements_from_genes( genes ):
    all_elements = defaultdict( lambda: defaultdict(set) )
    for gene in genes:
        for key, val in gene.extract_elements().iteritems():
            all_elements[(gene.chrm, gene.strand)][key].update(val)

    
    return convert_elements_to_arrays( all_elements )

def build_fl_dists( elements, reads ):
    from frag_len import estimate_fl_dists, analyze_fl_dists, \
        estimate_normal_fl_dist_from_reads
    from transcript import iter_nonoverlapping_exons
    from files.gtf import GenomicInterval
    
    def iter_good_exons():
        num = 0
        for (chrm, strand), exons in sorted( elements.iteritems()):
            for start,stop in iter_nonoverlapping_exons(exons['internal_exon']):
                num += 1
                yield GenomicInterval(chrm, strand, start, stop)
            if config.DEBUG_VERBOSE: 
                config.log_statement("FL ESTIMATION: %s %s" % ((chrm, strand), num ))
        return
    
    good_exons = list(iter_good_exons())
    fl_dists, fragments = estimate_fl_dists( reads, good_exons )
    # if we can't estiamte it from the good exons, then use all reads to 
    # estiamte the fragment length distribution
    if len( fragments ) == 0:
        fl_dists, fragments = estimate_normal_fl_dist_from_reads( reads )
    #if False and None != fragments and  None != analyze_pdf_fname:
    #    analyze_fl_dists( fragments, analyze_pdf_fname )
    return fl_dists

def worker( data, 
            promoter_reads,
            rnaseq_reads,
            polya_reads,
            fl_dists,
            write_design_matrices, 
            estimate_confidence_bounds):
    # reload all the reads, to make sure that this process has its
    # own file pointer
    for reads in ( promoter_reads, rnaseq_reads, polya_reads ):
        if reads != None: reads.reload()

    # perform a maximum of 50 operations before restarting the process. This
    # is to stop python from leaking memory. We also respawn any process that
    # has been running for longer than a minute
    gene_id = None
    start_time = time.time()
    for i in xrange(50):
        # if we've been in this worker longer than an hour, then 
        # quit and re-spawn to free any un-reclaimed memory
        if time.time() - start_time > 3600:
            config.log_statement("")
            return
        # get the data to process
        try: work_type, gene_id, trans_indices = data.get_queue_item(
            gene_id, raise_error_if_no_matching_gene=True)
        except Queue.Empty: 
            # if less than 5 minutes has passed, get a new gene
            if time.time() - start_time < 300:
                try: work_type, gene_id, trans_indices = data.get_queue_item()
                except Queue.Empty:
                    config.log_statement("")
                    return
            else:
                config.log_statement("")
                return
        
        if work_type == 'gene':
            # build the gene with transcripts, and optionally call orfs
            gene = build_genes_worker(
                gene_id, data.output_dict_lock, data.output_dict)
            if gene == None: 
                config.log_statement("")
                continue
            # if there too many transcripts, then we are done
            if len(gene.transcripts) > config.MAX_NUM_CANDIDATE_TRANSCRIPTS:
                config.log_statement("")
                continue
            
            data.set_gene(gene)
            # if we are only building candidate transcripts, then we are done
            if config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS:
                data.set_gene_to_finished(gene.id)
            else:
                data.add_build_design_matrices_to_queue(gene.id)
        elif work_type == 'design_matrices':
            config.log_statement("Finding design matrix for Gene %s" % gene_id)
            gene = data.get_gene(gene_id)         
            
            # otherwise, we build the design matrices
            f_mat = build_design_matrices_worker(
                gene, fl_dists,
                (rnaseq_reads, promoter_reads, polya_reads) )
            # if we didn't get an f-matrix in return, then the error
            # should already have been logged so continue
            if f_mat == None: 
                config.log_statement("")
                continue
            data.set_design_matrix(gene.id, f_mat)
            data.add_estimate_mle_to_work_queue(gene.id, f_mat)
        elif work_type == 'mle':
            gene = data.get_gene(gene_id)
            f_mat = data.get_design_matrix(gene.id)
            num_reads_in_bams = data.get_num_reads_in_bams(gene.id)
            
            mle = estimate_mle_worker(gene, fl_dists, f_mat, num_reads_in_bams)
            if mle == None: 
                config.log_statement("")
                continue
            data.set_mle(gene, mle)
            if estimate_confidence_bounds:
                data.add_estimate_cbs_to_work_queue(gene.id)
            else:
                data.set_gene_to_finished(gene.id)
        else:
            assert work_type in ('ub', 'lb')

            gene = data.get_gene(gene_id)
            f_mat = data.get_design_matrix(gene.id)
            num_reads_in_bams = data.get_num_reads_in_bams(gene.id)
            mle = data.get_mle(gene.id)

            cb_estimates = find_confidence_bounds_worker( 
                gene, fl_dists, num_reads_in_bams,
                f_mat, mle, 
                trans_indices, work_type, cb_alpha=config.CB_SIG_LEVEL)

            data.set_cbs(gene_id, work_type, cb_estimates)
        
        config.log_statement("")
    
    return

def spawn_and_manage_children( data,
                               promoter_reads,
                               rnaseq_reads,
                               polya_reads,
                               write_design_matrices, 
                               estimate_confidence_bounds):    
    ps = [None]*config.NTHREADS
    args = ( data,
             promoter_reads, rnaseq_reads, polya_reads,
             data.get_fl_dists(),
             write_design_matrices, 
             estimate_confidence_bounds)
    
    while True:
        # sleep until we have a free process index
        while all( p != None and p.is_alive() for p in ps ):
            config.log_statement(
                "Waiting for free children (%i genes in queue)" 
                % len(data.input_queue) )
            time.sleep(1.0)
        
        # check to see if the queue is empty and all processes
        # are finished. If so, then we are done
        if len(data.input_queue) == 0:
            if config.DEBUG_VERBOSE: config.log_statement( 
                "Input queue is empty - acquiring input queue lock" )
            with data.input_queue_lock:
                if config.DEBUG_VERBOSE: config.log_statement( 
                    "Acquired input queue lock" )
                if len(data.input_queue) == 0 and \
                        all( p == None or not p.is_alive() for p in ps ):
                    # populate the input queue from the expression queue
                    if len(data.expression_queue) > 0:
                        data.migrate_expression_to_input_queue()
                        continue
                    else:
                        return
            if config.DEBUG_VERBOSE: config.log_statement( 
                "Waiting for children to finish" )
            time.sleep(1.)
            continue
        
        # start a worker in this slot
        if config.NTHREADS > 1:
            config.log_statement( "Spawning new worker child" )
            while True:
                # find the empty process slots
                proc_is = [ i for i, p in enumerate(ps) 
                            if p == None or not p.is_alive() ]
                # if we can't find new empty processes or the input queue is 
                # empty, then continue
                if len(proc_is) == 0 or len(data.input_queue) == 0: break
                
                for proc_i in proc_is:
                    p = multiprocessing.Process(target=worker, args=args)
                    p.start()
                    if ps[proc_i] != None: ps[proc_i].join()
                    ps[proc_i] = p
        
        # if we're running in single thread mode, there is no need
        # to spawn a process
        else:
            worker(*args)
    
    config.log_statement( "" )
    return

def build_and_quantify_transcripts(
    promoter_reads, rnaseq_reads, polya_reads,
    exons_bed_fp, transcripts_gtf_fp, 
    ofprefix, fasta=None,
    estimate_confidence_bounds=True ):
    """Build transcripts
    """
    write_design_matrices=False

    # make sure that we're starting from the start of the 
    # elements files
    if exons_bed_fp != None: exons_bed_fp.seek(0)
    if transcripts_gtf_fp != None: transcripts_gtf_fp.seek(0)
    
    gtf_ofp = ThreadSafeFile("%s.gtf" % ofprefix, "w")
    gtf_ofp.write("track name=%s useScore=1\n" % ofprefix)
    
    expression_ofp = ThreadSafeFile("%s.expression.csv" % ofprefix, "w")
    expression_ofp.write("\t".join(
            ["tracking_id", "class_code", "nearest_ref_id", "gene_id", 
             "gene_short_name", "tss_id", "locus", "length", "coverage", 
             "FPKM", "FPKM_conf_lo","FPKM_conf_hi","FPKM_status"])+"\n")
    
    try:
        # initialize a process safe data strcture to pass commands between
        # processes
        data = SharedData()
        
        elements, genes = None, None
        if exons_bed_fp != None:
            config.log_statement( "Loading %s" % exons_bed_fp.name )
            elements = load_elements( exons_bed_fp )
            config.log_statement( "Finished Loading %s" % exons_bed_fp.name )
        else:
            assert transcripts_gtf_fp != None
            config.log_statement( "Loading %s" % transcripts_gtf_fp.name )
            genes = load_gtf( transcripts_gtf_fp )
            elements = extract_elements_from_genes(genes)
            config.log_statement( "Finished Loading %s" % transcripts_gtf_fp.name )
        
        if not config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS:
            # estimate the fragment length distribution
            config.log_statement("Estimating the fragment length distribution")
            fl_dists = build_fl_dists( elements, rnaseq_reads )
            config.log_statement("Finished estimating the fragment length distribution" )
        else:
            fl_dists, rnaseq_reads, promoter_reads, polya_reads \
                = None, None, None, None
                
        config.log_statement( "Initializing processing data" )    
        data.init_processing_data(             
            elements, genes, 
            fasta, fl_dists,
            rnaseq_reads, promoter_reads, polya_reads )
        config.log_statement( "Finished initializing processing data" )
        
        write_p = multiprocessing.Process(
            target=write_finished_data_to_disk, args=(
                data, gtf_ofp, expression_ofp,
                estimate_confidence_bounds, write_design_matrices ) )
        
        write_p.start()    
        
        spawn_and_manage_children( data,
                                   
                                   promoter_reads,
                                   rnaseq_reads,
                                   polya_reads,
                                   
                                   write_design_matrices, 
                                   estimate_confidence_bounds)
        
                                                              
        data.finished_queue.put( ('FINISHED', None) )
        write_p.join()
    except Exception, inst:
        config.log_statement(traceback.format_exc(), log=True, display=False)
        raise
    finally:
        gtf_ofp.close()
        expression_ofp.close()
        config.log_statement("Finished building transcripts")
    
    return
