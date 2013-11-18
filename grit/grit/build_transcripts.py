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
from files.reads import RNAseqReads, CAGEReads, RAMPAGEReads, PolyAReads, \
    fix_chrm_name_for_ucsc
from transcript import cluster_exons, build_transcripts
from proteomics.ORF import find_cds_for_gene

from f_matrix import build_design_matrices
import frequency_estimation
from frag_len import load_fl_dists, FlDist, build_normal_density

MAX_NUM_TRANSCRIPTS = 500

from lib.logging import Logger
# log statement is set in the main init, and is a global
# function which facilitates smart, ncurses based logging
log_statement = None

log_fp = sys.stderr
num_threads = 1

DEBUG = True
DEBUG_VERBOSE = True

ALPHA = 0.025

def log(text):
    if VERBOSE: log_fp.write(  text + "\n" )
    return

class ThreadSafeFile( file ):
    def __init__( *args ):
        file.__init__( *args )
        args[0].lock = multiprocessing.Lock()

    def write( self, string ):
        with self.lock:
            file.write( self, string )
            self.flush()

def calc_fpkm( gene, fl_dists, freqs, 
               num_reads_in_bam, num_reads_in_gene, 
               bound_alpha=0.5 ):
    fl_dist = fl_dists['mean']
    # account for paired end reads
    corrected_num_reads_in_gene = int(num_reads_in_bam*scipy.stats.beta.ppf(
            bound_alpha, 
            num_reads_in_gene+1e-6, 
            num_reads_in_bam-num_reads_in_gene+1e-6))

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
        num_reads_in_t = corrected_num_reads_in_gene*freq
        effective_t_len, t_scale = calc_effective_length_and_scale_factor(t)
        if effective_t_len <= 0: 
            fpkms.append( 0 )
            continue
        assert effective_t_len > 0
        assert num_reads_in_bam > 0
        fpk = num_reads_in_t/(effective_t_len/1000.)
        fpkm = fpk/(num_reads_in_bam/1000000.)
        fpkms.append( fpkm )
    return fpkms

class MaxIterError( ValueError ):
    pass

def write_gene_to_gtf( ofp, gene, mles=None, lbs=None, ubs=None, fpkms=None,
                       unobservable_transcripts=set()):
    if mles != None:
        assert len(gene.transcripts) == len(mles)+len(unobservable_transcripts)
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
        
        if FIX_CHRM_NAMES_FOR_UCSC:
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
        line = ['-']*12
        line[0] = str(gene.id)
        line[5] = str(transcript.id)
        line[6] = '%s:%i-%i' % ( gene.chrm, transcript.start, transcript.stop)
        line[7] = str( transcript.calc_length() )
        line[11] = 'OK'
        
        if fpkms != None:
            line[8] =  "%.2e" % fpkms[index-n_skipped_ts]
        if lbs != None:
            line[9] = "%.2e" % lbs[index-n_skipped_ts]
        if ubs != None:
            line[10] = "%.2e" % ubs[index-n_skipped_ts]
        
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

def pre_filter_design_matrices(
        observed_array, expected_array, unobservable_transcripts):
    # find the indices of the various starts
    low_expression_ts = set()
    """
    # find the intervals where the multinomial change occurs
    all_mult_regions = []
    for trans_i, t in enumerate(expected_array.T):
        mult_regions = []
        for i, val in enumerate(expected_array[0,:].cumsum()):
            if abs(val - int(val + 1e-12)) < 1e-12:
                mult_regions.append(i)
        all_mult_regions.append( mult_regions )

    assert all( len(all_mult_regions[0]) == len(x) 
                for x in all_mult_regions )
    merged_mult_regions = [0,]
    for i in xrange(len(all_mult_regions[0])):
        merged_mult_regions.append(max(x[i] for x in all_mult_regions)+1)
    """
    mult_regions = [(0,1), (1,len(observed_array))]

    for start, stop in mult_regions:
        N = observed_array[start:stop].sum()
        for trans_i, t in enumerate(expected_array.T):
            t = t[start:stop]
            non_zero_expected = t.nonzero()
            cnts = observed_array[start:stop][non_zero_expected]
            N_t = float(cnts.sum())
            ps = t[non_zero_expected]
            rv = scipy.stats.binom(N_t, 0.05*(N/N_t)*ps[start:stop])
            if (rv.ppf(0.10) > cnts[start:stop]).any():
                low_expression_ts.add(trans_i)

    new_unobservable_transcripts = [] + list(unobservable_transcripts)
    for low_exp_t_i in low_expression_ts:
        new_unobservable_transcripts.append( low_exp_t_i + sum(
                x <= low_exp_t_i for x in unobservable_transcripts ))
    unobservable_transcripts = set( new_unobservable_transcripts )
    high_exp_ts = numpy.array(sorted(set(range(expected_array.shape[1]))
                                     - low_expression_ts), dtype=int)
    expected_array = expected_array[:,high_exp_ts]
    return observed_array, expected_array, unobservable_transcripts

def estimate_gene_expression_worker( work_type, (gene_id,sample_id,trans_index),
                                     input_queue, input_queue_lock,
                                     op_lock, output, 
                                     estimate_confidence_bounds,
                                     cb_alpha=ALPHA):
    try:
        if work_type == 'gene':
            log_statement("Building transcript and ORFs for Gene %s" % gene_id)
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
                fasta_fn = output[ (gene_id, 'fasta_fn') ]
            
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
            
            gene_min = min( min(e) for e in chain(
                    tss_exons, tes_exons, se_transcripts))
            gene_max = max( max(e) for e in chain(
                    tss_exons, tes_exons, se_transcripts))
            gene = Gene(gene_id, contig,strand, gene_min, gene_max, transcripts)

            if fasta_fn != None:
                fasta = Fastafile( fasta_fn )
                gene.transcripts = find_cds_for_gene( 
                    gene, fasta, only_longest_orf=True )
            
            with op_lock:
                output[(gene_id, 'gene')] = gene

            # only try and build the design matrix if we were able to build full
            # length transcripts
            if ONLY_BUILD_CANDIDATE_TRANSCRIPTS:
                input_queue.append(('FINISHED', (gene_id, None, None)))
            elif len( gene.transcripts ) > 0:
                with input_queue_lock:
                    input_queue.append( (
                            'design_matrices', (gene_id, None, None)) )
            log_statement("")

        elif work_type == 'design_matrices':
            log_statement( "Finding design matrix for Gene %s" % gene_id  )
            with op_lock:
                gene = output[(gene_id, 'gene')]
                log_statement( 
                    "Finding design matrix for Gene %s(%s:%s:%i-%i) - %i transcripts"\
                        % (gene_id, gene.chrm, gene.strand, 
                           gene.start, gene.stop, len(gene.transcripts) ) )
                fl_dists = output[(gene_id, 'fl_dists')]
                promoter_reads_init_data = output[(gene_id, 'promoter_reads')]
                rnaseq_reads_init_data = output[(gene_id, 'rnaseq_reads')]
                polya_reads_init_data = output[(gene_id, 'polya_reads')]
            
            rnaseq_reads = [ RNAseqReads(fname).init(**kwargs) 
                             for fname, kwargs in rnaseq_reads_init_data ][0]
            promoter_reads = [ readsclass(fname).init(**kwargs) 
                             for readsclass, fname, kwargs 
                               in promoter_reads_init_data ]
            polya_reads = [ readsclass(fname).init(**kwargs) 
                             for readsclass, fname, kwargs 
                               in polya_reads_init_data ]
            try:
                expected_array, observed_array, unobservable_transcripts \
                    = build_design_matrices( gene, rnaseq_reads, fl_dists, 
                                             chain(promoter_reads, polya_reads),
                                             MAX_NUM_TRANSCRIPTS)
            except ValueError, inst:
                error_msg = "%i: Skipping %s (%s:%s:%i-%i): %s" % (
                    os.getpid(), gene_id, 
                    gene.chrm, gene.strand, gene.start, gene.stop, inst)
                log_statement( error_msg )
                with input_queue_lock:
                    input_queue.append(
                        ('ERROR', ((gene_id, trans_index), 
                                   traceback.format_exc())))
                return
            except MemoryError, inst:
                error_msg = "%i: Skipping %s (%s:%s:%i-%i): %s" % (
                    os.getpid(), gene_id, 
                    gene.chrm, gene.strand, gene.start, gene.stop, inst)
                log_statement( error_msg )
                with input_queue_lock:
                    input_queue.append(
                        ('ERROR', ((gene_id, trans_index), error_msg)))
                
                return
            
            log_statement( "FINISHED DESIGN MATRICES %s" % gene_id )
            log_statement( "" )

            with op_lock:
                try:
                    output[(gene_id, 'design_matrices')] = ( 
                       observed_array, expected_array, unobservable_transcripts)
                except SystemError, inst:                
                    error_msg =  "SYSTEM ERROR: %i: Skipping %s: %s" % ( 
                        os.getpid(), gene_id, inst )
                    log_statement( error_msg )
                    with input_queue_lock:
                        input_queue.append(
                            ('ERROR', ((gene_id, trans_index), error_msg)))
                    return

            with input_queue_lock:
                input_queue.append( ('mle', (gene_id, None, None)) )
        elif work_type == 'mle':
            log_statement( "Finding MLE for Gene %s" % gene_id  )
            with op_lock:
                observed_array, expected_array, unobservable_transcripts = \
                    output[(gene_id, 'design_matrices')]
                gene = output[(gene_id, 'gene')]
                fl_dists = output[(gene_id, 'fl_dists')]
                promoter_reads = output[(gene_id, 'promoter_reads')]
                rnaseq_reads_init_data = output[(gene_id, 'rnaseq_reads')]
            
            log_statement(
                "Finding MLE for Gene %s(%s:%s:%i-%i) - %i transcripts" \
                    % (gene_id, gene.chrm, gene.strand, 
                       gene.start, gene.stop, len(gene.transcripts) ) )
            
            rnaseq_reads = [ RNAseqReads(fname).init(args) 
                             for fname, args in rnaseq_reads_init_data ][0]

            try:
                mle = frequency_estimation.estimate_transcript_frequencies( 
                    observed_array, expected_array)
                num_reads_in_gene = observed_array.sum()
                num_reads_in_bam = NUMBER_OF_READS_IN_BAM
                fpkms = calc_fpkm( gene, fl_dists, mle, 
                                   num_reads_in_bam, num_reads_in_gene )
            except ValueError, inst:
                error_msg = "Skipping %s: %s" % ( gene_id, inst )
                log_statement( error_msg )
                with input_queue_lock:
                    input_queue.append(('ERROR', (
                                (gene_id, trans_index), 
                                error_msg)))
                return
            
            log_lhd = frequency_estimation.calc_lhd( 
                mle, observed_array, expected_array)
            log_statement( "FINISHED MLE %s\t%.2f - updating queues" % ( 
                    gene_id, log_lhd ) )
            
            with op_lock:
                output[(gene_id, 'mle')] = mle
                output[(gene_id, 'fpkm')] = fpkms
            
            if estimate_confidence_bounds:
                with op_lock:
                    output[(gene_id, 'ub')] = [None]*len(mle)
                    output[(gene_id, 'lb')] = [None]*len(mle)
                
                NUM_TRANS_IN_GRP = 10
                grouped_indices = []
                for i in xrange(expected_array.shape[1]):
                    if i%NUM_TRANS_IN_GRP == 0:
                        grouped_indices.append( [] )
                    grouped_indices[-1].append( i )

                with input_queue_lock:
                    for indices in grouped_indices:
                        input_queue.append( ('lb', (gene_id, None, indices)) )
                        input_queue.append( ('ub', (gene_id, None, indices)) )
            else:
                with input_queue_lock:
                    input_queue.append(('FINISHED', (gene_id, None, None)))
            log_statement("")

        elif work_type in ('lb', 'ub'):
            with op_lock:
                observed_array, expected_array, unobservable_transcripts = \
                    output[(gene_id, 'design_matrices')]
                mle_estimate = output[(gene_id, 'mle')]
            
            bnd_type = 'LOWER' if work_type == 'lb' else 'UPPER'

            if type(trans_index) == int:
                trans_indices = [trans_index,]
            else:
                assert isinstance( trans_index, list )
                trans_indices = trans_index

            res = []
            log_statement( 
                "Estimating %s confidence bound for gene %s transcript %i-%i/%i" % ( 
                    bnd_type,gene_id,trans_indices[0]+1, trans_indices[-1]+1, 
                    mle_estimate.shape[0]))
            for trans_index in trans_indices:
                if DEBUG_VERBOSE: log_statement( 
                    "Estimating %s confidence bound for gene %s transcript %i/%i" % ( 
                    bnd_type,gene_id,trans_index+1,mle_estimate.shape[0]))
                p_value, bnd = frequency_estimation.estimate_confidence_bound( 
                    observed_array, expected_array, 
                    trans_index, mle_estimate, bnd_type, cb_alpha )
                if DEBUG_VERBOSE: log_statement( 
                    "FINISHED %s BOUND %s\t%s\t%i/%i\t%.2e\t%.2e" % (
                    bnd_type, gene_id, None, 
                    trans_index+1, mle_estimate.shape[0], 
                    bnd, p_value ), do_log=True )
                res.append((trans_index, bnd))
            log_statement( 
                "FINISHED Estimating %s confidence bound for gene %s transcript %i-%i/%i" % ( 
                    bnd_type,gene_id,trans_indices[0]+1, trans_indices[-1]+1, 
                    mle_estimate.shape[0]))
            
            with op_lock:
                bnds = output[(gene_id, work_type+'s')]
                for trans_index, bnd in res:
                    bnds[trans_index] = bnd
                output[(gene_id, work_type+'s')] = bnds
                ubs = output[(gene_id, 'ubs')]
                lbs = output[(gene_id, 'lbs')]
                mle = output[(gene_id, 'mle')]
                if len(ubs) == len(lbs) == len(mle):
                    gene = output[(gene_id, 'gene')]
                    fl_dists = output[(gene_id, 'fl_dists')]
                    num_reads_in_gene = observed_array.sum()
                    num_reads_in_bam = NUMBER_OF_READS_IN_BAM
                    ub_fpkms = calc_fpkm( gene, fl_dists, 
                                          [ ubs[i] for i in xrange(len(mle)) ], 
                                          num_reads_in_bam, num_reads_in_gene,
                                          1.0 - cb_alpha)
                    output[(gene_id, 'ubs')] = ub_fpkms
                    lb_fpkms = calc_fpkm( gene, fl_dists, 
                                          [ lbs[i] for i in xrange(len(mle)) ], 
                                          num_reads_in_bam, num_reads_in_gene,
                                          cb_alpha )
                    output[(gene_id, 'lbs')] = lb_fpkms
                    with input_queue_lock:
                        input_queue.append(('FINISHED', (gene_id, None, None)))
            
            log_statement("")
    
    except Exception, inst:
        with input_queue_lock:
            input_queue.append(
                ('ERROR', ((gene_id, trans_index), traceback.format_exc())))
    
    return

def write_finished_data_to_disk( output_dict, output_dict_lock, 
                                 finished_genes_queue, 
                                 gtf_ofp, expression_ofp,
                                 compute_confidence_bounds, 
                                 write_design_matrices ):
    log_statement("Initializing background writer")
    while True:
        try:
            write_type, key = finished_genes_queue.get(timeout=0.1)
            if write_type == 'FINISHED':
                break
        except Queue.Empty:
            log_statement( "Waiting for write queue to fill." )
            time.sleep( 1 )
            continue
        
        # write out the design matrix
        try:            
            if write_type == 'design_matrix' and write_design_matrices:
                if DEBUG_VERBOSE: 
                    log_statement("Writing design matrix mat to '%s'" % ofname)
                observed,expected,missed = output_dict[(key,'design_matrices')]
                ofname = "./%s_%s.mat" % ( key[0], os.path.basename(key[1]) )
                if DEBUG_VERBOSE: log_statement("Writing mat to '%s'" % ofname)
                savemat( ofname, {'observed': observed, 'expected': expected}, 
                         oned_as='column' )
                ofname = "./%s_%s.observed.txt" % ( 
                    key[0], os.path.basename(key[1]) )
                with open( ofname, "w" ) as ofp:
                    ofp.write("\n".join( "%e" % x for x in  observed ))
                ofname = "./%s_%s.expected.txt" % ( 
                    key[0], os.path.basename(key[1]) )
                with open( ofname, "w" ) as ofp:
                    ofp.write("\n".join( "\t".join( "%e" % y for y in x ) 
                                         for x in expected ))
                log_statement("" % ofname)
            if write_type == 'gtf':
                log_statement( "Writing GENE %s to gtf" % key )

                with output_dict_lock:
                    gene = output_dict[(key, 'gene')]
                    unobservable_transcripts = output_dict[
                        (key, 'design_matrices')][2]if (
                        not ONLY_BUILD_CANDIDATE_TRANSCRIPTS ) else []
                        
                    mles = output_dict[(key, 'mle')] \
                        if not ONLY_BUILD_CANDIDATE_TRANSCRIPTS else None
                    fpkms = output_dict[(key, 'fpkm')] \
                        if not ONLY_BUILD_CANDIDATE_TRANSCRIPTS else None
                    lbs = output_dict[(key, 'lbs')] \
                        if compute_confidence_bounds else None
                    ubs = output_dict[(key, 'ubs')] \
                        if compute_confidence_bounds else None

                write_gene_to_gtf(gtf_ofp, gene, mles, lbs, ubs, fpkms, 
                                  unobservable_transcripts=unobservable_transcripts)

                if expression_ofp != None:
                    write_gene_to_fpkm_tracking( 
                        expression_ofp, gene, lbs, ubs, fpkms, 
                        unobservable_transcripts=unobservable_transcripts)

                with output_dict_lock:
                    del output_dict[(key, 'gene')]
                    if not ONLY_BUILD_CANDIDATE_TRANSCRIPTS:
                        del output_dict[(key, 'mle')]
                        del output_dict[(key, 'fpkm')]
                        del output_dict[(key, 'design_matrices')]
                    if compute_confidence_bounds:
                        del output_dict[(key, 'lbs')]
                        del output_dict[(key, 'ubs')]
                
                log_statement( "" )
        except Exception, inst:
            log_statement( "FATAL ERROR" )
            log_statement( traceback.format_exc() )
        
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

def build_fl_dists( elements, rnaseq_reads,
                    analyze_pdf_fname=None ):
    from frag_len import estimate_fl_dists, analyze_fl_dists, \
        estimate_normal_fl_dist_from_reads
    from transcript import iter_nonoverlapping_exons
    from files.gtf import GenomicInterval
    assert len( rnaseq_reads ) == 1
    reads = rnaseq_reads[0]
    
    def iter_good_exons():
        num = 0
        for (chrm, strand), exons in sorted( 
                elements.iteritems(), 
                key=lambda x: reads.contig_len(x[0][0]) ):
            for start,stop in iter_nonoverlapping_exons(exons['internal_exon']):
                num += 1
                yield GenomicInterval(chrm, strand, start, stop)
            if DEBUG_VERBOSE: 
                log_statement("FL ESTIMATION: %s %s" % ((chrm, strand), num ))
        return
    
    good_exons = list(iter_good_exons())
    fl_dists, fragments = estimate_fl_dists( reads, good_exons )
    # if we can't estiamte it from the good exons, then use all reads to 
    # estiamte the fragment length distribution
    if len( fragments ) == 0:
        x = reads.filename
        tmp_reads = Samfile( x )
        fl_dists, fragments = estimate_normal_fl_dist_from_reads( tmp_reads )
        tmp_reads.close()
    if False and None != fragments and  None != analyze_pdf_fname:
        analyze_fl_dists( fragments, analyze_pdf_fname )
    
    return fl_dists

def add_universal_processing_data((contig, strand), gene_id, output_dict, 
                                  rnaseq_reads, promoter_reads, polya_reads, 
                                  fl_dists, fasta):
    """Add stuff we need to provide whether we havea  list of 
       already built genes or not.
    """
    output_dict[ (gene_id, 'contig') ] = contig
    output_dict[ (gene_id, 'strand') ] = strand

    output_dict[ (gene_id, 'rnaseq_reads') ] = (
        [(x.filename, x._init_kwargs) for x in rnaseq_reads] 
        if rnaseq_reads != None else None )
    output_dict[ (gene_id, 'promoter_reads') ] = (
        [(type(x), x.filename, x._init_kwargs) for x in promoter_reads]
        if promoter_reads != None else None )
    output_dict[ (gene_id, 'polya_reads') ] = (
        [(type(x), x.filename, x._init_kwargs) for x in polya_reads]
        if polya_reads != None else None )
    output_dict[ (gene_id, 'fasta_fn') ] = ( None 
        if fasta == None else fasta.name )

    output_dict[ (gene_id, 'fl_dists') ] = fl_dists
    output_dict[ (gene_id, 'lbs') ] = {}
    output_dict[ (gene_id, 'ubs') ] = {}
    output_dict[ (gene_id, 'mle') ] = None
    output_dict[ (gene_id, 'fpkm') ] = None
    output_dict[ (gene_id, 'design_matrices') ] = None


def add_elements_for_contig_and_strand((contig, strand), grpd_exons, 
                                       input_queue_lock, input_queue,
                                       output_dict_lock, output_dict,
                                       rnaseq_reads,promoter_reads,polya_reads,
                                       fl_dists, fasta ):
    gene_id_num = 1
    log_statement( "Clustering elements into genes for %s:%s" % ( contig, strand ) )
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
        
        with input_queue_lock:
            input_queue.append(('gene', (gene_id, None, None)))

        with output_dict_lock:
            output_dict[ (gene_id, 'tss_exons') ] = tss_es
            output_dict[ (gene_id, 'internal_exons') ] = internal_es
            output_dict[ (gene_id, 'tes_exons') ] = tes_es
            output_dict[ (gene_id, 'se_transcripts') ] = se_ts
            output_dict[ (gene_id, 'promoters') ] = promoters
            output_dict[ (gene_id, 'polyas') ] = polyas
            # XXX - BUG - FIXME
            output_dict[ (gene_id, 'introns') ] = grpd_exons['intron']

            output_dict[ (gene_id, 'gene') ] = None

            add_universal_processing_data(
                (contig, strand), gene_id, output_dict, 
                rnaseq_reads, promoter_reads, polya_reads, 
                fl_dists, fasta )
    
    log_statement("")
    return
    

def initialize_processing_data( elements, genes, fl_dists,
                                rnaseq_reads, promoter_reads,
                                polya_reads, fasta,
                                input_queue, input_queue_lock, 
                                output_dict, output_dict_lock ):    
    if genes != None:
        for gene in genes:
            with output_dict_lock:
                output_dict[(gene.id, 'gene')] = gene
                add_universal_processing_data(
                    (gene.chrm, gene.strand), gene.id, output_dict, 
                    rnaseq_reads, promoter_reads, polya_reads, 
                    fl_dists, fasta )
            with input_queue_lock:
                input_queue.append(('design_matrices', (gene.id, None, None)))
    else:
        args_template = [input_queue_lock, input_queue,
                         output_dict_lock, output_dict,
                         rnaseq_reads,promoter_reads,polya_reads,
                         fl_dists, fasta]
        all_args = []
        for (contig, strand), grpd_exons in elements.iteritems():
            all_args.append( [(contig, strand), grpd_exons] + args_template )

        p = Pool(num_threads)
        p.apply( add_elements_for_contig_and_strand, all_args )
    
    return

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='Determine valid transcripts and estimate frequencies.')
    parser.add_argument( '--ofname', help='Output filename.', 
                         default="discovered.transcripts.gtf")
    parser.add_argument( '--expression-ofname', 
                         help='Output filename for expression levels.', 
                         default="discovered.isoforms.fpkm_tracking")

    parser.add_argument( '--elements', type=file,
        help='Bed file containing elements')
    parser.add_argument( '--transcripts', type=file,
        help='GTF file containing transcripts for which to estimate expression')
    
    parser.add_argument( '--rnaseq-reads', 
                         type=argparse.FileType('rb'), nargs='+',
        help='BAM files containing mapped RNAseq reads ( must be indexed ).')
    parser.add_argument( '--rnaseq-read-type',
        choices=["forward", "backward"],
        help='Whether or not the first RNAseq read in a pair needs to be reversed to be on the correct strand.')
    
    parser.add_argument( '--cage-reads', type=file, default=[], nargs='*', 
        help='BAM files containing mapped cage reads.')
    parser.add_argument( '--rampage-reads', type=file, default=[], nargs='*',
        help='BAM files containing mapped rampage reads.')

    parser.add_argument( '--polya-reads', type=file, default=[], nargs='*', 
        help='BAM files containing mapped poly(A)-seq reads.')
    
    parser.add_argument( '--fasta', type=file,
        help='Fasta file containing the genome sequence - if provided the ORF finder is automatically run.')
    
    parser.add_argument( '--only-build-candidate-transcripts', default=False,
        action="store_true",
        help='If set, we will output all possible transcripts without expression estimates.')
    parser.add_argument( '--estimate-confidence-bounds', '-c', default=False,
        action="store_true",
        help='Whether or not to calculate confidence bounds ( this is slow )')
    parser.add_argument( '--write-design-matrices', default=False,
        action="store_true",
        help='Write the design matrices out to a matlab-style matrix file.')
    parser.add_argument( '--num-mapped-rnaseq-reads', type=int,
        help="The total number of mapped rnaseq reads ( needed to calculate the FPKM ). This only needs to be set if it isn't found by a call to samtools idxstats." )
    
    parser.add_argument( '--threads', '-t', type=int , default=1,
        help='Number of threads spawn for multithreading (default=1)')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',
                             help='Whether or not to print status information.')
    parser.add_argument( '--debug-verbose', default=False, action='store_true',
                             help='Prints the optimization path updates.')

    parser.add_argument( '--ucsc', default=False, action='store_true',
        help='Try to format contig names in the ucsc format (typically by prepending a chr).')    
    parser.add_argument( '--batch-mode', '-b', 
        default=False, action='store_true',
        help='Disable the ncurses frontend, and just print status messages to stderr.')
    
    args = parser.parse_args()
        
    if args.elements == None and args.transcripts == None:
        raise ValueError, "--elements or --transcripts must be set"

    if args.elements != None and args.transcripts != None:
        raise ValueError, "--elements and --transcripts must not both be set"

    if args.transcripts != None  and args.rnaseq_reads == None:
        raise ValueError, "--rnaseq-reads must be set if --transcripts is set"

    if args.only_build_candidate_transcripts == True \
            and args.elements == None:
        raise ValueError, "need --elements if --only-build-candidate-transcripts is set"
    if args.only_build_candidate_transcripts == True \
            and args.rnaseq_reads != None:
        raise ValueError, "It doesn't make sense to set --rnaseq-reads if --only-build-candidate-transcripts is not set"
    if args.only_build_candidate_transcripts == True \
            and args.estimate_confidence_bounds == True:
        raise ValueError, "--only-build-candidate-transcripts and --estimate-confidence-bounds may not both be set"
    
    reverse_rnaseq_strand = ( 
        True if args.rnaseq_read_type == 'backward' else False )
    
    global DEBUG_VERBOSE
    DEBUG_VERBOSE = args.debug_verbose
    frequency_estimation.DEBUG_VERBOSE = DEBUG_VERBOSE
    
    global VERBOSE
    VERBOSE = ( args.verbose or DEBUG_VERBOSE )
    frequency_estimation.VERBOSE = VERBOSE
        
    global PROCESS_SEQUENTIALLY
    if args.threads == 1:
        PROCESS_SEQUENTIALLY = True

    global FIX_CHRM_NAMES_FOR_UCSC
    FIX_CHRM_NAMES_FOR_UCSC = args.ucsc
    
    global NUMBER_OF_READS_IN_BAM
    NUMBER_OF_READS_IN_BAM = ( None if args.num_mapped_rnaseq_reads == None 
                               else args.num_mapped_rnaseq_reads )
    
    global ONLY_BUILD_CANDIDATE_TRANSCRIPTS
    ONLY_BUILD_CANDIDATE_TRANSCRIPTS = args.only_build_candidate_transcripts
    if not ONLY_BUILD_CANDIDATE_TRANSCRIPTS and args.rnaseq_reads == None:
        raise ValueError, "Must provide RNAseq data to estimate transcript frequencies"

    if args.rnaseq_reads != None and args.rnaseq_read_type == None:
        raise ValueError, "--rnaseq-read-type must be set if --rnaseq-reads is set"

    if args.rnaseq_reads == None and args.rnaseq_read_type != None:
        raise ValueError, "It doesn't make sense to set --rnaseq-read-type if --rnaseq-reads is not set"
    
    global num_threads
    num_threads = args.threads
    
    gtf_ofp = ThreadSafeFile( args.ofname, "w" )
    if args.ofname == None:
        track_name = ".".join( os.path.basename(x.name) 
                               for x in args.rnaseq_reads)
    else:
        track_name = args.ofname
    
    gtf_ofp.write( "track name=transcripts.%s useScore=1\n" % track_name )

    expression_ofp = ThreadSafeFile( args.expression_ofname, "w" )
    columns = [ "tracking_id", "class_code", "nearest_ref_id", "gene_id", 
                "gene_short_name", "tss_id", "locus", "length", "coverageFPKM", 
                "FPKM_conf_lo", "FPKM_conf_hi", "FPKM_status" ]

    expression_ofp.write( "\t".join(columns) + "\n" )
    
    return ( args.elements, args.transcripts, args.rnaseq_reads, 
             args.cage_reads, args.rampage_reads, args.polya_reads,
             gtf_ofp, expression_ofp, args.fasta, reverse_rnaseq_strand, 
             args.estimate_confidence_bounds, args.write_design_matrices, 
             not args.batch_mode )

def worker( input_queue, input_queue_lock,
            output_dict_lock, output_dict,
            finished_queue,
            write_design_matrices, 
            estimate_confidence_bounds):
    for i in xrange(50):
        # get the data to process
        try:
            with input_queue_lock:
                work_type, work_data = input_queue.pop()
        except IndexError, inst:
            if len(input_queue) == 0:
                break
            else:
                continue
        
        if work_type == 'ERROR':
            ( gene_id, trans_index ), msg = work_data
            log_statement( str(gene_id) + "\tERROR\t" + msg, only_log=True ) 
            continue
        else:
            gene_id, bam_fn, trans_index = work_data

        if work_type == 'FINISHED':
            finished_queue.put( ('gtf', gene_id) )
            continue

        if work_type == 'mle':
            if write_design_matrices:
                finished_queue.put( ('design_matrix', gene_id) )
        
        args = (work_type, (gene_id, bam_fn, trans_index),
                input_queue, input_queue_lock, 
                output_dict_lock, output_dict, 
                estimate_confidence_bounds )
        estimate_gene_expression_worker(*args)
    
    return

def spawn_and_manage_children( input_queue, input_queue_lock,
                               output_dict_lock, output_dict,
                               finished_queue,
                               write_design_matrices, 
                               estimate_confidence_bounds):
    ps = [None]*num_threads
    args = ( input_queue, input_queue_lock,
             output_dict_lock, output_dict,
             finished_queue,
             write_design_matrices, 
             estimate_confidence_bounds)
    time.sleep(0.1)
    log_statement( "Waiting on children" )
    while True:        
        # if there is nothign in the queue, and no children are alive, we're
        # done so exit. Otherwise, wait until they finish
        if len(input_queue) == 0:
            if all( p == None or not p.is_alive() for p in ps ): 
                break        
            else:
                time.sleep(1.)
                continue
        
        # sleep until we have a free process index
        while all( p != None and p.is_alive() for p in ps ):
            time.sleep(1.0)
        
        proc_i = min( i for i, p in enumerate(ps) 
                      if p == None or not p.is_alive() )
        
        if num_threads > 1:
            p = multiprocessing.Process(target=worker, args=args)
            p.start()
            if ps[proc_i] != None: ps[proc_i].join()
            ps[proc_i] = p
        else:
            worker(*args)
    
    log_statement( "" )
    return

def main():
    # Get file objects from command line
    (exons_bed_fp, transcripts_gtf_fp, 
     rnaseq_bams, cage_bams, rampage_bams, polya_bams,
     gtf_ofp, expression_ofp, fasta, reverse_rnaseq_strand,
     estimate_confidence_bounds, write_design_matrices, 
     use_ncurses) = parse_arguments()
    
    global log_statement
    # add an extra thread for the background writer
    log_fp = open( gtf_ofp.name + ".log", "w" )
    log_statement = Logger(num_threads+1, 
                           use_ncurses=use_ncurses, 
                           log_ofstream=log_fp )
    frequency_estimation.log_statement = log_statement
    
    try:
        manager = multiprocessing.Manager()
        input_queue = manager.list()
        input_queue_lock = multiprocessing.Lock()
        finished_queue = manager.Queue()
        output_dict_lock = multiprocessing.Lock()    
        output_dict = manager.dict()

        elements, genes = None, None
        if exons_bed_fp != None:
            log_statement( "Loading %s" % exons_bed_fp.name )
            elements = load_elements( exons_bed_fp )
            log_statement( "Finished Loading %s" % exons_bed_fp.name )
        else:
            assert transcripts_gtf_fp != None
            log_statement( "Loading %s" % transcripts_gtf_fp.name )
            genes = load_gtf( transcripts_gtf_fp )
            elements = extract_elements_from_genes(genes)
            log_statement( "Finished Loading %s" % transcripts_gtf_fp.name )
        
        if not ONLY_BUILD_CANDIDATE_TRANSCRIPTS:
            log_statement( "Loading data files." )
            rnaseq_reads = [ RNAseqReads(fp.name).init(
                             reverse_read_strand=reverse_rnaseq_strand) 
                             for fp in rnaseq_bams ]
            for fp in rnaseq_bams: fp.close()    

            global NUMBER_OF_READS_IN_BAM
            if NUMBER_OF_READS_IN_BAM == None:
                NUMBER_OF_READS_IN_BAM = 0
                NUMBER_OF_READS_IN_BAM += sum( 
                    x.mapped/2 for x in rnaseq_reads if x.reads_are_paired )
                NUMBER_OF_READS_IN_BAM += sum( 
                    x.mapped for x in rnaseq_reads if not x.reads_are_paired )
                if NUMBER_OF_READS_IN_BAM == 0:
                    raise ValueError, "Can't determine the number of reads in the RNASeq BAM (by samtools idxstats). Please set --num-mapped-rnaseq-reads"
            assert NUMBER_OF_READS_IN_BAM > 0
            
            cage_reads = [ CAGEReads(fp.name).init(reverse_read_strand=True) 
                           for fp in cage_bams ]    
            for fp in cage_bams: fp.close()
            rampage_reads = [ 
                RAMPAGEReads(fp.name).init(reverse_read_strand=True) 
                for fp in rampage_bams ]
            for fp in rampage_bams: fp.close()
            promoter_reads = [] + cage_reads + rampage_reads
            assert len(promoter_reads) <= 1    
            
            polya_reads = [
                PolyAReads(fp.name).init(
                    reverse_read_strand=True, pairs_are_opp_strand=True)
                for fp in polya_bams  ]
            assert len(polya_reads) <= 1
            log_statement( "Finished loading data files." )
            
            # estimate the fragment length distribution
            log_statement( "Estimating the fragment length distribution" )
            fl_dists = build_fl_dists( 
                elements, rnaseq_reads, log_fp.name + ".fldist.pdf" )
            log_statement( "Finished estimating the fragment length distribution" )
        else:
            fl_dists, rnaseq_reads, promoter_reads, polya_reads \
                = None, None, None, None
                
        log_statement( "Initializing processing data" )    
        initialize_processing_data(             
            elements, genes, fl_dists,
            rnaseq_reads, promoter_reads,
            polya_reads, fasta,
            input_queue, input_queue_lock, 
            output_dict, output_dict_lock )    
        log_statement( "Finished initializing processing data" )

        write_p = multiprocessing.Process(
            target=write_finished_data_to_disk, args=(
                output_dict, output_dict_lock, 
                finished_queue, gtf_ofp, expression_ofp,
                estimate_confidence_bounds, write_design_matrices ) )

        write_p.start()    

        spawn_and_manage_children( input_queue, input_queue_lock,
                                   output_dict_lock, output_dict, 
                                   finished_queue,
                                   write_design_matrices, 
                                   estimate_confidence_bounds)
        
        finished_queue.put( ('FINISHED', None) )
        write_p.join()
    except Exception, inst:
        log_statement(traceback.format_exc())
        log_statement.close()
        raise
    else:
        log_statement.close()
    finally:
        gtf_ofp.close()
        log_fp.close()
        expression_ofp.close()
    
    return

if __name__ == "__main__":
    main()
