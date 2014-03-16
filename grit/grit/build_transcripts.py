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

import cPickle as pickle

import config

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
    num_reads_in_bam = num_reads_in_bam[1]
    num_reads_in_gene = num_reads_in_gene[1]
    
    fl_dist = fl_dists['mean']
    # account for paired end reads
    #### Fix this in the optimization stage by using real read counts and then
    #### a pseudo bin XXX BUG
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

    return gene

def build_design_matrices_worker(gene_id, 
                                 input_queue, input_queue_lock,
                                 op_lock, output,
                                 (rnaseq_reads, promoter_reads, polya_reads)):
    config.log_statement( "Finding design matrix for Gene %s" % gene_id  )
    with op_lock:
        gene = output[(gene_id, 'gene')]
        fl_dists = output[(gene.id, 'fl_dists')]
    config.log_statement( 
        "Finding design matrix for Gene %s(%s:%s:%i-%i) - %i transcripts"\
            % (gene.id, gene.chrm, gene.strand, 
               gene.start, gene.stop, len(gene.transcripts) ) )
    
    try:
        f_mat = DesignMatrix(gene, fl_dists, 
                             rnaseq_reads, promoter_reads, polya_reads,
                             config.MAX_NUM_TRANSCRIPTS)
    except ValueError, inst:
        error_msg = "%i: Skipping %s (%s:%s:%i-%i): %s" % (
            os.getpid(), gene.id, 
            gene.chrm, gene.strand, gene.start, gene.stop, inst)
        config.log_statement( error_msg )
        with input_queue_lock:
            input_queue.append(
                ('ERROR', ((gene.id, None), 
                           traceback.format_exc())))
        return
    except MemoryError, inst:
        error_msg = "%i: Skipping %s (%s:%s:%i-%i): %s" % (
            os.getpid(), gene.id, 
            gene.chrm, gene.strand, gene.start, gene.stop, inst)
        config.log_statement( error_msg )
        with input_queue_lock:
            input_queue.append(
                ('ERROR', ((gene.id, None), error_msg)))

        return

    config.log_statement( "FINISHED DESIGN MATRICES %s" % gene.id )
    config.log_statement( "" )
    return f_mat


def estimate_expression_worker( work_type, (gene_id,sample_id,trans_index),
                                     input_queue, input_queue_lock,
                                     op_lock, output, 
                                     estimate_confidence_bounds,
                                     cb_alpha=config.CB_SIG_LEVEL):
    try:
        config.log_statement( "Loading estimation data for Gene %s" % gene_id  )
        with op_lock:
            gene = output[(gene_id, 'gene')]
            f_mat = output[(gene_id, 'design_matrices')]
            fl_dists = output[(gene_id, 'fl_dists')]
            num_rnaseq_reads = output['num_rnaseq_reads']
            num_cage_reads = output['num_cage_reads']
            num_polya_reads = output['num_polya_reads']
        num_reads_in_bams = (num_cage_reads, num_rnaseq_reads, num_polya_reads)
        num_reads_in_gene = (
            f_mat.num_fp_reads, f_mat.num_rnaseq_reads, f_mat.num_tp_reads)
        
        if work_type == 'mle':
            config.log_statement(
                "Finding MLE for Gene %s(%s:%s:%i-%i) - %i transcripts" \
                    % (gene_id, gene.chrm, gene.strand, 
                       gene.start, gene.stop, len(gene.transcripts) ) )
            
            try:
                expected_array, observed_array = f_mat.expected_and_observed()
                mle = frequency_estimation.estimate_transcript_frequencies( 
                    observed_array, expected_array)
                fpkms = calc_fpkm( 
                    gene, fl_dists, mle, 
                    num_reads_in_bams, num_reads_in_gene )
                    
            except ValueError, inst:
                error_msg = "Skipping %s: %s" % ( gene_id, inst )
                config.log_statement( error_msg )
                with input_queue_lock:
                    input_queue.append(('ERROR', (
                                (gene_id, trans_index), 
                                error_msg)))
                return
            
            log_lhd = frequency_estimation.calc_lhd( 
                mle, observed_array, expected_array)
            config.log_statement( "FINISHED MLE %s\t%.2f - updating queues" % ( 
                    gene_id, log_lhd ) )
            
            with op_lock:
                output[(gene_id, 'mle')] = mle
                output[(gene_id, 'fpkm')] = fpkms
            
            if estimate_confidence_bounds:
                with op_lock:
                    output[(gene_id, 'ub')] = [None]*len(mle)
                    output[(gene_id, 'lb')] = [None]*len(mle)
                
                grouped_indices = []
                for i in xrange(expected_array.shape[1]):
                    if i%config.NUM_TRANS_IN_GRP == 0:
                        grouped_indices.append( [] )
                    grouped_indices[-1].append( i )

                with input_queue_lock:
                    for indices in grouped_indices:
                        input_queue.append( ('lb', (gene_id, None, indices)) )
                        input_queue.append( ('ub', (gene_id, None, indices)) )
            else:
                with input_queue_lock:
                    input_queue.append(('FINISHED', (gene_id, None, None)))
            config.log_statement("")

        elif work_type in ('lb', 'ub'):
            with op_lock:
                mle_estimate = output[(gene_id, 'mle')]
            expected_array, observed_array = f_mat.expected_and_observed()

            bnd_type = 'LOWER' if work_type == 'lb' else 'UPPER'

            if type(trans_index) == int:
                trans_indices = [trans_index,]
            else:
                assert isinstance( trans_index, list )
                trans_indices = trans_index

            res = []
            config.log_statement( 
                "Estimating %s confidence bound for gene %s transcript %i-%i/%i" % ( 
                    bnd_type,gene_id,trans_indices[0]+1, trans_indices[-1]+1, 
                    mle_estimate.shape[0]))
            for trans_index in trans_indices:
                if config.DEBUG_VERBOSE: config.log_statement( 
                    "Estimating %s confidence bound for gene %s transcript %i/%i" % ( 
                    bnd_type,gene_id,trans_index+1,mle_estimate.shape[0]))
                p_value, bnd = frequency_estimation.estimate_confidence_bound( 
                    observed_array, expected_array, 
                    trans_index, mle_estimate, bnd_type, cb_alpha )
                if config.DEBUG_VERBOSE: config.log_statement( 
                    "FINISHED %s BOUND %s\t%s\t%i/%i\t%.2e\t%.2e" % (
                    bnd_type, gene_id, None, 
                    trans_index+1, mle_estimate.shape[0], 
                    bnd, p_value ), do_log=True )
                res.append((trans_index, bnd))
            config.log_statement( 
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
                    ub_fpkms = calc_fpkm( gene, fl_dists, 
                                          [ ubs[i] for i in xrange(len(mle)) ], 
                                          num_reads_in_bams, num_reads_in_gene,
                                          1.0 - cb_alpha)
                    output[(gene_id, 'ubs')] = ub_fpkms
                    lb_fpkms = calc_fpkm( gene, fl_dists, 
                                          [ lbs[i] for i in xrange(len(mle)) ], 
                                          num_reads_in_bams, num_reads_in_gene,
                                          cb_alpha )
                    output[(gene_id, 'lbs')] = lb_fpkms
                    with input_queue_lock:
                        input_queue.append(('FINISHED', (gene_id, None, None)))
            
            config.log_statement("")
    
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
    config.log_statement("Initializing background writer")
    while True:
        try:
            write_type, key = finished_genes_queue.get(timeout=0.1)
            if write_type == 'FINISHED':
                break
        except Queue.Empty:
            config.log_statement( "Waiting for write queue to fill." )
            time.sleep( 1 )
            continue
        
        # write out the design matrix
        try:            
            if write_type == 'design_matrices' and write_design_matrices:
                raise NotImplemented, "This was for debugging, and hasn't been maintained. Althouhg it would be relatively easy to fix if necessary"
                # to fix this, I would add a 'write' method to the DesignMatrix class and then
                # call it from here
                if config.DEBUG_VERBOSE: 
                    config.log_statement("Writing design matrix mat to '%s'" % ofname)
                f_mat = output_dict[(key,'design_matrices')]
                ofname = "./%s_%s.mat" % ( key[0], os.path.basename(key[1]) )
                if config.DEBUG_VERBOSE: config.log_statement("Writing mat to '%s'" % ofname)
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
                config.log_statement("" % ofname)
            if write_type == 'gtf':
                config.log_statement( "Writing GENE %s to gtf" % key )

                with output_dict_lock:
                    gene = output_dict[(key, 'gene')]
                    unobservable_transcripts = []
                    if not config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS:
                        f_mat = output_dict[(key, 'design_matrices')]
                        unobservable_transcripts = sorted(
                            f_mat.filtered_transcripts)
                                            
                    mles = output_dict[(key, 'mle')] \
                        if not config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS else None
                    fpkms = output_dict[(key, 'fpkm')] \
                        if not config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS else None
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
                    if not config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS:
                        del output_dict[(key, 'mle')]
                        del output_dict[(key, 'fpkm')]
                        del output_dict[(key, 'design_matrices')]
                    if compute_confidence_bounds:
                        del output_dict[(key, 'lbs')]
                        del output_dict[(key, 'ubs')]
                
                config.log_statement( "" )
        except Exception, inst:
            config.log_statement( "FATAL ERROR" )
            config.log_statement( traceback.format_exc() )
        
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

def add_universal_processing_data((contig, strand), gene_id, output_dict, 
                                  fl_dists, fasta):
    """Add stuff we need to provide whether we havea  list of 
       already built genes or not.
    """
    output_dict[ (gene_id, 'contig') ] = contig
    output_dict[ (gene_id, 'strand') ] = strand
    
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
    config.log_statement( "Clustering elements into genes for %s:%s" % ( contig, strand ) )
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
                fl_dists, fasta )
    
    config.log_statement("")
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

        p = Pool(config.NTHREADS)
        p.apply( add_elements_for_contig_and_strand, all_args )

    output_dict['num_rnaseq_reads'] = 0
    output_dict['num_cage_reads'] = 0
    output_dict['num_polya_reads'] = 0
    
    return

def worker( input_queue, input_queue_lock,
            expression_queue, 
            output_dict_lock, output_dict,
            promoter_reads,
            rnaseq_reads,
            polya_reads,
            finished_queue,
            write_design_matrices, 
            estimate_confidence_bounds):
    # reload all the reads, to make sure that this process has its
    # own file pointer
    for reads in ( promoter_reads, rnaseq_reads, polya_reads ):
        if reads != None: reads.reload()

    # perform a maximum of 50 operations before restarting the process. This
    # is to stop python from leaking memory. We also respawn any process that
    # has been running for longer than a minute
    start_time = time.time()
    for i in xrange(50):
        # if we've been in this worker longer than a minute, then 
        # so quit and re-spawn to free any unused memory
        if time.time() - start_time > 60: return
        # get the data to process
        try:
            with input_queue_lock:
                work_type, work_data = input_queue.pop()
        except IndexError, inst:
            if len(input_queue) == 0: break
            else: continue
        
        if work_type == 'ERROR':
            ( gene_id, trans_index ), msg = work_data
            config.log_statement( str(gene_id) + "\tERROR\t" + msg, only_log=True ) 
            continue

        # unpack the work data
        gene_id, bam_fn, trans_index = work_data
        
        if work_type == 'FINISHED':
            finished_queue.put( ('gtf', gene_id) )
        elif work_type == 'gene':
            # build the gene with transcripts, and optionally call orfs
            gene = build_genes_worker(gene_id, output_dict_lock, output_dict)
            # add the gene to the output structure
            with output_dict_lock: output_dict[(gene_id, 'gene')] = gene
            # if we are only building candidate transcripts, then we are done
            if config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS:
                input_queue.append(('FINISHED', (gene_id, None, None)))
            else:
                with input_queue_lock:
                    input_queue.append( ('design_matrices', (gene.id, None, None)) )
        elif work_type == 'design_matrices':
            # otherwise, we build the design matrices
            f_mat = build_design_matrices_worker(
                gene_id, 
                input_queue, input_queue_lock,
                output_dict_lock, output_dict,
                (rnaseq_reads, promoter_reads, polya_reads) )
            if f_mat == None: continue
            with output_dict_lock: 
                output_dict[(gene_id, 'design_matrices')] = f_mat
                output_dict['num_rnaseq_reads'] += f_mat.num_rnaseq_reads
                if f_mat.num_fp_reads != None:
                    output_dict['num_cage_reads'] += f_mat.num_fp_reads
                if f_mat.num_tp_reads != None:
                    output_dict['num_polya_reads'] += f_mat.num_tp_reads
            with input_queue_lock:
                expression_queue.append( ('mle', (gene_id, None, None)) )
            if write_design_matrices:
                finished_queue.put( ('design_matrices', gene_id) )
        else:
            assert work_type in ('mle', 'ub', 'lb')
            estimate_expression_worker(
                work_type, (gene_id, bam_fn, trans_index),
                input_queue, input_queue_lock, 
                output_dict_lock, output_dict, 
                estimate_confidence_bounds)
    
    return

def spawn_and_manage_children( input_queue, expression_queue, input_queue_lock,
                               output_dict_lock, output_dict,
                               promoter_reads,
                               rnaseq_reads,
                               polya_reads,
                               finished_queue,
                               write_design_matrices, 
                               estimate_confidence_bounds):    
    ps = [None]*config.NTHREADS
    args = ( input_queue, input_queue_lock,
             expression_queue, 
             output_dict_lock, output_dict,
             promoter_reads, rnaseq_reads, polya_reads,
             finished_queue,
             write_design_matrices, 
             estimate_confidence_bounds)
    time.sleep(0.1)
    config.log_statement( "Waiting on children" )
    while True:        
        # if there is nothign in the queue, and no children are alive, we're
        # done so exit. Otherwise, wait until they finish
        if len(input_queue) == 0:
            # if any proccesses are still working, continue
            if any( p != None and p.is_alive() for p in ps ):
                time.sleep(1.)
                continue
            # we think we're completely done, but take a lock to make sure
            with input_queue_lock:
                if ( len(input_queue) == 0 and
                     len(expression_queue) == 0 and
                     all( p == None or not p.is_alive() for p in ps )):
                         break
            # if we need to migrate data from the expression queue into the
            # input queue, then do it
            with input_queue_lock:
                if ( len(input_queue) == 0 and
                     all( p == None or not p.is_alive() for p in ps )):
                    assert len(expression_queue) > 0
                    input_queue.extend(expression_queue)
                    del expression_queue[:]
        
        # sleep until we have a free process index
        while all( p != None and p.is_alive() for p in ps ):
            time.sleep(1.0)
        
        proc_i = min( i for i, p in enumerate(ps) 
                      if p == None or not p.is_alive() )
        
        if config.NTHREADS > 1:
            p = multiprocessing.Process(target=worker, args=args)
            p.start()
            if ps[proc_i] != None: ps[proc_i].join()
            ps[proc_i] = p
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
        manager = multiprocessing.Manager()
        # store commands to pass onto children
        input_queue = manager.list()
        # initialize queues to store the expression commands. Since we dont
        # know the total number of reads until the design matrices have been
        # built, we do that first and then start every commands that's been 
        # stored into expression queue
        expression_queue = manager.list()
        input_queue_lock = multiprocessing.Lock()
        # store data to be written out by the writer process
        finished_queue = manager.Queue()
        # store data that all children need to be able to access
        output_dict_lock = multiprocessing.Lock()    
        output_dict = manager.dict()
        
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
            config.log_statement( "Estimating the fragment length distribution" )
            fl_dists = build_fl_dists( elements, rnaseq_reads )
            config.log_statement( "Finished estimating the fragment length distribution" )
        else:
            fl_dists, rnaseq_reads, promoter_reads, polya_reads \
                = None, None, None, None
                
        config.log_statement( "Initializing processing data" )    
        initialize_processing_data(             
            elements, genes, fl_dists,
            rnaseq_reads, promoter_reads,
            polya_reads, fasta,
            input_queue, input_queue_lock, 
            output_dict, output_dict_lock )    
        config.log_statement( "Finished initializing processing data" )
        
        write_p = multiprocessing.Process(
            target=write_finished_data_to_disk, args=(
                output_dict, output_dict_lock, 
                finished_queue, gtf_ofp, expression_ofp,
                estimate_confidence_bounds, write_design_matrices ) )
        
        write_p.start()    
        
        spawn_and_manage_children( input_queue, expression_queue, input_queue_lock,
                                   output_dict_lock, output_dict, 
                                   promoter_reads,
                                   rnaseq_reads,
                                   polya_reads,
                                   finished_queue,
                                   write_design_matrices, 
                                   estimate_confidence_bounds)
        
        finished_queue.put( ('FINISHED', None) )
        write_p.join()
    except Exception, inst:
        config.log_statement(traceback.format_exc())
        raise
    finally:
        gtf_ofp.close()
        expression_ofp.close()
        config.log_statement("Finished building transcripts")
    
    return
