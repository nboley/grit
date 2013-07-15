import sys, os
import numpy

import time

from pysam import Fastafile, Samfile

from itertools import izip, chain
from collections import defaultdict
import Queue

import multiprocessing

from files.gtf import load_gtf, Transcript, Gene
from files.reads import RNAseqReads, CAGEReads, RAMPAGEReads
from transcript import cluster_exons, build_transcripts
from proteomics.ORF_finder import find_cds_for_gene

from f_matrix import build_design_matrices
import frequency_estimation
from frag_len import load_fl_dists, FlDist, build_normal_density

MAX_NUM_TRANSCRIPTS = 50000

from lib.logging import Logger
# log statement is set in the main init, and is a global
# function which facilitates smart, ncurses based logging
log_statement = None

log_fp = sys.stderr
num_threads = 1

DEBUG = True
DEBUG_VERBOSE = False

def log(text):
    if VERBOSE: log_fp.write(  text + "\n" )
    return

class ThreadSafeFile( file ):
    def __init__( *args ):
        file.__init__( *args )
        args[0].lock = multiprocessing.Lock()

    def write( self, string ):
        self.lock.acquire()
        file.write( self, string )
        self.flush()
        self.lock.release()

def calc_fpkm( gene, fl_dist, freqs, num_reads_in_bam, num_reads_in_gene ):
    fpkms = []
    for t, freq in izip( gene.transcripts, freqs ):
        num_reads_in_t = num_reads_in_gene*freq
        t_len = sum( e[1] - e[0] + 1 for e in t.exons )
        fpk = num_reads_in_t/(t_len/1000.)
        fpkm = fpk/(num_reads_in_bam/1000000.)
        fpkms.append( fpkm )
    return fpkms

class MaxIterError( ValueError ):
    pass

def write_gene_to_gtf( ofp, gene, mles=None, lbs=None, ubs=None, fpkms=None,
                       unobservable_transcripts=set()):
    if mles != None:
        assert len(gene.transcripts) == len(mles) + len(unobservable_transcripts)
    n_skipped_ts = 0
    for index, transcript in enumerate(gene.transcripts):
        if index in unobservable_transcripts:
            n_skipped_ts += 1
            continue
        meta_data = {}
        if mles != None:
            transcript.score = max(1,int((1000.*mles[index-n_skipped_ts])/max(mles)))
            meta_data["frac"] = ("%.2e" % mles[index-n_skipped_ts])
        if lbs != None:
            meta_data["conf_lo"] = "%.2e" % lbs[index-n_skipped_ts]
        if ubs != None:
            meta_data["conf_hi"] = "%.2e" % ubs[index-n_skipped_ts]
        if fpkms != None:
            meta_data["FPKM"] = "%.2e" % fpkms[index-n_skipped_ts]
        
        ofp.write( transcript.build_gtf_lines(
                gene.id, meta_data, source="grit") + "\n" )
    
    return

def estimate_gene_expression_worker( work_type, (gene_id,sample_id,trans_index),
                                     input_queue, input_queue_lock,
                                     op_lock, output, 
                                     estimate_confidence_bounds ):
    
    if work_type == 'gene':
        op_lock.acquire()
        contig = output[ (gene_id, 'contig') ]
        strand = output[ (gene_id, 'strand') ]
        tss_exons = output[ (gene_id, 'tss_exons') ]
        internal_exons = output[(gene_id, 'internal_exons')]
        tes_exons = output[ (gene_id, 'tes_exons') ]
        se_transcripts = output[ (gene_id, 'se_transcripts') ]
        introns = output[ (gene_id, 'introns') ]
        fasta_fn = output[ (gene_id, 'fasta_fn') ]
        op_lock.release()
                
        transcripts = []
        for i, exons in enumerate( build_transcripts( 
                tss_exons, internal_exons, tes_exons,
                se_transcripts, introns, strand, MAX_NUM_TRANSCRIPTS ) ):
            transcripts.append( Transcript(
                    "%s_%i" % ( gene_id, i ), contig, strand, 
                    exons, cds_region=None, gene_id=gene_id) )
        
        gene_min = min( min(e) for e in chain(
                tss_exons, tes_exons, se_transcripts))
        gene_max = max( max(e) for e in chain(
                tss_exons, tes_exons, se_transcripts))
        gene = Gene( gene_id, contig, strand, gene_min, gene_max, transcripts )

        if fasta_fn != None:
            fasta = Fastafile( fasta_fn )
            gene.transcripts = find_cds_for_gene( gene, fasta, only_longest_orf=True )
        
        op_lock.acquire()
        output[(gene_id, 'gene')] = gene
        op_lock.release()
        
        # only try and build the design matrix if we were able to build full 
        # length transcripts
        if len( gene.transcripts ) > 0:
            input_queue_lock.acquire()
            input_queue.append( ('design_matrices', (gene_id, None, None)) )
            input_queue_lock.release()
        
    elif work_type == 'design_matrices':
        log_statement( "Finding design matrix for Gene %s" % gene_id  )
        op_lock.acquire()
        gene = output[(gene_id, 'gene')]
        fl_dists = output[(gene_id, 'fl_dists')]
        promoter_reads_init_data = output[(gene_id, 'promoter_reads')]
        rnaseq_reads_init_data = output[(gene_id, 'rnaseq_reads')]
        op_lock.release()
        rnaseq_reads = [ RNAseqReads(fname).init(**kwargs) 
                         for fname, kwargs in rnaseq_reads_init_data ][0]
        promoter_reads = [ readsclass(fname).init(**kwargs) 
                         for readsclass, fname, kwargs in promoter_reads_init_data ]
        try:
            expected_array, observed_array, unobservable_transcripts \
                = build_design_matrices( gene, rnaseq_reads, 
                                         fl_dists, promoter_reads )
        except ValueError, inst:
            error_msg = "%i: Skipping %s: %s" % ( os.getpid(), gene_id, inst )
            log_statement( error_msg )
            if DEBUG: raise
            input_queue_lock.acquire()
            input_queue.append(
                ('ERROR', ((gene_id, rnaseq_reads.filename, trans_index), error_msg)))
            input_queue_lock.release()
            return
        except MemoryError, inst:
            error_msg =  "%i: Skipping %s: %s" % ( os.getpid(), gene_id, inst )
            log_statement( error_msg )
            if DEBUG: raise
            input_queue_lock.acquire()
            input_queue.append(
                ('ERROR', ((gene_id, rnaseq_reads.filename, trans_index), error_msg)))
            input_queue_lock.release()
            return
        
        #log_statement( "FINISHED DESIGN MATRICES %s" % gene_id )
        log_statement( "" )

        op_lock.acquire()
        try:
            output[(gene_id, 'design_matrices')] = \
                ( observed_array, expected_array, unobservable_transcripts )
        except SystemError, inst:
            op_lock.release()
            error_msg =  "SYSTEM ERROR: %i: Skipping %s: %s" % ( 
                os.getpid(), gene_id, inst )
            log_statement( error_msg )
            input_queue_lock.acquire()
            input_queue.append(
                ('ERROR', ((gene_id, rnaseq_reads.filename, trans_index), error_msg)))
            input_queue_lock.release()
            return
        else:
            op_lock.release()
        
        input_queue_lock.acquire()
        if ONLY_BUILD_CANDIDATE_TRANSCRIPTS:
            input_queue.append(('FINISHED', (gene_id, None, None)))
        else:
            input_queue.append( ('mle', (gene_id, None, None)) )
        input_queue_lock.release()
    elif work_type == 'mle':
        log_statement( "Finding MLE for Gene %s" % gene_id  )
        op_lock.acquire()
        observed_array, expected_array, unobservable_transcripts = \
            output[(gene_id, 'design_matrices')]
        gene = output[(gene_id, 'gene')]
        fl_dists = output[(gene_id, 'fl_dists')]
        promoter_reads = output[(gene_id, 'promoter_reads')]
        rnaseq_reads_init_data = output[(gene_id, 'rnaseq_reads')]
        op_lock.release()
        
        rnaseq_reads = [ RNAseqReads(fname).init(args) 
                         for fname, args in rnaseq_reads_init_data ][0]
        
        try:
            mle_estimate = frequency_estimation.estimate_transcript_frequencies( 
                observed_array, expected_array)
            num_reads_in_gene = observed_array.sum()
            num_reads_in_bam = NUMBER_OF_READS_IN_BAM
            fpkms = calc_fpkm( gene, fl_dists, mle_estimate, 
                               num_reads_in_bam, num_reads_in_gene )
        except ValueError, inst:
            error_msg = "Skipping %s: %s" % ( gene_id, inst )
            log_statement( error_msg )
            if DEBUG: raise
            input_queue_lock.acquire()
            input_queue.append(('ERROR', (
                        (gene_id, rnaseq_reads.filename, trans_index), 
                        error_msg)))
            input_queue_lock.release()
            return
        
        log_lhd = frequency_estimation.calc_lhd( 
            mle_estimate, observed_array, expected_array)
        #log_statement( "FINISHED MLE %s\t%.2f" % ( gene_id, log_lhd ) )
        log_statement("")

        op_lock.acquire()
        output[(gene_id, 'mle')] = mle_estimate
        output[(gene_id, 'fpkm')] = fpkms
        op_lock.release()

        input_queue_lock.acquire()
        if estimate_confidence_bounds:
            op_lock.acquire()
            output[(gene_id, 'ub')] = [None]*len(mle_estimate)
            output[(gene_id, 'lb')] = [None]*len(mle_estimate)
            op_lock.release()        

            for i in xrange(expected_array.shape[1]):
                if i != 1: continue
                input_queue.append( ('lb', (gene_id, None, i)) )
                input_queue.append( ('ub', (gene_id, None, i)) )
        else:
            input_queue.append(('FINISHED', (gene_id, None, None)))
        input_queue_lock.release()

    elif work_type in ('lb', 'ub'):
        op_lock.acquire()
        observed_array, expected_array, unobservable_transcripts = \
            output[(gene_id, 'design_matrices')]
        mle_estimate = output[(gene_id, 'mle')]
        op_lock.release()

        bnd_type = 'LOWER' if work_type == 'lb' else 'UPPER'

        log_statement( "Estimating %s confidence bound for gene %s transcript %i" % ( 
            bnd_type, gene_id, trans_index ) )
        p_value, bnd = frequency_estimation.estimate_confidence_bound( 
            observed_array, expected_array, 
            trans_index, mle_estimate, bnd_type )
        log_statement( "FINISHED %s BOUND %s\t%s\t%i\t%.2e\t%.2e" % ( 
            bnd_type, gene_id, None, trans_index, bnd, p_value ), do_log=True )
        
        op_lock.acquire()
        bnds = output[(gene_id, work_type+'s')]
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
            ub_fpkms = calc_fpkm( gene, fl_dists, [ ubs[i] for i in xrange(len(mle)) ], 
                                  num_reads_in_bam, num_reads_in_gene )
            output[(gene_id, 'ubs')] = ub_fpkms
            lb_fpkms = calc_fpkm( gene, fl_dists, [ lbs[i] for i in xrange(len(mle)) ], 
                                  num_reads_in_bam, num_reads_in_gene )
            output[(gene_id, 'lbs')] = lb_fpkms
            input_queue_lock.acquire()
            input_queue.append(('FINISHED', (gene_id, None, None)))
            input_queue_lock.release()
        
        op_lock.release()        
        

    return

def write_finished_data_to_disk( output_dict, output_dict_lock, 
                                 finished_genes_queue, ofp,
                                 compute_confidence_bounds=True, 
                                 write_design_matrices=False ):
    log_statement("Initializing background writer")
    while True:
        try:
            write_type, key = finished_genes_queue.get(timeout=1.0)
            if write_type == 'FINISHED':
                break
        except Queue.Empty:
            log_statement( "Waiting for write queue to fill." )
            time.sleep( 1 )
            continue
        
        # write out the design matrix
        if write_type == 'design_matrix':
            if write_design_matrices:
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
            continue
        elif write_type == 'gtf':
            output_dict_lock.acquire()
            log_statement( "Writing GENE %s to gtf" % key )
            
            gene = output_dict[(key, 'gene')]
            unobservable_transcripts = output_dict[(key, 'design_matrices')][2]
            mles = output_dict[(key, 'mle')] \
                if not ONLY_BUILD_CANDIDATE_TRANSCRIPTS else None
            fpkms = output_dict[(key, 'fpkm')] \
                if not ONLY_BUILD_CANDIDATE_TRANSCRIPTS else None
            lbs = output_dict[(key, 'lbs')] \
                if compute_confidence_bounds else None
            ubs = output_dict[(key, 'ubs')] \
                if compute_confidence_bounds else None
            write_gene_to_gtf( ofp, gene, mles, lbs, ubs, fpkms, 
                               unobservable_transcripts=unobservable_transcripts )
            
            del output_dict[(key, 'gene')]
            del output_dict[(key, 'mle')]
            del output_dict[(key, 'design_matrices')]
            del output_dict[(key, 'lbs')]
            del output_dict[(key, 'ubs')]
            
            output_dict_lock.release()
            log_statement( "" )
        
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
        # subtract 1 from stop becausee beds are closed open
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
    
    good_exons = iter_good_exons()
    fl_dists, fragments = estimate_fl_dists( reads, good_exons )
    if len( fragments ) == 0:
        x = reads.filename
        tmp_reads = Samfile( x )
        fl_dists, fragments = estimate_normal_fl_dist_from_reads( tmp_reads )
        tmp_reads.close()
    if False and None != fragments and  None != analyze_pdf_fname:
        analyze_fl_dists( fragments, analyze_pdf_fname )
    
    return fl_dists

def initialize_processing_data( elements, genes, fl_dists,
                                rnaseq_reads, promoter_reads,
                                fasta,
                                input_queue, input_queue_lock, 
                                output_dict, output_dict_lock ):
    def add_universal_data(output_dict, gene_id, contig, strand):
        """Add stuff we need to provide whether we havea  list of 
           already built genes or not.
        """
        output_dict[ (gene_id, 'contig') ] = contig
        output_dict[ (gene_id, 'strand') ] = strand

        output_dict[ (gene_id, 'rnaseq_reads') ] = \
            [(x.filename, x._init_kwargs) for x in rnaseq_reads]
        output_dict[ (gene_id, 'promoter_reads') ] = \
            [(type(x), x.filename, x._init_kwargs) for x in promoter_reads]
        output_dict[ (gene_id, 'fasta_fn') ] = None \
            if fasta == None else fasta.name

        output_dict[ (gene_id, 'fl_dists') ] = fl_dists
        output_dict[ (gene_id, 'lbs') ] = {}
        output_dict[ (gene_id, 'ubs') ] = {}
        output_dict[ (gene_id, 'mle') ] = None
        output_dict[ (gene_id, 'fpkm') ] = None
        output_dict[ (gene_id, 'design_matrices') ] = None
    
    
    gene_id = 0
    if genes != None:
        for gene in genes:
            gene_id += 1
            
            output_dict[(gene_id, 'gene')] = gene
            add_universal_data(output_dict, gene_id, gene.chrm, gene.strand)
            
            input_queue_lock.acquire()
            input_queue.append(('design_matrices', (gene_id, None, None)))
            input_queue_lock.release()            
    else:
        for (contig, strand), grpd_exons in elements.iteritems():
            for tss_es, tes_es, internal_es, se_ts in cluster_exons( 
                    set(map(tuple, grpd_exons['tss_exon'].tolist())), 
                    set(map(tuple, grpd_exons['internal_exon'].tolist())), 
                    set(map(tuple, grpd_exons['tes_exon'].tolist())), 
                    set(map(tuple, grpd_exons['single_exon_gene'].tolist())),
                    set(map(tuple, grpd_exons['intron'].tolist())), 
                    strand):
                # skip genes without all of the element types
                if len(se_ts) == 0 and (
                        len(tes_es) == 0 
                        or len( tss_es ) == 0 ):
                    continue
                
                gene_id += 1
                
                input_queue_lock.acquire()
                input_queue.append(('gene', (gene_id, None, None)))
                input_queue_lock.release()
                
                output_dict[ (gene_id, 'tss_exons') ] = tss_es
                output_dict[ (gene_id, 'internal_exons') ] = internal_es
                output_dict[ (gene_id, 'tes_exons') ] = tes_es
                output_dict[ (gene_id, 'se_transcripts') ] = se_ts
                output_dict[ (gene_id, 'introns') ] = grpd_exons['intron']

                output_dict[ (gene_id, 'gene') ] = None
                
                add_universal_data(output_dict, gene_id, contig, strand)
    
    return

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='Determine valid transcripts and estimate frequencies.')
    parser.add_argument( '--ofname', help='Output filename.', 
                         default="discovered_transcripts.gtf")
    parser.add_argument( '--elements', type=file,
        help='Bed file containing elements')

    parser.add_argument( '--transcripts', type=file,
        help='GTF file containing transcripts for which to estimate expression')
    parser.add_argument( '--rnaseq-reads', 
                         type=argparse.FileType('rb'), nargs='+',
        help='BAM files containing mapped RNAseq reads ( must be indexed ).')
    
    parser.add_argument( '--cage-reads', type=file, default=[], nargs='*', 
        help='BAM files containing mapped cage reads.')
    parser.add_argument( '--rampage-reads', type=file, default=[], nargs='*',
        help='BAM files containing mapped rampage reads.')

    parser.add_argument( '--fasta', type=file,
        help='Fasta file containing the genome sequence - if provided the ORF finder is automatically run.')
    
    parser.add_argument( '--reverse-rnaseq-strand', 
                         default=False, action="store_true",
        help='Whether to reverse the RNAseq read strand (default False).')
    parser.add_argument( '--only-build-candidate-transcripts', default=False,
        action="store_true",
        help='If set, we will output all possible transcripts without expression estimates.')
    parser.add_argument( '--estimate-confidence-bounds', '-c', default=False,
        action="store_true",
        help='Whether or not to calculate confidence bounds ( this is slow )')
    parser.add_argument( '--write-design-matrices', default=False,
        action="store_true",
        help='Write the design matrices out to a matlab-style matrix file.')

    
    parser.add_argument( '--threads', '-t', type=int , default=1,
        help='Number of threads spawn for multithreading (default=1)')

    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',
                             help='Whether or not to print status information.')
    parser.add_argument( '--debug-verbose', default=False, action='store_true',
                             help='Prints the optimization path updates.')
    parser.add_argument( '--dont-use-ncurses', default=False, action='store_true',
                             help='Whether or not to use the ncurses display.')
    
    args = parser.parse_args()
        
    if args.elements == None and args.transcripts == None:
        raise ValueError, "--elements or --transcripts must be set"

    if args.elements != None and args.transcripts != None:
        raise ValueError, "--elements and --transcripts must not both be set"

    if args.transcripts != None  and args.rnaseq_reads == None:
        raise ValueError, "--rnaseq-reads must be set if --transcripts is set"

    if args.only_build_candidate_transcripts == True and args.elements == None:
        raise ValueError, "--elements must be set if --only-build-transcripts is set"
    if args.only_build_candidate_transcripts == True and args.rnaseq_reads == None:
        raise ValueError, "--rnaseq-reads and --only-build-transcripts must not both be set"
    
    global DEBUG_VERBOSE
    DEBUG_VERBOSE = args.debug_verbose
    frequency_estimation.DEBUG_VERBOSE = DEBUG_VERBOSE
    
    global VERBOSE
    VERBOSE = ( args.verbose or DEBUG_VERBOSE )
    frequency_estimation.VERBOSE = VERBOSE
        
    global PROCESS_SEQUENTIALLY
    if args.threads == 1:
        PROCESS_SEQUENTIALLY = True

    global ONLY_BUILD_CANDIDATE_TRANSCRIPTS
    ONLY_BUILD_CANDIDATE_TRANSCRIPTS = args.only_build_candidate_transcripts
    if not ONLY_BUILD_CANDIDATE_TRANSCRIPTS and len( args.rnaseq_reads ) == 0:
        raise ValueError, "Must provide RNAseq data to estimate transcript frequencies"
    
    global num_threads
    num_threads = args.threads
    
    ofp = ThreadSafeFile( args.ofname, "w" )
    track_name = "." + os.path.basename(args.rnaseq_reads[0].name) \
        if args.rnaseq_reads != None else ""
    ofp.write( "track name=transcripts.%s useScore=1\n" % track_name )
    
    return args.elements, args.transcripts, \
        args.rnaseq_reads, args.cage_reads, args.rampage_reads, \
        ofp, args.fasta, args.reverse_rnaseq_strand, \
        args.estimate_confidence_bounds, args.write_design_matrices, \
        not args.dont_use_ncurses

def main():
    # Get file objects from command line
    exons_bed_fp, transcripts_gtf_fp, rnaseq_bams, cage_bams, rampage_bams, \
        ofp, fasta, reverse_rnaseq_strand, \
        estimate_confidence_bounds, write_design_matrices, \
        use_ncurses \
        = parse_arguments()
        
    global log_statement
    # add an extra thread for the background writer
    log_statement = Logger(num_threads+1, use_ncurses=use_ncurses )
    frequency_estimation.log_statement = log_statement
    
    manager = multiprocessing.Manager()
    input_queue = manager.list()
    input_queue_lock = manager.Lock()
    finished_queue = manager.Queue()
    output_dict_lock = manager.Lock()    
    output_dict = manager.dict()
        
    rnaseq_reads = [ RNAseqReads(fp.name).init(
                     reverse_read_strand=reverse_rnaseq_strand) 
                     for fp in rnaseq_bams ]
    for fp in rnaseq_bams: fp.close()    
    
    global NUMBER_OF_READS_IN_BAM
    NUMBER_OF_READS_IN_BAM = sum( x.mapped for x in rnaseq_reads )
    assert ONLY_BUILD_CANDIDATE_TRANSCRIPTS or len(rnaseq_reads) == 1

    cage_reads = [ CAGEReads(fp.name).init(reverse_read_strand=True) 
                   for fp in cage_bams ]    
    for fp in cage_bams: fp.close()
    rampage_reads = [ RAMPAGEReads(fp.name).init(reverse_read_strand=True) 
                      for fp in rampage_bams ]
    for fp in rampage_bams: fp.close()
    promoter_reads = [] + cage_reads + rampage_reads
    assert len(promoter_reads) <= 1
    
    log_statement( "Finished loading data files." )
        
    elements, genes = None, None
    if exons_bed_fp != None:
        elements = load_elements( exons_bed_fp )
        log_statement( "Finished Loading %s" % exons_bed_fp.name )
    else:
        assert transcripts_gtf_fp != None
        genes = load_gtf( transcripts_gtf_fp )
        elements = extract_elements_from_genes(genes)
        log_statement( "Finished Loading %s" % transcripts_gtf_fp.name )
    
    # estimate the fragment length distribution
    fl_dists = build_fl_dists( 
        elements, rnaseq_reads, log_fp.name + ".fldist.pdf" )
    log_statement( "Finished estimating the fragment length distribution" )
    
    initialize_processing_data(             
        elements, genes, fl_dists,
        rnaseq_reads, promoter_reads,
        fasta,
        input_queue, input_queue_lock, 
        output_dict, output_dict_lock )    
    log_statement( "Finished initializing processing data" )
    
    write_p = multiprocessing.Process(target=write_finished_data_to_disk, args=(
            output_dict, output_dict_lock, 
            finished_queue, ofp,
            estimate_confidence_bounds, write_design_matrices ) )
    
    write_p.start()    
    
    ps = [None]*num_threads
    time.sleep(0.1)
    log_statement( "Waiting on children" )
    while True:        
        # get the data to process
        try:
            input_queue_lock.acquire()
            work_type, key = input_queue.pop()
        except IndexError:
            if len(input_queue) == 0 and all( 
                    p == None or not p.is_alive() for p in ps ): 
                input_queue_lock.release()
                break
            input_queue_lock.release()
            # if the queue is empty but processing is still going on,
            # then just sleep
            time.sleep(1.0)
            continue
        
        input_queue_lock.release()
        
        if work_type == 'ERROR':
            ( gene_id, bam_fn, trans_index ), msg = key
            log_statement( str(gene_id) + "\tERROR\t" + msg, only_log=True ) 
            continue
        else:
            gene_id, bam_fn, trans_index = key

        if work_type == 'FINISHED':
            finished_queue.put( ('gtf', gene_id) )
            continue

        if work_type == 'mle':
            if write_design_matrices:
                finished_queue.put( ('design_matrix', gene_id) )
        
        # sleep until we have a free process index
        while True:
            if all( p != None and p.is_alive() for p in ps ):
                time.sleep(0.1)
                continue
            break
        
        proc_i = min( i for i, p in enumerate(ps) 
                      if p == None or not p.is_alive() )
        
        # find a finished process index
        args = (work_type, (gene_id, bam_fn, trans_index),
                input_queue, input_queue_lock, 
                output_dict_lock, output_dict, 
                estimate_confidence_bounds )
        if num_threads > 1:
            p = multiprocessing.Process(
                target=estimate_gene_expression_worker, args=args )
            p.start()
            if ps[proc_i] != None: ps[proc_i].join()
            ps[proc_i] = p
        else:
            estimate_gene_expression_worker(*args)
    
    finished_queue.put( ('FINISHED', None) )
    write_p.join()

    log_statement.close()
    ofp.close()
    
    return

if __name__ == "__main__":
    main()
