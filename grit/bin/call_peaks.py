import os, sys
import cPickle

from collections import namedtuple

try: import grit
except ImportError: sys.path.insert(0, "/home/nboley/grit/grit/")

from grit.lib.multiprocessing_utils import ProcessSafeOPStream
from grit import config

from grit.files.reads import (
    CAGEReads, RAMPAGEReads, RNAseqReads, PolyAReads, 
    get_contigs_and_lens, fix_chrm_name_for_ucsc)
from grit.files.gtf import load_gtf
from grit.find_elements import find_all_gene_segments, get_contigs_and_lens

from grit import peaks

import multiprocessing
import Queue

def shift_and_write_narrow_peak(region, peaks, ofp):
    """
    track name=CAGE.pan type=narrowPeak
    chr4    89932   89933   .    0    +    2    -1    -1    -1
   """
    if config.FIX_CHRM_NAMES_FOR_UCSC: 
        chrm = fix_chrm_name_for_ucsc(region['chrm'])
    else:
        chrm = region['chrm']
    for start, stop, value in peaks:
        ofp.write( "\t".join(
                 ( chrm, str(region['start'] + start), 
                   str(region['start'] + stop + 1), 
                   ".", "1000", region['strand'], 
                   "%e" % value, 
                   "-1", "-1", "-1")) + "\n")
    return

def process_genes(
        genes_queue, promoter_reads, rnaseq_reads, ofp):
    promoter_reads.reload()
    rnaseq_reads.reload()
    num_genes = genes_queue.qsize()
    while True:
        try: gene = genes_queue.get(timeout=1.0)
        except Queue.Empty: break
        
        if config.VERBOSE: config.log_statement(
                "Processing %s (%i\tremain)" % (
                    str(gene).ljust(30), genes_queue.qsize()))
        region_tuple = ( gene.chrm, gene.strand, 
                         max(0, gene.start-1000), gene.stop+1000)
        region = dict(zip(('chrm', 'strand', 'start', 'stop'), 
                          region_tuple))
        called_peaks = peaks.estimate_read_cov_and_call_peaks(
            region, promoter_reads, 'promoter', rnaseq_reads)
        shift_and_write_narrow_peak(region, called_peaks, ofp)

    return

def parse_arguments():
    allowed_assays = ['cage', 'rampage', 'rnaseq', 'polya']
    
    import argparse
    parser = argparse.ArgumentParser(
        description='Call peaks from a RAMPAGE/CAGE experiment and matching RNASeq.')

    parser.add_argument( '--reference', type=file, required=True,
        help='GTF file containing genes to extract gene boundaries from.')
    
    parser.add_argument( '--rnaseq-reads', type=argparse.FileType('rb'), 
        help='BAM file containing mapped RNAseq reads.')
    parser.add_argument( '--rnaseq-read-type', 
                         choices=["forward", "backward", "auto"],
                         default='auto',
        help="If 'forward' then the first RNAseq read in a pair that maps to the genome without being reverse complemented is assumed to be on the correct strand. default: auto")
    parser.add_argument( '--num-mapped-rnaseq-reads', type=int,
        help="The total number of mapped rnaseq reads ( needed to calculate the FPKM ). This only needs to be set if it isn't found by a call to samtools idxstats." )
    
    parser.add_argument( '--cage-reads', type=argparse.FileType('rb'),
        help='BAM file containing mapped cage reads.')
    parser.add_argument( '--cage-read-type', 
                         choices=["forward", "backward", "auto"],
                         default='auto',
        help="If 'forward' then the reads that maps to the genome without being reverse complemented are assumed to be on the '+'. default: auto")

    parser.add_argument( '--rampage-reads', type=argparse.FileType('rb'),
        help='BAM file containing mapped rampage reads.')
    parser.add_argument( '--rampage-read-type', 
                         choices=["forward", "backward", "auto"],
                         default='auto',
        help="If 'forward' then the first read in a pair that maps to the genome without being reverse complemented are assumed to be on the '+' strand. default: auto")
    
    parser.add_argument( '--out-fname', '-o', 
                         help='Output file name. (default stdout)')
    
    parser.add_argument( '--ucsc', default=False, action='store_true', 
                         help='Format the contig names to work with the UCSC genome browser.')

    parser.add_argument( '--min-merge-distance', default=50, type=int,
                         help='The distance in basepairs under whihc peaks will be merged .')
    parser.add_argument( '--min-relative-merge-distance', default=0.5, type=float,
                         help='The distance as a fraction of a peak size under which peaks will be merged .')
    parser.add_argument( '--trim-fraction', default=0.01, type=float,
                         help='The fraction of reads that will be trimmed from merged reads.')
    parser.add_argument( '--exp-filter-fraction', default=0.10, type=float,
                         help='Peaks with a relative expression fraction under this amount will be filtered.')
        
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', 
                         help='Whether or not to print status information.')
    parser.add_argument( '--threads', '-t', default=1, type=int,
                         help='The number of threads to run.')

        
    args = parser.parse_args()
    config.VERBOSE = args.verbose
    config.FIX_CHRM_NAMES_FOR_UCSC = args.ucsc
    config.NTHREADS = args.threads
    
    peaks.MIN_MERGE_SIZE = args.min_merge_distance
    peaks.MIN_REL_MERGE_SIZE = args.min_relative_merge_distance
    
    peaks.TRIM_FRACTION = args.trim_fraction
    peaks.MAX_EXP_FRACTION = args.exp_filter_fraction
    ref_genes = load_gtf(args.reference)

    RefElementsToInclude = namedtuple(
        'RefElementsToInclude', 
        ['genes', 'junctions', 'TSS', 'TES', 'promoters', 'polya_sites'])
    ref_elements_to_include = RefElementsToInclude(
        True, False, False, False, False, False )
    
    if args.cage_reads != None:
        assert args.rampage_reads == None, "Can not use RAMPAGE and CAGE reads"
        if config.VERBOSE: config.log_statement( "Loading %s" % args.cage_reads.name )
        rev_reads = {'forward':False, 'backward':True, 'auto': None}[
            args.cage_read_type]
        promoter_reads = CAGEReads(args.cage_reads.name, "rb").init(
            reverse_read_strand=rev_reads, ref_genes=ref_genes)
    elif args.rampage_reads != None:
        assert args.cage_reads == None, "Can not use RAMPAGE and CAGE reads"
        if config.VERBOSE: 
            config.log_statement( "Loading %s" % args.rampage_reads.name )
        rev_reads = {'forward':False, 'backward':True, 'auto': None}[
            args.rampage_read_type]
        promoter_reads = RAMPAGEReads(args.rampage_reads.name, "rb").init(
            reverse_read_strand=rev_reads, ref_genes=ref_genes)
    else:
        assert False, "RAMPAGE or CAGE reads must be set"

    rev_reads = {'forward':False, 'backward':True, 'auto': None}[
        args.rnaseq_read_type]
    rnaseq_reads = RNAseqReads(args.rnaseq_reads.name, "rb").init(
        reverse_read_strand=rev_reads, ref_genes=ref_genes)

    output_stream = ( open(args.out_fname, "w") 
                      if args.out_fname != None
                      else sys.stdout )
    
    return ( ref_genes, ref_elements_to_include, 
             promoter_reads, rnaseq_reads, 
             output_stream )

def main():
    ( ref_genes, ref_elements_to_include, promoter_reads, rnaseq_reads, output_stream 
      ) = parse_arguments()
    try:
        ofp = ProcessSafeOPStream(output_stream)
        track_name = ( 
            "peaks" if output_stream == sys.stdout else output_stream.name )
        print >> ofp, "track name=%s type=narrowPeak" % track_name

        contigs, contig_lens = get_contigs_and_lens( 
            [promoter_reads, rnaseq_reads] )
        contig_lens = dict(zip(contigs, contig_lens))
        
        gene_segments = find_all_gene_segments( 
            contig_lens, 
            rnaseq_reads, promoter_reads, None,
            ref_genes, ref_elements_to_include, 
            region_to_use=("20", (0, 62000000)) )

        queue = multiprocessing.Queue()
        for gene in gene_segments:
            queue.put(gene)

        args = [queue, promoter_reads, rnaseq_reads, ofp]
        if config.NTHREADS == 1:
            process_genes(*args)
        else:
            ps = []
            for i in xrange(config.NTHREADS):
                p = multiprocessing.Process(target=process_genes, args=args)
                p.daemon=True
                p.start()
                ps.append(p)
            for p in ps: p.join()
    finally:
        if output_stream != sys.stdout:
            ofp.close()
    
    return

if __name__ == '__main__':
    main()
