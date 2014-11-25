import os, sys
import cPickle

from collections import namedtuple

try: import grit
except ImportError: sys.path.insert(0, "/home/nboley/grit/grit/")

from grit.lib.multiprocessing_utils import ProcessSafeOPStream
from grit import config

from grit.files.reads import (
    CAGEReads, RAMPAGEReads, RNAseqReads, PolyAReads, 
    get_contigs_and_lens, fix_chrm_name_for_ucsc, clean_chr_name)
from grit.files.gtf import load_gtf
from grit.find_elements import find_all_gene_segments, get_contigs_and_lens, load_gene_bndry_bins
from grit.elements import RefElementsToInclude

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
        genes_queue, promoter_reads, rnaseq_reads, ofp,
        call_peaks_tuning_params):
    promoter_reads.reload()
    rnaseq_reads.reload()
    num_genes = genes_queue.qsize()
    while True:
        try: gene = genes_queue.get(timeout=1.0)
        except Queue.Empty: break
        
        if config.VERBOSE: config.log_statement(
                "Processing %s (%i\tremain)" % (
                    str(gene).ljust(30), genes_queue.qsize()))
        
        signal_cov, control_cov = peaks.estimate_read_and_control_cov_in_gene(
            gene, promoter_reads, 'promoter', rnaseq_reads )

        called_peaks = peaks.call_peaks(
            signal_cov, control_cov, 'promoter', **call_peaks_tuning_params)
        
        region = {'chrm': gene.chrm, 'strand':gene.strand, 
                  'start':gene.start, 'stop':gene.stop}
        shift_and_write_narrow_peak(region, called_peaks, ofp)

    return

def parse_arguments():
    allowed_assays = ['cage', 'rampage', 'rnaseq', 'polya']
    
    import argparse
    parser = argparse.ArgumentParser(
        description='Call peaks from a RAMPAGE/CAGE experiment and matching RNASeq.')

    parser.add_argument( '--reference', type=file, 
        help='GTF file containing transcripts to extract reference elements to use.')
    parser.add_argument( '--use-reference-genes', 
                         default=False, action='store_true', 
        help='Use reference gene boundaries.')
    
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
    
    parser.add_argument( '--outfname', '-o', 
                         help='Output file name. (default stdout)')
    parser.add_argument( '--gene-regions-ofname', 
                         help='Output bed file name to write gene regions to. (default: do not save gene regions)')
    
    parser.add_argument( '--ucsc', default=False, action='store_true', 
                         help='Format the contig names to work with the UCSC genome browser.')
    parser.add_argument( '--region', 
        help='Only use the specified region (contig_name:start-stop).')

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
    
    config.NTHREADS = args.threads

    config.VERBOSE = args.verbose
    if config.VERBOSE:
        assert args.outfname != None, "--outfname must be set if --verbose is set"
        from grit.lib.logging import Logger
        log_ofstream = open( "log.txt", "a" )
        log_statement = Logger(
            nthreads=config.NTHREADS,
            use_ncurses=True,
            log_ofstream=log_ofstream)
        config.log_statement = log_statement

    config.FIX_CHRM_NAMES_FOR_UCSC = args.ucsc

    call_peaks_tuning_params = {
        'alpha': 1e-2, 
        'min_noise_frac': 0.01, 
        'min_merge_size': args.min_merge_distance, 
        'min_rel_merge_size': args.min_relative_merge_distance,
        'min_rd_cnt': 5,
        'min_peak_size': 5, 
        'max_peak_size': 1000,
        'trim_fraction': args.trim_fraction,
        'max_exp_sum_fraction': args.exp_filter_fraction, 
        'max_exp_mean_cvg_fraction': args.exp_filter_fraction/10
    }
    
    if args.reference != None:
        if config.VERBOSE:
            log_statement("Loading reference genes")
        ref_genes = load_gtf(args.reference) 
        if args.use_reference_genes:
            ref_elements_to_include = RefElementsToInclude(
                True, False, False, False, False, False, exons=False )
        else:
            ref_elements_to_include = RefElementsToInclude(
                False, True, False, False, False, False, exons=True )
    else:
        assert args.use_reference_genes == False, \
            "You must set --reference to --use-reference-genes"
        ref_genes = None
        ref_elements_to_include = RefElementsToInclude(
            False, False, False, False, False, False, False )
    
    if args.region != None:
        region_data = args.region.strip().split(":")
        args.region = ( clean_chr_name(region_data[0]), 
                        [int(x) for x in region_data[1].split("-")])
    
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

    output_stream = ( open(args.outfname, "w") 
                      if args.outfname != None
                      else sys.stdout )
    
    return ( ref_genes, ref_elements_to_include, 
             promoter_reads, rnaseq_reads, 
             output_stream, 
             args.gene_regions_ofname,
             args.region,
             call_peaks_tuning_params )

def main():
    ( ref_genes, ref_elements_to_include, 
      promoter_reads, rnaseq_reads, 
      output_stream, gene_regions_ofname,
      region_to_use,
      call_peaks_tuning_params
      ) = parse_arguments()
    try:
        ofp = ProcessSafeOPStream(output_stream)
        track_name = ( 
            "peaks" if output_stream == sys.stdout else output_stream.name )
        print >> ofp, "track name=%s type=narrowPeak" % track_name

        contigs, contig_lens = get_contigs_and_lens( 
            [promoter_reads, rnaseq_reads] )
        contig_lens = dict((ctg, ctg_len)
                           for ctg, ctg_len in zip(contigs, contig_lens)
                           if ctg not in ('M',) 
                           and not ctg.startswith('Un'))

        # if we are supposed to use the annotation genes
        gene_segments = []
        if ref_elements_to_include.genes == True:
            for contig, contig_len in contig_lens.iteritems():
                for strand in '+-':
                    contig_gene_bndry_bins = load_gene_bndry_bins(
                        ref_genes, contig, strand, contig_len)
                    gene_segments.extend( contig_gene_bndry_bins )
        else:
            gene_segments = find_all_gene_segments( 
                contig_lens, 
                rnaseq_reads, promoter_reads, None,
                ref_genes, ref_elements_to_include, 
                region_to_use=region_to_use)

        if gene_regions_ofname != None:
            with open(gene_regions_ofname, "w") as genes_ofp:
                print >> genes_ofp, "track name={} type=bed".format(
                    gene_regions_ofname)
                for gene in gene_segments:
                    gene.write_elements_bed(genes_ofp)
        
        queue = multiprocessing.Queue()
        for gene in gene_segments:
            queue.put(gene)

        args = [queue, promoter_reads, rnaseq_reads, ofp, 
                call_peaks_tuning_params]
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
    try: main()
    finally: 
        try:config.log_statement.close()
        except: pass
