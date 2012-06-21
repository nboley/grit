# Copyright (c) 2011-2012 Nathan Boley

#!/usr/bin/python

import sys
import os

VERBOSE = True
MIN_VERBOSE = True

# parameters for negative control thresholding
THRESHOLD_WIG = False
SCALE_WIG = False

# set parameters for distal exon finding algorithm
import distal_exons
distal_exons.EXON_FILTER_TYPE = 'region'
distal_exons.REGION_FILTER_LEN = 50
distal_exons.MIN_SCORE_RATIO = 0.01
distal_exons.CVRG_FRACTION = 0.95
distal_exons.SMOOTHING_WINDOW_SIZE = None
distal_exons.RATIO_THRESH = None

#############################################################

# import slide modules
from distal_exons import find_all_distal_exons, convert_to_genomic_intervals

sys.path.append( os.path.join( os.path.dirname(__file__), ".." ) )
import slide
slide.VERBOSE = VERBOSE
from slide import get_raw_transcripts
from gene_models import GenomicInterval

sys.path.append( os.path.join( os.path.dirname(__file__), "..", 'file_types' ) )
from wiggle import wiggle
from clustered_exons_file import parse_clustered_exons_file
from exons_file import parse_exons_file


def create_tss_gff_line( tss ):
    gff_line = [ 'chr' + tss.chr, 'elements', 'TSS', str(tss.start), \
                     str(tss.stop), '.', tss.strand, '.' ]
    return '\t'.join( gff_line )

def write_TSSs( tss_s, out_fp ):
    for tss in sorted( tss_s ):
        out_fp.write( create_tss_gff_line( tss ) + '\n' )
    out_fp.close()

def get_tss_exons_from_trans( trans_fp ):
    tss_s = set()
    
    # get the gene -> trans ->  exons structure from the gtf file
    # verification of transcripts is done from within get_raw_transcripts
    raw_gene_transcripts = get_raw_transcripts( trans_fp )
    for gene_name, transcripts in raw_gene_transcripts.iteritems():
        for trans_name, exons in transcripts.iteritems():
            if exons[0].strand == '+':
                tss = GenomicInterval( \
                    exons[0].chr, '+', exons[0].start, exons[0].stop )
            else:
                assert exons[0].strand == '-'
                tss = GenomicInterval( \
                    exons[-1].chr, '-', exons[-1].start, exons[-1].stop )
            tss_s.add( tss )
    
    return tss_s

def get_tss_wig( chrm_sizes_fp, tss_reads_fps, tss_reads_nc_wig_fps ):
    if VERBOSE:
        print '\tParsing TSS coverage wiggle(s).'
    tss_cvrg = wiggle( chrm_sizes_fp )
    for tss_wig in tss_reads_fps:
        tss_cvrg.load_data_from_fp( tss_wig )
        tss_wig.close()
    
    if THRESHOLD_WIG and tss_reads_nc_wig_fps != None:
        if VERBOSE:
            print '\tParsing negative control wiggle(s).'
        
        chrm_sizes_fp.seek(0)
        threshold_wig( tss_cvrg, tss_reads_nc_wig_fps, chrm_sizes_fp )
        
        if WRITE_TSS_WIGGLE:
            print >> sys.stderr, 'Writing thresholded tss wiggles...'
            plus_wig = open( 'thresh_tss.plus.wig', 'w' )
            minus_wig = open( 'thresh_tss.minus.wig', 'w' )
            tss_cvrg.write_wiggles( plus_wig, minus_wig )
    
    elif SCALE_WIG and tss_reads_nc_wig_fps != None:
        if VERBOSE:
            print '\tParsing negative control wiggle(s).'
        
        chrm_sizes_fp.seek(0)
        scale_wig( tss_cvrg, tss_reads_nc_wig_fps, chrm_sizes_fp )
        
        if WRITE_TSS_WIGGLES:
            print >> sys.stderr, 'Writing scaled tss wiggles...'
            plus_wig = open( 'scaled_tss.plus.wig', 'w' )
            minus_wig = open( 'scaled_tss.minus.wig', 'w' )
            tss_cvrg.write_wiggles( plus_wig, minus_wig )
    
    return tss_cvrg
    
def get_tss_exons_from_reads( tss_reads_wigs, clustered_exons, nc_wig_fps, chrm_sizes_fp ):
    tss_cvrg = get_tss_wig( chrm_sizes_fp, tss_reads_wigs, nc_wig_fps )
    
    if VERBOSE:
        print '\tFinding TSS exons from reads.'
    all_tss_exons = find_all_distal_exons( clustered_exons, tss_cvrg )
    all_tss_exons = convert_to_genomic_intervals( all_tss_exons )
    
    return all_tss_exons

def get_tss_s( trans_gtfs, tss_gffs, tss_reads_wigs, clustered_exons_gtf, nc_wig_fps, chrm_sizes_fp ):
    tss_exons = set()
    
    if tss_gffs:
        if VERBOSE:
            print '\tParsing TSS exons from TSS file(s).'
        for tss_fp in tss_gffs:
            tss_exons_from_file = parse_exons_file( tss_fp )
            tss_fp.close()
            converted_tss_exons = convert_to_genomic_intervals( tss_exons_from_file )
            tss_exons.update( converted_tss_exons )
    
    if trans_gtfs:
        if VERBOSE:
            print '\tFinding TSS exons from transcripts.'
        # process transcript gtf files
        for trans_fp in trans_gtfs:
            tss_exons_from_trans = get_tss_exons_from_trans( trans_fp )
            trans_fp.close()
            tss_exons.update( tss_exons_from_trans )
    
    if tss_reads_wigs:
        # get exons associated with gene_id
        clustered_exons = parse_clustered_exons_file( clustered_exons_gtf )
        
        tss_exons_from_reads = get_tss_exons_from_reads( \
            tss_reads_wigs, clustered_exons, nc_wig_fps, chrm_sizes_fp )
        
        tss_exons.update( tss_exons_from_reads )
    
    return tss_exons

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description='Get transcript start sites(TSS) from a variety of sources.')
    parser.add_argument( '--transcripts_gtfs', '-t', type=file, nargs='*', \
        help='GTF file which contains exons associated with transcripts. TSS exons will be inferred from exon structure. (eg. flybase, ginormous)')
    parser.add_argument( '--tss_gffs', '-g', type=file, nargs='*', \
        help='GFF file which contains TSS exons. Usually cached from a previous run to this script.')
    parser.add_argument( '--tss_reads_wigs', '-r', type=file, nargs='*', \
        help='WIG file with TSS read coverage. (eg. CAGE, RACE) Note: strand will be infered from filename.')
    parser.add_argument( '--clustered_exons_gtf', '-e', type=file, \
        help="GTF file of exons associated with gene_id's.'")
    parser.add_argument( '--tss_reads_nc_wigs', '-n', type=file, nargs='*', \
        help='Wig files to use as a filtering negative control for CAGE reads.')
    parser.add_argument( '--chrm_sizes_fname', '-c', type=file, \
        help='File with chromosome names and sizes.')
    
    parser.add_argument( '--out_fname', '-o', \
        help='Output file will be written to default. default: stdout')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', \
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    out_fp = open( args.out_fname, "w" ) if args.out_fname else sys.stdout

    global VERBOSE
    VERBOSE = args.verbose

    # ensure that at least one file was provided
    if not any( item != None and len(item) > 0 for item in \
                    [args.transcripts_gtfs, args.tss_gffs, args.tss_reads_wigs] ):
        print 'Must provide one of transcripts_gtfs, tss_gffs or tss_reads_wigs.'
        parser.print_help()
        sys.exit()
    
    if args.tss_reads_wigs != None:
        if args.clustered_exons_gtf == None or args.chrm_sizes_fname == None or \
                ( (THRESHOLD_WIG or SCALE_WIG) and args.tss_reads_nc_wigs == None ):
            print 'WARNING: If tss_reads_gffs file is provided ' + \
                'clustered_exons_gtf and chrm_sizes_fname must also be provided.'
            print '\tNote: if THRESHOLD_WIG or SCALE_WIG flags are on tss_reads_nc_wigs is also required.'
            parser.print_help()
            sys.exit()
    
    return args.transcripts_gtfs, args.tss_gffs, args.tss_reads_wigs, args.clustered_exons_gtf, \
        args.tss_reads_nc_wigs, args.chrm_sizes_fname, out_fp

def main():
    trans_gtfs, tss_gffs, tss_reads_wigs, clustered_exons_gtf, nc_wig_fps, \
        chrm_sizes_fp, out_fp = parse_arguments()
    
    tss_s = get_tss_s( \
        trans_gtfs, tss_gffs, tss_reads_wigs, clustered_exons_gtf, nc_wig_fps, chrm_sizes_fp )
    
    write_TSSs( tss_s, out_fp )
    
    return

if __name__ == "__main__":
    main()

