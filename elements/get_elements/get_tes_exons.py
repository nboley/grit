#!/usr/bin/python

import sys
import os

VERBOSE = False
WRITE_TES_WIGGLES = False

# parameters for negative control thresholding or scaling
THRESHOLD_WIG = False
SCALE_WIG = False

# set parameters for distal exon finding algorithm
import distal_exons
distal_exons.EXON_FILTER_TYPE = 'region'
distal_exons.REGION_FILTER_LEN = 200
distal_exons.MIN_SCORE_RATIO = 0.01
distal_exons.CVRG_FRACTION = 0.95
distal_exons.SMOOTHING_WINDOW_SIZE = 70
distal_exons.RATIO_THRESH = 50.0
distal_exons.VERBOSE = VERBOSE

#############################################################

from collections import defaultdict

# import slide modules
from distal_exons import \
    find_all_distal_exons, convert_to_genomic_intervals, threshold_wig

sys.path.append( os.path.join( os.path.dirname(__file__), ".." ) )
import slide
slide.VERBOSE = VERBOSE
from slide import get_raw_transcripts
from gene_models import parse_gff_line, GenomicInterval

sys.path.append( os.path.join( os.path.dirname(__file__), "..", 'file_types' ) )
from wiggle import wiggle
from clustered_exons_file import parse_clustered_exons_file
from exons_file import parse_exons_file


def create_tes_gff_line( tes ):
    gff_line = [ 'chr' + tes.chr, 'elements', 'TES', str(tes.start), \
                     str(tes.stop), '.', tes.strand, '.' ]
    return '\t'.join( gff_line )

def write_TESs( tes_s, out_fp ):
    for tes in sorted( tes_s ):
        out_fp.write( create_tes_gff_line( tes ) + '\n' )
    out_fp.close()

def get_tes_exons_from_trans( trans_fp ):
    tes_s = set()
    
    # get the gene -> trans ->  exons structure from the gtf file
    # verification of transcripts is done from within get_raw_transcripts
    raw_gene_transcripts = get_raw_transcripts( trans_fp )
    for gene_name, transcripts in raw_gene_transcripts.iteritems():
        for trans_name, exons in transcripts.iteritems():
            if exons[0].strand == '.': continue
            if exons[0].strand == '+':
                tes = GenomicInterval( \
                    exons[0].chr, '+', exons[-1].start, exons[-1].stop )
            else:
                assert exons[0].strand == '-'
                tes = GenomicInterval( \
                    exons[0].chr, '-', exons[0].start, exons[0].stop )
            tes_s.add( tes )
    
    return tes_s

def build_coverage_wig_from_polya_starts( tes_reads_fps, chrm_sizes_fp ):
    # build lists of the tes locations
    locs = defaultdict( list )
    def add_read_locs( tes_reads_fp ):
        for line in tes_reads_fp:
            gene_name, trans_name, feature_type, region = parse_gff_line( line )
            if region != None:
                if region.strand == '+':
                    locs[ (region.chr, '+') ].append( region.stop )
                else:
                    if region.strand != '-':
                        print 'GTF line does not have valid strand: ' + line
                        continue
                    locs[ (region.chr, '-') ].append( region.start )
        tes_reads_fp.close()
        
        return
    
    
    # process tes reads gff files (i.e. poly-A)
    for tes_reads_fp in tes_reads_fps:
        add_read_locs( tes_reads_fp )
    
    # load locs data into a wiggle object
    tes_cvrg = wiggle( chrm_sizes_fp )
    tes_cvrg.load_data_from_positions( locs )
    
    if WRITE_TES_WIGGLES:
        print >> sys.stderr, 'Writing unthresholded tes wiggles out...'
        plus_wig = open( 'unthresh_tes.plus.wig', 'w' )
        minus_wig = open( 'unthresh_tes.minus.wig', 'w' )
        tes_cvrg.write_wiggles( plus_wig, minus_wig )
    
    return tes_cvrg

def get_tes_wig( chrm_sizes_fp, tes_reads_fps, tes_reads_nc_wig_fps ):
    if VERBOSE:
        print '\tParsing TES coverage wiggle(s).'
    tes_cvrg = build_coverage_wig_from_polya_starts( \
        tes_reads_fps, chrm_sizes_fp )
    
    if THRESHOLD_WIG and tes_reads_nc_wig_fps != None:
        if VERBOSE:
            print '\tParsing negative control wiggle(s).'
        
        chrm_sizes_fp.seek(0)
        threshold_wig( tes_cvrg, tes_reads_nc_wig_fps, chrm_sizes_fp )
        
        if WRITE_TES_WIGGLES:
            print >> sys.stderr, 'Writing thresholded tes wiggles...'
            plus_wig = open( 'thresh_tes.plus.wig', 'w' )
            minus_wig = open( 'thresh_tes.minus.wig', 'w' )
            tes_cvrg.write_wiggles( plus_wig, minus_wig )
    
    elif SCALE_WIG and tes_reads_nc_wig_fps != None:
        if VERBOSE:
            print '\tParsing negative control wiggle(s).'
        
        chrm_sizes_fp.seek(0)
        scale_wig( tes_cvrg, tes_reads_nc_wig_fps, chrm_sizes_fp )
        
        if WRITE_TES_WIGGLES:
            print >> sys.stderr, 'Writing scaled tes wiggles...'
            plus_wig = open( 'scaled_tes.plus.wig', 'w' )
            minus_wig = open( 'scaled_tes.minus.wig', 'w' )
            tes_cvrg.write_wiggles( plus_wig, minus_wig )
    
    return tes_cvrg

def get_tes_exons_from_reads( \
        tes_reads_fps, clustered_exons, nc_wig_fps, chrm_sizes_fp ):
    tes_cvrg = get_tes_wig( chrm_sizes_fp, tes_reads_fps, nc_wig_fps )
    
    if VERBOSE:
        print '\tFinding TES exons from reads.'
    tes_exons = find_all_distal_exons( clustered_exons, tes_cvrg, False )
    if VERBOSE:
        print '\tFound distal exons.'    
    tes_exons = convert_to_genomic_intervals( tes_exons )
    if VERBOSE:
        print '\Converted distal exons to genomic intervals.'
    return tes_exons

def get_tes_s( trans_gtfs, tes_gffs, tes_reads_gffs, \
                   clustered_exons_fp, nc_wig_fps, chrm_sizes_fp ):
    tes_exons = set()
    
    if tes_gffs:
        if VERBOSE:
            print '\tParsing TES exons from TES file(s).'
        for tes_fp in tes_gffs:
            tes_exons_from_file = parse_exons_file( tes_fp )
            tes_fp.close()
            converted_tes_exons = convert_to_genomic_intervals( tes_exons_from_file )
            tes_exons.update( converted_tes_exons )
    
    if trans_gtfs:
        if VERBOSE:
            print '\tFinding TES exons from transcripts.'
        # process transcript gtf files
        for trans_fp in trans_gtfs:
            tes_exons_from_trans = get_tes_exons_from_trans( trans_fp )
            trans_fp.close()
            tes_exons.update( tes_exons_from_trans )
    
    if tes_reads_gffs:
        # get exons associated with gene_id
        clustered_exons = parse_clustered_exons_file( clustered_exons_fp )
        
        tes_exons_from_reads = get_tes_exons_from_reads( \
            tes_reads_gffs, clustered_exons, nc_wig_fps, chrm_sizes_fp )
        
        tes_exons.update( tes_exons_from_reads )
    
    return tes_exons

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description='Get transcript end sites(TES) from a variety of sources.')
    parser.add_argument( '--transcripts_gtfs', '-t', type=file, nargs='*', \
        help='GTF file which contains exons associated with transcripts. TESs will be inferred from exon structure. (eg. flybase, ginormous)')
    parser.add_argument( '--tes_gffs', '-g', type=file, nargs='*', \
        help='GFF file which contains TES exons. Usually cached from a previous run to this script.')
    parser.add_argument( '--tes_reads_gffs', '-r', type=file, nargs='*', \
        help='GFF file which contains reads ending at a TES. (eg. poly-A reads)')
    parser.add_argument( '--clustered_exons_gtf', '-e', type=file, \
        help="GTF file of exons associated with gene_id's.'")
    parser.add_argument( '--tes_reads_nc_wigs', '-n', type=file, nargs='*', \
        help='Wig files to use as a filtering negative control for poly-A reads.')
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
                    [args.transcripts_gtfs, args.tes_gffs, args.tes_reads_gffs] ):
        print 'Must provide one of transcripts_gtfs, tes_gffs or tes_reads_wigs.'
        parser.print_help()
        sys.exit()
    
    if args.tes_reads_gffs != None:
        if args.clustered_exons_gtf == None or args.chrm_sizes_fname == None or \
                ( ( THRESHOLD_WIG or SCALE_WIG ) and args.tes_reads_nc_wigs == None):
            print 'WARNING: If tes_reads_gffs file is provided ' + \
                'clustered_exons_gtf and chrm_sizes_fname must also be provided.'
            print '\tNote: if THRESHOLD_WIG or SCALE_WIG flags are on tes_reads_nc_wigs is also required.'
            parser.print_help()
            sys.exit()
    
    return args.transcripts_gtfs, args.tes_gffs, args.tes_reads_gffs, args.clustered_exons_gtf, \
        args.tes_reads_nc_wigs, args.chrm_sizes_fname,  out_fp

def main():
    trans_gtfs, tes_gffs, tes_reads_gffs, exons_fp, tes_reads_nc_wigs, \
        chrm_sizes_fp, out_fp = parse_arguments()
    
    if VERBOSE:
        print >> sys.stderr, 'Getting TESs...'
    tes_s = get_tes_s( \
        trans_gtfs, tes_gffs, tes_reads_gffs, \
            exons_fp, tes_reads_nc_wigs, chrm_sizes_fp )
    
    if VERBOSE:
        print >> sys.stderr, 'Writing TESs...'
    write_TESs( tes_s, out_fp )
    
    return

if __name__ == "__main__":
    main()

