#!/usr/bin/python

import sys
import os

VERBOSE = False

### Set default parameters ###
CVRG_FRACTION = 0.99
WINDOW_LEN = 50
MIN_CVG_FRAC = 0.50
SMOOTHING_WINDOW_SIZE = 50

from collections import defaultdict

# import slide modules
import distal_exons
from distal_exons import refine_exon_bndrys, convert_to_genomic_intervals, \
    filter_exons_by_max_read_coverage, get_unique_internal_splice_exons

sys.path.append( os.path.join(
        os.path.dirname(__file__), "..", "..", 'file_types' ) )
import wiggle
from clustered_exons_file import parse_clustered_exons_file
from gtf_file import parse_gff_line, iter_gff_lines

def write_TESs( tes_s, out_fp, track_name="tes_exons" ):
    gff_iter = iter_gff_lines( \
        sorted(tes_s), source="disc_tes_from_polya", feature="TES_exon")
    
    out_fp.write( "track name=" + track_name + "\n" )    
    out_fp.write( "\n".join(gff_iter) + "\n" )
    
    return

def build_coverage_wig_from_polya_reads( tes_reads_fps, chrm_sizes_fp ):
    # build lists of the tes locations
    def add_reads_to_locs( tes_reads_fp, locs ):
        for line in tes_reads_fp:
            gff_l = parse_gff_line( line )
            # skip lines that weren't able to be parsed
            if gff_l == None: continue
            
            locs[ (gff_l.region.chr, gff_l.region.strand) ].append(
                gff_l.region.stop )
        
        return
    
    # process tes reads gff files (i.e. poly-A)
    locs = defaultdict( list )
    for tes_reads_fp in tes_reads_fps:
        add_reads_to_locs( tes_reads_fp, locs )
    
    # load locs data into a wiggle object
    tes_cvrg = wiggle.Wiggle( chrm_sizes_fp )
    tes_cvrg.load_data_from_positions( locs )
        
    return tes_cvrg

def parse_arguments():
    # global variables that can be set by arguments
    global CVRG_FRACTION
    
    import argparse
    parser = argparse.ArgumentParser(
        description='Get transcript end sites (TESs) from polya reads.')
    parser.add_argument(
        '--chrm-sizes', '-c', required=True, type=file,
        help='File with chromosome names and sizes.')
    parser.add_argument(
        '--clustered-exons-gtf', '-e', required=True, type=file,
        help="GTF file of exons associated with gene_id's.'")
    parser.add_argument(
        '--polya-read-gffs', '-r', required=True, type=file, nargs='+',
        help='GFF file which contains reads ending at a TES.')
    parser.add_argument(
        '--coverage-fraction', type=float, default=CVRG_FRACTION,
        help='Fraction of TES coverage to include in the TES refined exon '
        + 'boundaries. Default: %(default)f' )    
    parser.add_argument(
        '--out-fname', '-o',
        help='Output file will be written to default. default: stdout')
    parser.add_argument(
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    out_fp = open( args.out_fname, "w" ) if args.out_fname else sys.stdout
    
    global VERBOSE
    VERBOSE = args.verbose
    distal_exons.VERBOSE = VERBOSE
    wiggle.VERBOSE = VERBOSE
    
    # Set parameter arguments
    CVRG_FRACTION = args.coverage_fraction
    
    return args.clustered_exons_gtf, args.chrm_sizes,\
        args.polya_read_gffs, out_fp

def main():
    # parse arguments
    clustered_exons_fp, chrm_sizes_fp, polya_read_gff_fps, out_fp, \
        = parse_arguments()
    
    if VERBOSE: print >> sys.stderr, "Loading clustered exons."
    clustered_exons = parse_clustered_exons_file( clustered_exons_fp )
    
    if VERBOSE: print >> sys.stderr, "Loading TES coverage wiggle."
    tes_cvg = build_coverage_wig_from_polya_reads(
        polya_read_gff_fps, chrm_sizes_fp )
            
    smooth_tes_cvg = tes_cvg.get_smoothed_copy( SMOOTHING_WINDOW_SIZE )

    tes_exons = {}
    for (chrm, strand) in clustered_exons:
        if (chrm, strand) not in tes_cvg:
            print >> sys.stderr, \
                "WARNING: '%s:%s' does not exist in " % ( chrm, strand ) \
                + "the polya signal wiggle. Skipping that chr strand pair." 
            continue
        tes_exons[ (chrm, strand) ] = []
        
        exons_clusters = clustered_exons[ (chrm, strand) ]

        for exons in exons_clusters:
            unique_exons = get_unique_internal_splice_exons( \
                exons, strand, is_tss=False )
            
            try:
                filtered_exons = filter_exons_by_max_read_coverage( \
                    exons, tes_cvg[(chrm, strand)], \
                        WINDOW_LEN, MIN_CVG_FRAC )
            except ValueError, inst:
                if str( inst ) != "No signal.":
                    raise
                err_str = "WARNING: Region %s %s %s %s had no POLYA signal." \
                    % (chrm, strand, exons[0][0], exons[-1][1])
                print >> sys.stderr, err_str
                filtered_exons = []
            
            refined_exons = refine_exon_bndrys( \
                filtered_exons, tes_cvg, (chrm, strand), False, CVRG_FRACTION )
            
            tes_exons[ (chrm, strand) ].extend( refined_exons )
    
    # ocver the exons to genomic intervals
    tes_exons = convert_to_genomic_intervals( tes_exons )
    
    if VERBOSE: print >> sys.stderr, 'Writing TES exons to a gff file.'    
    write_TESs( tes_exons, out_fp )
    out_fp.close()
    
    return

if __name__ == "__main__":
    main()
