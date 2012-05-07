#!/usr/bin/python

import sys
import os
import numpy

VERBOSE = True
MIN_VERBOSE = True

# parameters for negative control thresholding
# HOW TO ESTIMATE THIS?
# this means threshold all bases with a score below
# X where X is the score at which THRESHOLD_FRAC of the cage 
# signal falls below
THRESHOLD_FRAC = 0.70
SMOOTHING_WINDOW_SIZE = 120

# this means filter exons where the ratio of the highest 50 bps of 
# cage signal in exon to the MAX exon is < MIN_CVG_FRAC. 
REGION_FILTER_LEN = 120
EXON_FILTER_TYPE = 'region'
MIN_CVG_FRAC = 0.10

# the fraction of CAGE reads to include in the boundary refined TSS exon
CVRG_FRACTION = 0.99

MIN_WIN_CVG = 1
#############################################################

# import slide modules
import distal_exons
from distal_exons import get_unique_internal_splice_exons, \
    refine_exon_bndrys, convert_to_genomic_intervals, \
    filter_exons_by_quantiles, filter_exons_by_max_read_coverage, \
    filter_clustered_exons

sys.path.append(os.path.join(os.path.dirname(__file__), "../../", 'file_types'))
import wiggle
from clustered_exons_file import parse_clustered_exons_file
from exons_file import parse_exons_file
from junctions_file import get_intron_bndry_sets

from gtf_file import parse_gff_line, iter_gff_lines

def write_TSSs( tss_s, out_fp, track_name="tss_exons" ):
    gff_iter = iter_gff_lines( \
        sorted(tss_s), source="disc_tss_from_CAGE", feature="TSS_exon")
    
    out_fp.write( "track name=" + track_name + "\n" )    
    out_fp.write( "\n".join(gff_iter) + "\n" )

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(
        description='Get transcript end sites (TESs) from polya reads.')
    parser.add_argument( '--chrm-sizes', '-c', required=True, type=file, \
        help='File with chromosome names and sizes.')
    parser.add_argument( '--clustered-exons-gtf', '-e', required=True, \
        type=file, help="GTF file of exons associated with gene_id's.'")
    parser.add_argument( '--cage-read-wigs', '-r', required=True, type=file, \
        nargs='+', help='Wiggle files containing CAGE read coverage.')
    parser.add_argument( '--junctions-gff', '-j', required=False, type=file,
        help='GFF file which contains junctions to discover TSS exons ' \
            + 'in the absense of CAGE data.')
    parser.add_argument( '--out-fname', '-o', \
        help='Output file will be written to default. default: stdout')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    out_fp = open( args.out_fname, "w" ) if args.out_fname else sys.stdout
    
    global VERBOSE
    VERBOSE = args.verbose
    distal_exons.VERBOSE = VERBOSE
    wiggle.VERBOSE = VERBOSE
    
    return args.clustered_exons_gtf, args.chrm_sizes, \
        args.cage_read_wigs, args.junctions_gff, out_fp

def main():
    # parse arguments
    clustered_exons_fp, chrm_sizes_fp, cage_wig_fps, jns_fp, out_fp, \
        = parse_arguments()
    
    if VERBOSE: print >> sys.stderr, "Loading clustered exons."
    clustered_exons = parse_clustered_exons_file( clustered_exons_fp )

    if jns_fp != None:
        upstream_bnds, downstream_bnds = get_intron_bndry_sets( jns_fp )
        jns = upstream_bnds
    else:
        jns = None
    
    if VERBOSE: print >> sys.stderr, "Loading CAGE coverage wiggle."
    cage_read_cvg = wiggle.Wiggle( chrm_sizes_fp, cage_wig_fps )
    smooth_cage_read_cvg = \
        cage_read_cvg.get_smoothed_copy( SMOOTHING_WINDOW_SIZE )

    tss_exons = {}
    for (chrm, strand) in clustered_exons:
        if (chrm, strand) not in cage_read_cvg:
            print >> sys.stderr, \
                "WARNING: '%s:%s' does not exist in " % ( chrm, strand ) \
                + "the cage signal wiggle. Skipping that chr strand pair." 
            continue
        
        tss_exons[ (chrm, strand) ] = []
        
        exons_clusters = clustered_exons[ (chrm, strand) ]

        for exons in exons_clusters:
            unique_exons = get_unique_internal_splice_exons( \
                exons, strand, is_tss=True, max_bndry_exp=200 )
            
            filtered_exons = filter_clustered_exons( \
                unique_exons, cage_read_cvg, smooth_cage_read_cvg, jns, \
                REGION_FILTER_LEN, MIN_CVG_FRAC, THRESHOLD_FRAC, \
                MIN_WIN_CVG, chrm, strand )
            
            refined_exons = refine_exon_bndrys( \
                filtered_exons, cage_read_cvg,  \
                    (chrm, strand), True, CVRG_FRACTION )
            tss_exons[ (chrm, strand) ].extend( refined_exons )
            
    tss_exons = convert_to_genomic_intervals( tss_exons )
    
    if VERBOSE: print >> sys.stderr, 'Writing TSS exons to a gff file.'
    write_TSSs( tss_exons, out_fp )
    out_fp.close()
    
    return

if __name__ == "__main__":
    main()
