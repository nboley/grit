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
THRESHOLD_FRAC = 0.90
SMOOTHING_WINDOW_SIZE = 50

# this means filter exons where the ratio of the highest 50 bps of 
# cage signal in exon to the MAX exon is < 0.20. 
REGION_FILTER_LEN = 50
EXON_FILTER_TYPE = 'region'
MIN_SCORE_RATIO = 0.20

# set parameters for distal exon finding algorithm
import distal_exons
distal_exons.MIN_SCORE_RATIO = 0.01
distal_exons.CVRG_FRACTION = 0.99
distal_exons.SMOOTHING_WINDOW_SIZE = None
distal_exons.RATIO_THRESH = None


#############################################################

# import slide modules
from distal_exons import get_unique_internal_splice_exons, \
    refine_exon_bndrys, convert_to_genomic_intervals

sys.path.append(os.path.join(os.path.dirname(__file__), "../../", 'file_types'))
import wiggle
from clustered_exons_file import parse_clustered_exons_file
from exons_file import parse_exons_file

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
    parser.add_argument( '--out-fname', '-o', \
        help='Output file will be written to default. default: stdout')
    parser.add_argument( '--write-intermediate-wiggles', default=False, \
                             action="store_true",\
        help='Whether or not to write intermediate wiggles ( mostly for debugging ).')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    out_fp = open( args.out_fname, "w" ) if args.out_fname else sys.stdout
    
    global VERBOSE
    VERBOSE = args.verbose
    distal_exons.VERBOSE = VERBOSE
    wiggle.VERBOSE = VERBOSE
    
    return args.clustered_exons_gtf, args.chrm_sizes, \
        args.cage_read_wigs, out_fp, \
        args.write_intermediate_wiggles

def filter_exons_by_quantiles( exons, cvg, smooth_cage_read_cvg, chrm, strand ):
    """The approach is as follows:

    1) Find promoter exons:
    - smooth cage signal with SMOOTH_SIZE bp smoothing window
    - order smoothed arrays by size
    - threshold every position in the lower RATIO_THRESH% ( by *signal* )
    - keep exons that still have cage signal
    
    """
    # concatinate the smmothed cage signal in all of the exons
    smooth_cvg_arrays = []
    nonoverlapping_signal = numpy.zeros(0)
    # we keep track of the previous stop position to 
    # make sure that we don't double count overlapping exons.
    # Since the starts are sorted, the only way the next exon
    # can overlap is if the stop of the previous exon is greater 
    # than the start of the next exon. So, we simply take the max of
    # the two, with a special case in case the exons are sub-sumed. 
    prev_stop = -1
    for start, stop in sorted(exons):
        smooth_cvg_arrays.append( \
            smooth_cage_read_cvg[(chrm, strand)][ start:stop+1 ])
        # if there is actually some new sequence, add it
        if stop > prev_stop:
            nonoverlapping_signal = numpy.append(
                nonoverlapping_signal,
                smooth_cage_read_cvg[(chrm, strand)][ \
                    max( start, prev_stop ):stop+1 ]
            )
        
        prev_stop = max( stop, prev_stop )
    
    nonoverlapping_signal.sort()
    thresh = THRESHOLD_FRAC*nonoverlapping_signal.sum()

    # if there are no cage reads, then skip this cluster
    if nonoverlapping_signal[-1] == 0: 
        assert thresh == 0
        return []

    # find how many bases we need to observe to account for at least
    # THRESHOLD_FRAC of the total CAGE coverage
    tot = 0
    for thresh_val in reversed(nonoverlapping_signal):
        tot += thresh_val
        if tot > thresh: 
            break

    # threshold bases that are below the threshold in each exon
    new_exons = []
    for i, array in enumerate( smooth_cvg_arrays ):
        # get the indices to threshold
        bad_indices = ( array < thresh_val ).nonzero()[0]
        shifted_bad_indices = start + bad_indices
        if len( shifted_bad_indices ) < len( array ):
            new_exons.append( exons[i] )
    
    return new_exons

def filter_exons_by_max_read_coverage( exons, cvg_array ):
    """Remove exons which contain less than MIN_SCORE_RATIO distal read coverage
    """
    def get_exon_score( exon ):
        start, stop = exon
        if stop - start + 1 <= REGION_FILTER_LEN or \
                EXON_FILTER_TYPE in [ 'total', 'avg' ]:
            # get total coverage from TSS_coverage
            score = cvg_array[start:stop+1].sum()
            
            if EXON_FILTER_TYPE == 'avg':
                # get average distal read coverage over exon
                score = score / (stop - start + 1)
        else:
            assert EXON_FILTER_TYPE == 'region'
            # get maximal tss coverage region accoss exon
            cumsum_cvg_array = \
                numpy.append(0, numpy.cumsum( cvg_array[start : stop + 1]))
            # get every coverage total over REGION_FILTER_LEN regions
            score = ( cumsum_cvg_array[REGION_FILTER_LEN:] \
                         - cumsum_cvg_array[:-REGION_FILTER_LEN] ).max()
        return score
    
    def get_exon_scores( exons ):
        max_score = 0
        exons_w_scores = []
        for exon in exons:
            score = get_exon_score( exon )
            max_score = max( max_score, score )
            exons_w_scores.append( (exon, score) )
        
        return exons_w_scores, max_score
    
    def get_candidate_exons( exons_w_scores, max_score ):
        """remove exons with sufficiently low distal read coverage
        """
        # the non_inclusive minimum score is scaled to the 
        # ceiling of the max_score * MIN_SCORE_RATIO
        # i.e. if the ratio is 0.01 and the max score is 30
        # a score greater than one will suffice for *high score*
        min_score = int( MIN_SCORE_RATIO * max_score ) + 1
        
        high_score_exons = []
        for exon, score in exons_w_scores:
            if score > min_score:
                high_score_exons.append( exon )
        
        return high_score_exons
    
    
    exons_w_scores, max_score = get_exon_scores( exons )
    candidate_exons = get_candidate_exons( exons_w_scores, max_score )
    
    return candidate_exons

def main():
    # parse arguments
    clustered_exons_fp, chrm_sizes_fp, cage_wig_fps, out_fp, \
        write_intermediate_wiggles = parse_arguments()
    
    if not write_intermediate_wiggles:
        intermediate_wiggles_prefix = None
    else:
        intermediate_wiggles_prefix = out_fname
    
    if VERBOSE: print >> sys.stderr, "Loading clustered exons."
    clustered_exons = parse_clustered_exons_file( clustered_exons_fp )

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
                exons, strand, is_tss=True )

            
            filtered_exons_1 = filter_exons_by_max_read_coverage( \
                unique_exons, cage_read_cvg[(chrm, strand)] )
            
            filtered_exons_2 = filter_exons_by_quantiles( \
                unique_exons, cage_read_cvg, smooth_cage_read_cvg, chrm, strand)
            
            #filtered_exons = sorted( set(filtered_exons_1).intersection(\
            #        filtered_exons_2) )
            filtered_exons = sorted( set(filtered_exons_1 + filtered_exons_2) )

            refined_exons = refine_exon_bndrys( \
                filtered_exons, cage_read_cvg, (chrm, strand), is_tss=True )

            tss_exons[ (chrm, strand) ].extend( refined_exons )
            
    tss_exons = convert_to_genomic_intervals( tss_exons )
    
    if VERBOSE: print >> sys.stderr, 'Writing TSS exons to a gff file.'
    write_TSSs( tss_exons, out_fp )
    out_fp.close()
    
    return

if __name__ == "__main__":
    main()
