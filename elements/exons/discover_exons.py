import sys, os
import copy
import numpy
from collections import defaultdict, namedtuple
from itertools import izip, chain, takewhile, count
from scipy import median
from build_genelets import cluster_exons
from math import log

# since some mismapped signal spills into intron
# this many bases are ignored when finding 
# ONE_SIDE_EXTENTION_FRACTION of cvrg
BASES_TO_IGNORE_IN_RETAIN = 15
# region adjacent to boundaries defined by junctions and empty space over 
# which to find collect cvrg stats for retained intron signal
BIN_BOUNDRY_SIZE = 20
# maximum ratio for retained intron bndry cvrg stats
MAX_RETAIN_RATIO = 10.0
# fraction of coverage within intron to include 
# in exon extention in a left/right retained intron
ONE_SIDE_EXTENTION_FRACTION = 0.95
# stop the intron when the average ( rather than the 
# total ) falls below some threshold
USE_AVERAGE_FOR_ONE_SIDE_RETAIN = False

SPLIT_EXON_BNDRY_COV_FRAC = 0.95
EXON_SPLIT_SIZE = 40
EXON_SPLIT_RATIO = 5

# length of contiguous space to define as an empty region
EMPTY_REGION_SPLIT_SIZE = 80
# how many bases to move the boundaries in,
# after we have dtermined that it is empty
# set to 0 to not do this
EMPTY_REGION_BNDRY_SHRINK = 0

MAX_EXON_SIZE = 50000
MIN_EXON_LEN = 15

MIN_SE_GENE_LEN = 500
MIN_SE_GENE_AVG_READCOV = 10

# LOC_THRESH_FRAC = 0.20
LOC_THRESH_REG_SZ = 50000

### CAGE TUUNING PARAMS
CAGE_WINDOW_LEN = 120
CAGE_MAX_SCORE_FRAC = 0.10
CAGE_MIN_SCORE = 5
CAGE_SUFFICIENT_SCORE = 200
# this is set in main
CAGE_TOT_FRAC = None


### PolyA TUNING PARAMS
POLYA_WINDOW_LEN = 10
POLYA_MAX_SCORE_FRAC = 0.001
POLYA_MIN_SCORE = 2
# this is set in main
POLYA_TOT_FRAC = None

NORMALIZE_BY_RNASEQ_COV = True

sys.path.append(os.path.join(os.path.dirname(__file__), "../../", 'file_types'))
from wiggle import Wiggle
from junctions_file import parse_junctions_file_dont_freeze
from gtf_file import parse_gff_line, iter_gff_lines, GenomicInterval

BinStats = namedtuple('BinStats', ['is_small', 'mean', 'lmn', 'rmn'] )

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
    tes_cvrg = Wiggle( chrm_sizes_fp )
    tes_cvrg.load_data_from_positions( locs )
        
    return tes_cvrg


class ReadCoverageData( object ):
    def __init__( self, zero_intervals, read_coverage_array ):
        self.zero_intervals = zero_intervals
        self.chrm_stop = len( read_coverage_array )
        
        self.read_coverage_array = read_coverage_array
        self.rca = self.read_coverage_array
        
        return
    
    def is_in_empty_region( self, pos ):
        # if the actual position is not zero, then this can't be
        # an empty region. ( necessary to be empty but not sufficient  )
        if self.read_coverage_array[ pos ] > 0: 
            return False
        
        # otherwise, we do a searchsorted into the zero_intervals array
        index = self.zero_intervals[:,0].searchsorted( pos, side='right' )
        
        # if the pos is before the first zero_interval or if the 
        # position does not fall within a zero region
        if index == 0 or pos > self.zero_intervals[index-1][1]:
            return False
        else:
            return True
    
    def iter_empty_regions( self ):
        # Note that zero_intervals is sorted
        # If the implimentation changes ensure that it stays sorted
        for start, stop in self.zero_intervals:
            yield start, stop
        
        return

def get_possible_exon_bndrys( labels, boundaries, chrm_stop ):
    """Get all of the possible exon boundaries corresponding to contiguous 'Ee's
    Invalid exons (starting with intron_start or ending with intron_stop) are 
    removed later.
    """
    # get all labels that are either 'Ee'
    exon_indices = []
    for index, label in enumerate( labels ):
        if label in ('TSS', 'TES', 'T_S', 'Exon', 'exon_extension'):
            exon_indices.append( index )
    
    if len( exon_indices ) == 0:
        return []
    
    # find and group contiguous 'Ee' regions
    exon_grps = [ [exon_indices[0],], ]
    for index in exon_indices[1:]:
        if index - exon_grps[-1][-1] == 1:
            exon_grps[-1].append( index )
        else:
            exon_grps.append( [ index,] )
    
    bndrys = []
    for exon_grp in exon_grps:
        # look at all combinations of E or e's in a row
        for start_i, start_exon in enumerate(exon_grp):
            for stop_i, stop_exon in enumerate(exon_grp[start_i:]):
                stop_i += start_i
                # skip grps that contain zero 'E''s
                if not any( labels[index] in ('Exon', 'TSS', 'T_S', 'TES') \
                                for index in exon_grp[start_i:stop_i+1] ):
                    continue
                # if this is the last region goes to the end of the chrm then
                # the last exon goes up to the end of the chromosome
                if stop_exon+1 == len(boundaries):
                    bndrys.append( ( boundaries[start_exon], chrm_stop ) )
                else:
                    bndrys.append( ( boundaries[start_exon], \
                                         boundaries[stop_exon+1]-1 ) )
    
    return bndrys

def get_possible_single_exon_genes( labels, boundaries, chrm_stop ):
    se_bndrys = []
    # add the single exon genes
    for index, label in enumerate( labels[:-1] ):
        if label != 'S':
            continue        
        se_bndrys.append( ( boundaries[index], \
                            boundaries[index+1]-1 ) )
            
    return se_bndrys

def build_bins_stats( boundaries, labels, read_coverages ):
    def build_bin_stats( region_cvrg, start, stop):
        size = len( region_cvrg )
        assert size > 0
        is_small = size < BIN_BOUNDRY_SIZE
        
        mean = region_cvrg.mean()
        
        # if the region is too small, try to calculate the 
        # lmn and rmn by extending past the boundaries if 
        # this does not extend past the end of the chromosome
        if is_small:
            lmn = read_coverages.rca[start:start+BIN_BOUNDRY_SIZE].mean() if \
                len( read_coverages.rca) > start+BIN_BOUNDRY_SIZE+1 else \
                mean
            rmn = read_coverages.rca[stop-BIN_BOUNDRY_SIZE:stop+1].mean() if \
                stop - BIN_BOUNDRY_SIZE >= 0 else \
                mean
        else:
            lmn = region_cvrg[:BIN_BOUNDRY_SIZE].mean()
            rmn = region_cvrg[-BIN_BOUNDRY_SIZE:].mean()
        
        return BinStats( is_small, mean, lmn, rmn )
    
    bins_stats = []
    for start, stop in izip( boundaries[:-1], boundaries[1:] ):
        # do not add one to stop b/c stop is the start of the next bndry
        # and since slicing is closed open start:stop gives the coverage
        # values for this region correctly
        region_cvrg = read_coverages.rca[start:stop].copy()
        bins_stats.append( build_bin_stats( region_cvrg, start, stop ) )
    
    if boundaries[-1] > len( read_coverages.rca ):
        # calculate stats for the final region
        last_region = read_coverages.rca[boundaries[-1]:].copy()
        bins_stats.append( build_bin_stats( last_region ) )
    else:
        # the last boundary is the end of the chrom
        bins_stats.append( BinStats(True,0,0,0) )
    
    return bins_stats

def is_right_retained( intron_stats, right_exon_stats ):
    if right_exon_stats.lmn/( intron_stats.rmn + 1e-8 ) \
            < MAX_RETAIN_RATIO:
        return True
    
    return False

def is_left_retained( left_exon_stats, intron_stats ):
    if left_exon_stats.rmn/( intron_stats.lmn + 1e-8 ) \
            < MAX_RETAIN_RATIO:
        return True
    
    return False

def calc_partially_retained_intron_regions(
    start, stop, read_coverages, is_left_retained, \
        one_side_retained_intron_coverage_frac, \
        edge_bases_to_ignore ):
    """If we have a partially retained intron ( ie code l or r )
    then find how many bases belong int he intron region, and 
    how many belong int he exon extension region.
    
    Args:
    start: the start of the region to split
    stop: the stop ( inclusive ) of the region to split
    
    Returns:
    a tuple of the form 
        ( start_1, stop_1, new_label_1 ), ( start_2, stop_2, new_label_2 ) )

    where stop_1 + 1 == start_2
    and new_label's are the new labels to assign to these subregions
    
    """
    # if this region is too small make it a potentially retained intron
    if stop - start + 1 < edge_bases_to_ignore:
        return (start,), 'exon_extension'
    
    # get an iter over the bases that are in this region. 
    # we make the intron portion the portion that contains
    # one_side_retained_intron_coverage_frac one the total
    # reads in the region. If the intron is left retained, 
    # that means that we want to start counting reads from 
    # the left: if it's right retained we want to count from
    # the right side. This means that ( for efficiency )
    # we can start counting from the opposite side,
    # and stop when we've counted 1 - FRAC of the reads. 

    # in addition, we ignore the first edge_bases_to_ignore
    # to account for read coverage spillover from mapping 
    # ( ie, the first 5 bases of a read may map in the intron
    # instead of the other side of the splice )
    region_cvg = read_coverages.rca[start:stop+1]

    # make the correction for right vs left retained
    if is_left_retained: region_cvg = region_cvg[::-1]
    # make the correction for mapping runover
    region_cvg = region_cvg[edge_bases_to_ignore:]
    
    if not USE_AVERAGE_FOR_ONE_SIDE_RETAIN:
        read_cutoff = region_cvg.sum() * \
            ( 1 - one_side_retained_intron_coverage_frac )
        cumsum = 0
        for num_intron_bases, cvg in enumerate(region_cvg):
            cumsum += cvg
            if cumsum > read_cutoff:
                break
    else:
        # we add 1 to the number of intron bases to accoutn for the 
        # fact that enumerate begins from 0
        region_len = len( region_cvg )
        region_cvg_cumsum = region_cvg.cumsum()
        total_cvg = region_cvg_cumsum[-1]
        for num_intron_bases, cvg in enumerate(region_cvg_cumsum):
            if cvg < 1e-8: continue
            left_avg_cvg = cvg/(num_intron_bases+1)
            right_avg_cvg = (total_cvg - cvg)/( \
                region_len - num_intron_bases + 1 )
            if left_avg_cvg/(right_avg_cvg+1e-6) \
                    < one_side_retained_intron_coverage_frac:
                break
        num_intron_bases += 1
    
    if is_left_retained:
        return ( start, stop - num_intron_bases + 1), ("exon_extension","I")
    else:
        return ( start, start + num_intron_bases + 1), ("I","exon_extension")
    
    return

def check_for_retained_intron( left, mid, right):
    # if we have and intron surrounded by low coverage regions,
    # then keep it an intron
    if left[0] == 'L' and right[0] == 'L':
        return 'I'
    
    # if this isn't next to an exon, it can't be retained
    #if left[0] in ('L', ) \
    #        and right[0] in ('L', ):
    #    return 'I'
    
    # if only the region to the right is an exon
    if left[0] == 'L':
        # hck if the signal ratios are high enough for this to be retained
        if is_right_retained( mid[1], right[1] ):
            return 'exon_extension'
        
        return 'I'
    
    # if only the region to the left is an exon
    if right[0] in 'L':
        # check if signal ratios are high enough for this to be retained
        if is_left_retained( left[1], mid[1] ):
            return 'exon_extension'
        
        return 'I'
    
    # if this is an intron surrounded by exons.
    if left[0] in ('exon_extension', 'Exon') \
            and right[0] in ('exon_extension', 'Exon'):
        # determine if this region is possibly a retained intron
        # with respect ot both left and right boundaries
        left_retained = is_left_retained( left[1], mid[1] )
        right_retained = is_right_retained( mid[1], right[1] )

        # if the intron is short, and it's right_retained *or*
        # left retained, then it is retained.
        # So change to exon extension
        if mid[1].is_small and ( right_retained or left_retained ):
            return 'exon_extension'
        # if it's right retained *and* left retained, then 
        # mark this as an exon extension
        elif right_retained and left_retained:
            return 'exon_extension'
        # otherwise, we try and split the exon if it's either 
        # right retained or left retained but not both
        elif left_retained:
            return 'left_retained'
        elif right_retained:
            return 'right_retained'
        else:
            return 'I'
    
    # we should never be able to get here
    assert False
    return

def refine_intron_calls( labels, stats, boundaries ):
    """Look for retained introns.
    
    """
    l = 0
    h = len( stats )
    def i_data( start_i, stop_i ):
        for label, stat, ( start, start_2 ) in \
                izip( labels[start_i:stop_i], \
                      stats[start_i:stop_i], \
                      izip( boundaries[start_i:stop_i-1], \
                            boundaries[start_i+1:stop_i] ) ):
            
            stop = start_2 - 1
            yield label, stat, ( start, stop )
    
    # first, figure out if there are any things that we can mark as 
    # definitely introns because of low signal. This will include 
    # the zero signal regions.
    for index, stat in enumerate( stats ):
        if stat.mean == 0:
            labels[ index ] = 'L'
        
    new_labels = [ labels[0], ]
    for index, (left, mid, right) in \
            enumerate( izip( i_data(0,h-2), i_data(1,h-1), i_data(2,h) ) ):
        # make sure that all of the empty regions were converted to 'L's
        assert mid[0] != 'N'
        
        # check if the intron should be retained
        if mid[0] == 'I':
            new_labels.append( check_for_retained_intron( left, mid, right ) )

        # otherwise, leave it the same
        else:            
            new_labels.append( mid[0] )

    new_labels.extend( labels[-2:] )
    
    return new_labels

def refine_one_side_retained_introns(            \
        labels, boundaries, read_coverages,      \
        one_side_retained_intron_coverage_frac,  \
        edge_bases_to_ignore ):
    
    # loop through all labels and refine one side
    # retained intron regions
    new_labels = []
    new_boundaries = []
    intervals = izip( boundaries[:-1], boundaries[1:] )
    for label, (start, stop) in izip( labels, intervals ):
        # if this is not a one side retained region, ignore it
        if label not in ('right_retained', 'left_retained'): 
            new_labels.append( label )
            new_boundaries.append( start )
        else:
            is_left_retained = ( label=='left_retained' )
            bs, ls = calc_partially_retained_intron_regions(
                start, stop, read_coverages, is_left_retained, \
                one_side_retained_intron_coverage_frac, \
                edge_bases_to_ignore )
        
            new_labels.extend( ls )
            new_boundaries.extend( bs )
    
    return new_labels, new_boundaries

def find_initial_label( prev_bt, bt, next_bt ):
    # if this is the start of an empty region this is always
    # an empty region
    if bt == 'empty':
        return 'N'
    
    # if this region is surrounded by empty regions, then it is an exon
    # TODO - give this a single exon label?
    if prev_bt == next_bt == 'empty':
        return 'S'
    
    # we know the next regions is not empty ( because we would have returned 
    # at the line above )
    if bt == 'after_empty':
        assert next_bt != 'empty'
        if next_bt == 'intron': 
            return 'Exon'
        # otherwise, if the next bndry is an exon, this is an extension
        else:
            return 'exon_extension'
    
    # If this is a proper exon i.e. this boundary is an exon_start and
    # the next boundary is an intron_start or empty region
    if bt == 'exon' and next_bt in ( 'intron', 'empty' ):
        return 'Exon'
    
    # if this is a proper intron i.e. this boundary is an intron_start 
    # and the next boundary is an exon_start or an empty region
    if bt == 'intron' and next_bt in ('exon', 'empty' ):
        return 'I'

    # otherwise, this is a possible exon extension, but should not be an
    # exon on its own
    # print >> sys.stderr, prev_bt, bt, next_bt
    return 'exon_extension'

def find_initial_boundaries_and_labels( read_cov_obj, jns ):
    boundary_types = {}
    # add junctions to boundary types
    for start, stop in jns:
        # move the stop so that boundaries correspond to region starts
        # ( ie, the start is an intron starty and the stop+1 is an 
        # exon start )
        stop += 1

        # skip junctions in empty regions
        if read_cov_obj.is_in_empty_region( start-1 ) \
                or read_cov_obj.is_in_empty_region( stop ):
            continue

        # add the start 
        boundary_types[ start ] = "intron"
        boundary_types[ stop ] = "exon"

    for start, stop in read_cov_obj.iter_empty_regions():
        boundary_types[ start ] = "empty"
        boundary_types[ stop+1 ] = "after_empty"

    # build the full set of boundaries and label them
    boundaries = sorted( boundary_types )
    labels = []
    # deal with the first label
    if boundary_types[ boundaries[0] ] == 'empty':
        labels.append( 'N' )
    elif boundary_types[ boundaries[1] ] == 'empty':
        labels.append( 'Exon' )
    else:
        labels.append( 'Exon' )

    for prev_bndry, bndry, next_bndry in izip( \
                boundaries[:-2], boundaries[1:-1], boundaries[2:] ):
        prev_bt = boundary_types[ prev_bndry ]
        bt = boundary_types[ bndry ]
        next_bt = boundary_types[ next_bndry ]

        labels.append( find_initial_label( prev_bt, bt, next_bt ) )
        # if this is an empty region start, this is always an empty region
        pass

    # deal with the last label
    if boundary_types[ boundaries[-1] ] == 'empty':
        labels.append( 'N' )
    elif boundary_types[ boundaries[-2] ] == 'empty':
        labels.append( 'Exon' )
    else:
        labels.append( 'exon_extension' )


    assert len( labels )  == len( boundaries )

    # shrink the empty regions
    if EMPTY_REGION_BNDRY_SHRINK > 0:
        for i, label in enumerate( labels[1:] ):
            if label == 'N':
                boundaries[ i+1 ] = boundaries[ i+1 ] \
                    + EMPTY_REGION_BNDRY_SHRINK
                boundaries[ i+2 ] = boundaries[ i+2 ] \
                    - EMPTY_REGION_BNDRY_SHRINK

    return labels, boundaries

def check_exon_for_gene_merge( strand, start, stop, \
                               prev_label, next_label, \
                               read_cov, cage_cov):
    # check to see if this exon merges 5' and 3' ends of a gene
    try:
        assert prev_label in ('exon_extension','I','L')
        assert next_label in ('exon_extension','I','L')
    except:
        #print >> sys.stderr, prev, curr, next
        return [ start, ], [ 'Exon', ]
    
    rca = read_cov.rca[start:(stop+1)]
    
    rca_cumsum = rca.cumsum()
    window_covs = rca_cumsum[EXON_SPLIT_SIZE:] \
        - rca_cumsum[:-EXON_SPLIT_SIZE]
    left_cov = window_covs[0]
    left_bndry = BIN_BOUNDRY_SIZE
    right_cov = window_covs[-1]
    right_bndry = len(rca)-BIN_BOUNDRY_SIZE
    global_min_cov = window_covs[BIN_BOUNDRY_SIZE:-BIN_BOUNDRY_SIZE].min() + 1

    # in the case where we have a gene followed by a single exon 
    # gene, to find the proper split ratio we can't look at the far
    # right because this is moving into a zero region
    if prev_label in ('exon_extension','I') and next_label == 'L' \
            and left_cov/global_min_cov > EXON_SPLIT_RATIO \
            and right_cov/global_min_cov <= EXON_SPLIT_RATIO:
        # moving in from the left boundary, find the location that the 
        # ratio of left to right coverage exceeds EXON_SPLIT_RATIO
        cr_low = lambda x: left_cov/(window_covs[x]+1)<EXON_SPLIT_RATIO
        first_lc_index = list(takewhile( cr_low, count(1)))[-1] + 1

        # now, find the max window to the right of this index. 
        right_cov_index = numpy.argmax(window_covs[first_lc_index:]) \
            + first_lc_index
        right_bndry = right_cov_index
        #print >> sys.stderr, right_cov, window_covs[ right_cov_index ] + 1
        right_cov = window_covs[ right_cov_index ] + 1

    # similar to above, but for an 'l' to the left
    elif prev_label == 'L' and next_label in ('exon_extension','I') \
            and right_cov/global_min_cov > EXON_SPLIT_RATIO \
            and left_cov/global_min_cov <= EXON_SPLIT_RATIO:
        # moving in from the left boundary, find the location that the 
        # ratio of left to right coverage exceeds EXON_SPLIT_RATIO
        cr_low = lambda x: right_cov/(window_covs[x]+1)<EXON_SPLIT_RATIO
        first_lc_index = list(takewhile( \
                cr_low, count(len(window_covs)-1, -1)))[-1] + 1

        # now, find the max window to the left of this index. 
        left_cov_index = numpy.argmax(window_covs[:first_lc_index])
        #print >> sys.stderr, left_cov, window_covs[ left_cov_index ] + 1
        left_cov = window_covs[ left_cov_index ] + 1
        left_bndry = left_cov_index

    if right_bndry - left_bndry < 1:
        return [ start, ], [ 'Exon', ]

    min_index = numpy.argmin( \
        window_covs[left_bndry:right_bndry+1] ) + left_bndry
    # add one to guard against divide zero, and remove low read cov
    min_cov = window_covs[ min_index ] + 1            
    
        
    # if this does not look like a gene merge exon...
    if left_cov/min_cov < EXON_SPLIT_RATIO \
            or right_cov/min_cov < EXON_SPLIT_RATIO:
        assert next != 'N'
        return [ start, ], [ 'Exon', ]
    
    new_intron_middle =  min_index + BIN_BOUNDRY_SIZE/2
    # find the new intron boundaries. use the 95% of signal rule...
    left_cumsum = rca_cumsum[:new_intron_middle]
    pos = left_cumsum.searchsorted(left_cumsum[-1]*0.999)
    new_intron_start = start + pos
    # new_intron_start = start + new_intron_middle - BIN_BOUNDRY_SIZE/2

    right_cumsum = rca_cumsum[new_intron_middle+1:] \
        - rca_cumsum[new_intron_middle]
    pos = right_cumsum.searchsorted(right_cumsum[-1]*(1-0.999))
    new_intron_stop = start + new_intron_middle + pos
    # new_intron_stop = start + new_intron_middle + BIN_BOUNDRY_SIZE/2
    
    # calc cage stats
    ds_read_cov = 0
    ds_cage_cov = 0
    if strand == '+':
        ds_read_cov = read_cov.rca[ new_intron_stop+1:stop+1 ].sum()
        ds_cage_cov = cage_cov[ new_intron_stop+1:stop+1 ].sum()
    else:
        ds_read_cov = read_cov.rca[ start:new_intron_start ].sum()
        ds_cage_cov = cage_cov[ start:new_intron_start ].sum()

    if ds_cage_cov < CAGE_MIN_SCORE \
            or (ds_cage_cov/ds_read_cov) < CAGE_TOT_FRAC:
        return [ start, ], [ 'Exon', ]

    """
    print strand, start, stop
    print prev_label, next_label, left_cov/min_cov, right_cov/min_cov
    print left_cov, right_cov, min_cov
    print 
    if start > 5070220: raw_input()
    """
    
    #print strand, start, new_intron_start, new_intron_stop+1, stop+1
    #print "%e" % ( cage_cov.sum()/read_cov.rca.sum()  )
    #print "%e" % (ds_cage_cov/ds_read_cov), ds_cage_cov, ds_read_cov
    #raw_input()

    #print >> sys.stderr, left_cov, min_cov, right_cov
    #print >> sys.stderr, start, new_intron_start, \
    #         new_intron_stop+1, stop
    assert start < new_intron_start < new_intron_stop+1 < stop
    new_bndrys = (start, new_intron_start, new_intron_stop+1)
    new_labels = ('Exon','L','Exon')
    return new_bndrys, new_labels

def refine_exon_extensions( strand, labels, bndrys, read_cov, cage_cov ):
    new_labels = [ labels[0], ]
    new_bndrys = [ bndrys[0], ]
    
    for i, (prev, curr, next) in \
            enumerate(izip( labels[:-2], labels[1:-1], labels[2:] )):
        # if it's not an exon extension continue
        if curr not in ('exon_extension', 'Exon'): 
            new_bndrys.append( bndrys[ i+1 ] )
            new_labels.append( curr )
            continue
        
        start, stop = bndrys[ i+1 ], bndrys[ i+2 ] - 1

        if curr == 'exon_extension':
            # if it's very short, and is next to a low coverage region
            # then always change it to low coverage
            # print >> sys.stderr, bndrys[ i+2 ], prev, curr, \
            #   bndrys[ i+2 ] - bndrys[ i+1 ], next
            region_len = stop - start + 1
            if region_len < BASES_TO_IGNORE_IN_RETAIN \
                    and (prev in 'L' or next in 'L'):
                new_bndrys.append( bndrys[ i+1 ] )
                new_labels.append( 'L' )
            else:
                new_bndrys.append( bndrys[ i+1 ] )
                new_labels.append( 'exon_extension' )
        elif curr == 'Exon':
            # ignore potentially single exon genes
            if prev == 'L' and next == 'L':
                new_bndrys.append( bndrys[ i+1 ] )
                new_labels.append( 'Exon' )
                continue
            
            assert not( prev == 'L' and next == 'L' )
            # if this is a short canonical exon, we don't permit a merge. 
            if stop - start + 1 < 2*BIN_BOUNDRY_SIZE + EXON_SPLIT_SIZE + 1:
                new_bndrys.append( bndrys[ i+1 ] )
                new_labels.append( 'Exon' )
                continue

            split_bndrys, split_labels, = check_exon_for_gene_merge( \
                strand, start, stop, prev, next, read_cov, cage_cov)
            new_bndrys.extend( split_bndrys )
            new_labels.extend( split_labels )

            #new_bndrys.append( bndrys[ i+1 ] )
            #new_labels.append( 'Exon' )
        else:
            assert False
    
    return new_labels, new_bndrys

def refine_single_exon_gene_label(start, stop, read_coverage_array, cage_array):
    # if it's too short
    if stop - start < MIN_SE_GENE_LEN+10:
        return 'L'

    seg_rc_array = read_coverage_array[ start:stop -1 ]
    seg_rc_array_cumsum = seg_rc_array.cumsum()
    scores = seg_rc_array_cumsum[MIN_SE_GENE_LEN:] \
                  - seg_rc_array_cumsum[:-MIN_SE_GENE_LEN]
    read_cov_score = scores.max()/MIN_SE_GENE_LEN
    # find if a region exists with sufficient read coverage. If it doesnt
    # then continue
    if read_cov_score < MIN_SE_GENE_AVG_READCOV:
        return 'L'
    
    # find the cage coverage
    cage_cov = cage_array[start:stop+1]
    cage_cov_score = score_distal_exon( cage_cov, CAGE_WINDOW_LEN  )
    if cage_cov_score < CAGE_MIN_SCORE:
        return 'L'

    """
    read_cov = seg_rc_array_cumsum[-1]
    cage_cov_sum = float( cage_cov.sum() )
    if cage_cov_sum/read_cov < CAGE_TOT_FRAC:
        return 'L'
    """
    
    return 'S'

def label_regions( strand, read_cov, bndrys, labels, cage_cov ):
    # calculate stats for each bin. 
    stats = build_bins_stats( bndrys, labels, read_cov )
    
    # change introns to retained introns when retention 
    # criterion are met
    labels = refine_intron_calls( labels, stats, bndrys )
    
    # split one side retained introns using coverage data
    labels, bndrys = refine_one_side_retained_introns( \
        labels, bndrys, read_cov, \
        ONE_SIDE_EXTENTION_FRACTION, BASES_TO_IGNORE_IN_RETAIN )
    
    # remove exon extensions that neighbor low coverage 
    # regions or exons. 
    labels, bndrys = refine_exon_extensions( \
        strand, labels, bndrys, read_cov, cage_cov )
    
    #refine_single_exon_genes( read_cov.rca, cage_cov, labels, bndrys )
    #iter_connected_regions( read_cov.read_coverage_array, labels, bndrys )

    return labels, bndrys

def find_internal_exons( exon_bndrys, intron_starts, intron_stops ):
    internal_exons = []
    for start, stop in exon_bndrys:
        if start in intron_starts: continue
        if stop in intron_stops: continue
        # make sure that both boundaries are spliced. Distal exons are dealt 
        # with seperately.
        if (not stop+1 in intron_starts) or (not start-1 in intron_stops):
            continue

        if stop == start: continue
        assert stop > start
        if stop - start < MIN_EXON_LEN: continue
        
        internal_exons.append( (start, stop) )
    
    return internal_exons

def find_T_S_exons(cov, labels, bndrys, intron_starts, intron_stops, \
                       strand, is_tss):
    assert is_tss in ( False, True )
    assert strand in '-+'
    if is_tss: distal_type_labels = ( 'TSS', 'T_S' )
    else: distal_type_labels = ( 'TES', 'T_S' )
    # find the TSS starts
    t_s_exons = []
    t_s_start_indices = [ i for i, label in enumerate(labels) \
                              if label in distal_type_labels ]
    for t_s_region_index in t_s_start_indices:
        if (is_tss and strand == '+') or ( not is_tss and strand == '-' ): 
            index_iter = xrange( t_s_region_index, len(bndrys)-1 )
        else:
            index_iter = xrange( t_s_region_index, -1, -1 )
        
        for index in index_iter:
            label = labels[ index ]
            # if we are not moving to an exon, break
            if index != t_s_region_index and \
                    label not in ('exon_extension', 'Exon'):
                break
                        
            stop = max( bndrys[ index+1 ], bndrys[ t_s_region_index+1 ] )-1
            start = min( bndrys[ index ], bndrys[ t_s_region_index ] )
            
            # print t_s_region_index, index, start, stop

            # check to see if it splices somewhere. SE genes are dealt with 
            # seperately.
            if ((is_tss and strand == '+') or ( not is_tss and strand == '-' ))\
                    and stop+1 in intron_starts:
                if start >= stop: continue
                assert start < stop
                t_s_exons.append( (start, stop) )                
            elif ((is_tss and strand == '-') or (not is_tss and strand == '+'))\
                    and start-1 in intron_stops:
                if start >= stop: continue
                assert start < stop
                t_s_exons.append( (start, stop) )
    
    # refine the eoxn bndrys
    refined_exons = []
    for exon in t_s_exons:
        exon_cvg = cov[exon[0]:exon[1]+1]
        total = exon_cvg.sum()
        if is_tss and strand == '+' or not is_tss and strand == '-':
            curr_tot = 0
            for i, val in enumerate( exon_cvg ):
                curr_tot += val
                if curr_tot > .05*total:
                    break
            refined_exons.append( (i+exon[0], exon[1]) )
        else:
            curr_tot = 0
            for i, val in enumerate( exon_cvg ):
                curr_tot += val
                if curr_tot > .95*total:
                    break
            refined_exons.append( (exon[0], exon[0]+i) )
        
        if refined_exons[-1][1] - refined_exons[-1][0] == 0:
            del refined_exons[-1]

    """
    print is_tss, strand
    print labels
    print bndrys
    print t_s_start_indices
    print t_s_exons
    print refined_exons
    #print sorted( intron_starts )
    #print sorted( intron_stops )
    raw_input()
    """

    return refined_exons

def filter_exon_bndrys( exon_bndrys, jn_bndrys ):
    """Filter out exons that start at a junction start or stop
    """
    # build the intron sets
    intron_starts = set( i[0] for i in jn_bndrys )
    intron_stops = set( i[1] for i in jn_bndrys )

    verified_bndrys = []
    for start, stop in exon_bndrys:
        if start in intron_starts: continue
        if stop in intron_stops: continue
        if stop - start + 1 > MAX_EXON_SIZE: continue
        # make sure that at least one bndry is spliced: single exon genes are 
        # dealt with seperately
        if (not stop+1 in intron_starts) and (not start-1 in intron_stops):
            continue

        verified_bndrys.append( (start, stop) )
    
    return verified_bndrys

def create_exon_gff_line( chrm, strand, start, stop, exon_id ):
    gff_line = ["chr" + chrm, 'exon_finder', 'exon', str(start), \
                    str(stop), '.', strand, '.', str(exon_id) ]
    return '\t'.join( gff_line )

def write_gff_lines( bndrys, chrm, strand, out_fp ):
    for exon_id, (start, stop) in enumerate(bndrys):
        strand_str = 'plus' if strand == '+' else 'minus'
        exon_id_str = "exon_%s_%s_%i" % ( strand_str, chrm, exon_id )
        out_fp.write(  create_exon_gff_line( \
                chrm, strand, start, stop, exon_id_str ) + "\n" )
    
    return

def find_bndry_indices( bndries, exon ):
    start_index = bndries.searchsorted( exon[0] )
    index = start_index
    while index < len( bndries ) and bndries[index] < exon[1]:
        index += 1

    return start_index, index-1

def build_cluster_labels_and_bndrys(cluster_indices, bndrys, labels, chrm_stop):
    if len( cluster_indices ) == 0:
        raise ValueError, "A cluster must have at least one non-zero region."
    
    new_bndrys = [1, bndrys[cluster_indices[0]]]
    new_labels = ['L', labels[cluster_indices[0]]]
        
    for prev_i, next_i in izip(cluster_indices[:-1], cluster_indices[1:]):
        # if there is a space between the labels, add them in
        if next_i - prev_i > 1:
            new_bndrys.append( bndrys[prev_i+1] )
            new_labels.append( 'L' )
        
        new_bndrys.append( bndrys[next_i] )
        new_labels.append( labels[next_i] )
    
    # finally, add the ending, empty label ( if we're not 
    # already at the last bndry )
    if bndrys[cluster_indices[-1]] < bndrys[-1]:
        new_bndrys.append( bndrys[cluster_indices[-1]+1] )
        new_labels.append( 'L' )
    elif bndrys[cluster_indices[-1]] < chrm_stop:
        new_bndrys.append( chrm_stop )
        new_labels.append( 'L' )
    else:
        print >> sys.stderr, "LAST BNDRY: ", \
            bndrys[cluster_indices[-1]], chrm_stop
    
    return new_bndrys, new_labels

def cluster_labels_and_bndrys( labels, bndrys, jns, strand, chrm_stop  ):
    # get exon boundaries from assigned labels and boundaries
    exon_bndrys = get_possible_exon_bndrys( \
        labels, bndrys, chrm_stop )
        
    # filter out any exon that start with a junction start or 
    # stops at an junction stop
    filtered_exon_bndrys = filter_exon_bndrys( \
        exon_bndrys, jns )
    filtered_exon_bndrys = exon_bndrys
    
    bndries_array = numpy.array( bndrys )
    exons = numpy.array( filtered_exon_bndrys )
    clusters = cluster_exons( exons, jns )
        
    clustered_labels = []
    clustered_bndrys = []
    for cluster in clusters:
        cluster_bndry_indices = set()
        for exon in cluster:
            start_index, stop_index = find_bndry_indices( bndries_array, exon )
            cluster_bndry_indices.update( xrange( start_index, stop_index+1 ) )

        cluster_bndry_indices = sorted( cluster_bndry_indices )
        
        cluster_bndrys, cluster_labels = build_cluster_labels_and_bndrys( 
            cluster_bndry_indices, bndries_array, labels, chrm_stop )
        
        clustered_labels.append( cluster_labels )
        clustered_bndrys.append( cluster_bndrys )
    
    return clustered_labels, clustered_bndrys
    
def score_distal_exon( cov, window_len, read_cov=None ):
    # get maximal tss coverage region accoss exon
    cumsum_cvg_array = \
        numpy.append(0, numpy.cumsum( cov ))

    # if the region is too short, return the total for the entire region
    # but scale so that the average is correct
    if len( cumsum_cvg_array ) <= window_len:
        return cumsum_cvg_array[-1]*( \
            float( window_len )/len( cumsum_cvg_array ) )

    # get the coverage total for every window in the interval of length
    # window_len
    scores = cumsum_cvg_array[window_len:] - cumsum_cvg_array[:-window_len]
    max_index = scores.argmax() + window_len
    score = scores[ max_index - window_len ]
    
    if NORMALIZE_BY_RNASEQ_COV and read_cov != None:
        #print read_cov
        #print max_index, max_index+window_len
        #print read_cov[max_index:(max_index+window_len)]
        #raw_input()
        #print score, read_cov[max_index:(max_index+window_len)].sum()+1, \
        #    score/(log(read_cov[max_index:(max_index+window_len)].sum()+1)+1)
        score /= (log(read_cov[max_index:(max_index+window_len)].sum()+1)+1)
    
    return score

def find_distal_exon_indices( cov, find_upstream_exons,  \
                              intron_starts, intron_stops, \
                              bs, ls, \
                              rna_seq_cov, \
                              min_score, window_len, 
                              max_score_frac ):
    # skip the first and the last labels because these are outside of the cluster
    assert ls[0] == 'L'
    if ls[-1] != 'L':
        print >> sys.stderr, ls
        print >> sys.stderr, bs
        return []
    
    # if we have signal data
    if cov != None:
        assert ls[-1] == 'L'
        scores = []
        for index, ((start, stop), label) in \
                enumerate(izip( izip(bs[1:-1], bs[2:]), ls[1:-1] )):
            if label == 'L': continue
            if rna_seq_cov != None:
                reg_rnaseq_cov = rna_seq_cov[start:stop+1]
            else:
                reg_rnaseq_cov = None
            
            scores.append( (index+1, score_distal_exon( \
                        cov[start:stop], window_len, reg_rnaseq_cov) ) )

        max_score = max( score for i, score in scores )
        if max_score > min_score:
            rvs = set( [ i for i, score in scores \
                           if float(score)/max_score > max_score_frac ] )
            #               or score >= CAGE_SUFFICIENT_SCORE] )
            
            # always add the first and/or last exon segment
            """
            if find_upstream_exons:
                rvs.add( 1 )
            else:
                rvs.add( len(bs) - 2 )
            """
            return sorted( rvs )
    
    # if we don't have signal data
    indices = []
    if True == find_upstream_exons:
        for i, (prev, cand) in enumerate(izip(ls[:-1], ls[1:])):
           if prev == 'L' and \
                   cand in ( 'Exon', 'exon_extension', 'TSS', 'T_S', 'TES' ) \
                   and bs[i+1]-1 not in intron_stops : 
               indices.append( i+1 )
    else:
        assert False == find_upstream_exons
        for i, (cand, prev) in reversed(list(enumerate(zip(ls[:-1], ls[1:])))):
           if prev == 'L' and \
                   cand in ( 'Exon', 'exon_extension', 'TSS', 'T_S', 'TES' ) \
                   and bs[i+1] not in intron_starts: 
               indices.append( i )
    
    return indices

def find_exons_in_contig( strand, read_cov_obj, jns, cage_cov, polya_cov ):
    labels, bndrys = find_initial_boundaries_and_labels( \
        read_cov_obj, jns )
    
    new_labels, new_bndrys = label_regions( \
        strand, read_cov_obj, bndrys, labels, cage_cov )
    
    clustered_labels, clustered_bndrys = cluster_labels_and_bndrys( 
        new_labels, new_bndrys, jns, strand, read_cov_obj.chrm_stop )

    intron_starts = set( jn[0] for jn in jns )
    intron_stops = set( jn[1] for jn in jns )

    all_exons = []
    all_seg_exons = []
    all_tss_exons = []
    all_internal_exons = []
    all_tes_exons = []
    
    for ls, bs in izip( clustered_labels, clustered_bndrys ):
        # search for single exon genes
        if len( ls ) == 3 and ls[1] in ( 'Exon', 'exon_extension' ):
            start, stop = bs[1], bs[2]-1
            new_label = refine_single_exon_gene_label(\
                start, stop, read_cov_obj.rca, cage_cov)
            if new_label == 'S': 
                all_seg_exons.append( (start, stop) )
            else:
                assert new_label == 'L'
            
            continue
        
        """
        # check for merege
        split_bndrys, split_labels, = check_exon_for_gene_merge( \
            strand, start, stop, prev, next, read_cov, cage_cov)
        new_bndrys.extend( split_bndrys )
        new_labels.extend( split_labels )
        """

        # find tss exons
        tss_indices = find_distal_exon_indices( \
            cage_cov, (strand=='+'), intron_starts, intron_stops,
            bs, ls, read_cov_obj.rca, \
            CAGE_MIN_SCORE, CAGE_WINDOW_LEN, CAGE_MAX_SCORE_FRAC )
        
        # if we can't find a TSS index, continue
        if len( tss_indices ) == 0:
            #print "HERE", tss_indices
            continue

        tes_indices = find_distal_exon_indices( \
            polya_cov, (strand=='-'), intron_starts, intron_stops, \
            bs, ls, None, \
            POLYA_MIN_SCORE, POLYA_WINDOW_LEN, POLYA_MAX_SCORE_FRAC )        
        # if we can't find a TSS index, continue
        if len( tes_indices ) == 0:
            #print "HERE2", tes_indices
            continue
            
        for tss_index in tss_indices:
            ls[ tss_index ] = 'TSS'
        for tes_index in tes_indices:
            if ls[ tes_index ] == 'TSS':
                ls[ tes_index ] = 'T_S'
            else:
                ls[ tes_index ] = 'TES'
        
        exon_bndrys = get_possible_exon_bndrys( \
            ls, bs, read_cov_obj.chrm_stop )
        
        internal_exons = find_internal_exons( \
            exon_bndrys, intron_starts, intron_stops )
        tss_exons = find_T_S_exons( \
            cage_cov, ls, bs, intron_starts, intron_stops, strand, True )
        tes_exons = find_T_S_exons( \
            polya_cov, ls, bs, intron_starts, intron_stops, strand, False )

        
        """
        print ls
        print bs
        print exon_bndrys
        print
        print tss_exons
        print internal_exons
        print tes_exons
        print
        #print sorted(intron_starts)
        #print sorted(intron_stops)
        raw_input()
        """
        
        all_tss_exons.extend( tss_exons )
        all_internal_exons.extend( internal_exons )
        all_tes_exons.extend( tes_exons )
        all_exons.extend( exon_bndrys )
    
    #se_genes = get_possible_single_exon_genes( \
    #    new_labels, new_bndrys, read_cov_obj.chrm_stop )
    
    return all_seg_exons, all_tss_exons, all_internal_exons, \
           all_tes_exons, all_exons
    
def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(\
        description='Find exons from wiggle and junctions files.')
    parser.add_argument( 'plus_wig', type=file, \
        help='wig file over which to search for exons on the plus strand.')
    parser.add_argument( 'minus_wig', type=file, \
        help='wig file over which to search for exons on the minus strand.')
    parser.add_argument( 'junctions', type=file, \
        help='GTF format file of junctions(introns).')
    parser.add_argument( 'chrm_sizes_fname', type=file, \
        help='File with chromosome names and sizes.')

    parser.add_argument( '--cage-wigs', type=file, nargs='+', \
        help='wig files with cage reads, to identify tss exons.')
    parser.add_argument( '--polya-reads-gffs', type=file, nargs='+', \
        help='files with polya reads, to identify tes exons.')

    
    parser.add_argument( '--out-file-prefix', '-o', default="discovered_exons",\
        help='Output file name. (default: discovered_exons)')
    parser.add_argument( '--labels-fname', '-l', type=file,\
        help='Output file name for the region labels.')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
        help='Whether or not to print status information.')
    parser.add_argument( '--make-debug-beds', '-d', default=False, \
                             action='store_true', \
        help="Produce bed files which show the intermediate label stages.")

    args = parser.parse_args()
    
    OutFPS = namedtuple( "OutFPS", [\
            "single_exon_genes", "tss_exons", \
            "internal_exons", "tes_exons", "all_exons"])
    fps = []
    for field_name in OutFPS._fields:
        fps.append(open("%s.%s.gff" % (args.out_file_prefix, field_name), "w"))
    ofps = OutFPS( *fps )
    
    # open debug filepointers if we chose that
    debug_fps = None
    if args.make_debug_beds:
        raise NotImplemented, "DEBUG output hasn't been implemented."
        debug_fps = { "+": [], "-": [] }
        for strand in "+-":
            strand_str = "plus" if strand == "+" else "minus"
            for loop in xrange( 1, 5 ):
                fname = args.out_fname + \
                    ".debug.{0}.stage{1}.bed".format(  strand_str, str( loop ))
                debug_fps[strand].append( open( fname, "w" ) )
    
    # set flag args
    global VERBOSE
    VERBOSE = args.verbose
    
    return (args.plus_wig, args.minus_wig), args.junctions, \
        args.chrm_sizes_fname, args.cage_wigs, args.polya_reads_gffs, \
       ofps, debug_fps, args.labels_fname


def main():
    (plus_wig_fp, minus_wig_fp), \
        jns_fp, chrm_sizes_fp,   \
        cage_wig_fps, tes_reads_fps, \
        out_fps, debug_fps, labels_fp = parse_arguments()
    
    # TUNING PARAMS
    empty_region_split_size = EMPTY_REGION_SPLIT_SIZE
    
    if VERBOSE: print >> sys.stderr, 'Loading read coverage data.'
    read_cov = Wiggle( chrm_sizes_fp, \
                                [plus_wig_fp, minus_wig_fp], ['+','-'] )
    
    if VERBOSE: print >> sys.stderr, 'Calculating zero intervals.'
    read_cov.calc_zero_intervals( empty_region_split_size )
    plus_wig_fp.close()
    minus_wig_fp.close()
    read_cov_sum = sum( array.sum() for array in read_cov.itervalues() )
    
    # open the cage data
    cage_cov = Wiggle( chrm_sizes_fp, cage_wig_fps )
    cage_sum = sum( cage_cov.apply( lambda a: a.sum() ).values() )
    global CAGE_TOT_FRAC
    CAGE_TOT_FRAC = (float(cage_sum)/read_cov_sum)*(1e-2)    
    for cage_fp in cage_wig_fps: cage_fp.close()
    
    # open the polya data
    polya_cov = build_coverage_wig_from_polya_reads( \
        tes_reads_fps, chrm_sizes_fp )
    polya_sum = sum( cage_cov.apply( lambda a: a.sum() ).values() )
    global POLYA_TOT_FRAC
    POLYA_TOT_FRAC = (float(polya_sum)/read_cov_sum)*(1e-2)    
    for tes_reads_fp in tes_reads_fps: tes_reads_fp.close()
    

    if VERBOSE: print >> sys.stderr,  'Loading junctions.'
    jns = parse_junctions_file_dont_freeze( jns_fp )
    
    # process each chrm, strand combination separately
    for out_fp, track_name in zip( out_fps, out_fps._fields ):
        out_fp.write( "track name=%s\n" % track_name )

    all_regions_iters = [ [], [], [], [], [] ]
    
    keys = sorted( set( jns ) )
    for chrm, strand in keys:        
        read_cov_obj = ReadCoverageData( \
            read_cov.zero_intervals[(chrm, strand)], read_cov[(chrm, strand)] )
        
        cage_cov_array = cage_cov[ (chrm, strand) ]
        polya_cov_array = polya_cov[ (chrm, strand) ]
        
        disc_grpd_exons = find_exons_in_contig( \
           strand, read_cov_obj, jns[(chrm, strand)], \
           cage_cov_array, polya_cov_array)
        
        for container, exons in zip( all_regions_iters, disc_grpd_exons ):
            regions_iter = ( GenomicInterval( chrm, strand, start, stop) \
                                 for start, stop in exons )
            container.extend( regions_iter )
    
    for regions_iter, ofp in zip( all_regions_iters, out_fps ):
        ofp.write( "\n".join(iter_gff_lines( sorted(regions_iter) )) + "\n")
    
    if debug_fps != None:
        for fp in [ fp for item in debug_fps.values() for fp in item ]:
            fp.close()
    
    return
        
if __name__ == "__main__":
    main()
