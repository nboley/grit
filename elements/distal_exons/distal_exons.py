import sys
import os
import numpy
from collections import defaultdict
from operator import itemgetter

VERBOSE = False

### Declare global default values ###
# variables to be set by script using algorithm
# the method by which each exon is scored
# should be one of [ 'total', 'avg', 'region' ]
EXON_FILTER_TYPE = 'region'
# length of region over which to calculate highest coverage within an exon
# only used if EXON_FILTER_TYPE is region
REGION_FILTER_LEN = 25
# the minimum score for an exon to be included compared 
# to the highest score at that locus/cluster
MIN_SCORE_RATIO = 0.05
# fraction of distal read coverage to extend confirmed exons out to
CVRG_FRACTION = 0.95
# Fraction of total coverage at a locus to consider the upper quartile and 
# remove exons that do not have a coverage value above this level
THRESHOLD_FRAC = 0.3


# add parent(slide) directory to sys.path and import SLIDE mods
sys.path.append( os.path.join( os.path.dirname(__file__), 
                               "..", "..", "sparsify" ) )
from gene_models import parse_junctions_file, parse_gff_line, GenomicInterval

sys.path.append( os.path.join( os.path.dirname(__file__), 
                               "..", "..", "file_types" ) )
from wiggle import Wiggle


def convert_to_genomic_intervals( all_exons ):
    converted_exons = set()
    for (chrm, strand), exons in all_exons.iteritems():
        for start, stop in exons:
            converted_exons.add( \
                GenomicInterval( chrm, strand, start, stop ) )
    
    return converted_exons

def get_exon_cvrg_scores( exons, cvrg, exon_filter_type=EXON_FILTER_TYPE, 
                          region_filter_len=REGION_FILTER_LEN ):
    """Get the scores for a group of exons
    """
    def get_exon_score( exon, cvrg ):
        # unpack exon
        start, stop = exon
        
        if stop - start + 1 <= region_filter_len or \
                exon_filter_type in [ 'total', 'avg' ]:
            # get total coverage from cvrg
            score = cvrg[start:stop+1].sum()
            
            if exon_filter_type == 'avg':
                # get average distal read coverage over exon
                score = score / (stop - start + 1)
        else:
            assert exon_filter_type == 'region', 'Invalid exon filter type.'
            # get maximal tss coverage region accoss exon
            cumsum_cvrg = numpy.append(0, numpy.cumsum( cvrg[start : stop + 1]))
            # get every coverage total over region_filter_len regions
            score = ( cumsum_cvrg[region_filter_len:]
                      - cumsum_cvrg[:-region_filter_len] ).max()
        return score
    
    
    max_score = 0
    total_score = 0
    exons_w_scores = []
    for exon in exons:
        score = get_exon_score( exon, cvrg )
        total_score += score
        if score > max_score: max_score = score
        exons_w_scores.append( (exon, score) )
    
    return exons_w_scores, max_score, total_score

def filter_exons_by_max_cvrg( exons, cvrg, exon_filter_type=EXON_FILTER_TYPE, 
                              region_filter_len=REGION_FILTER_LEN,
                              min_score_ratio=MIN_SCORE_RATIO):
    """Remove exons which contain less than MIN_SCORE_RATIO distal read coverage
    """
    def get_candidate_exons( exons_w_scores, max_score ):
        """remove exons with sufficiently low distal read coverage
        """
        # the non_inclusive minimum score is scaled to the 
        # ceiling of the max_score * min_score_ratio
        # i.e. if the ratio is 0.01 and the max score is 30
        # a score greater than one will suffice for *high score*
        min_score = int( min_score_ratio * max_score ) + 1
        
        high_score_exons = []
        for exon, score in exons_w_scores:
            if score > min_score:
                high_score_exons.append( exon )
        
        return high_score_exons
    
    
    exons_w_scores, max_score, total_score = \
        get_exon_cvrg_scores( exons, cvrg, exon_filter_type, region_filter_len )
    candidate_exons = get_candidate_exons( exons_w_scores, max_score )
    
    return candidate_exons

def get_exons_w_quantls( exons, all_cvrg ):
    """Get the highest coverage position within each exon
    1) collapse exons so as not to get duplicate coverage from one bp
    2) retvrieve cvrg overlapping exons from this gene
    3) sort cvrg and maintain list of original positions
    4) get the maximal quantile at which each exon would be included as a 
       valid tss exon
    """
    def get_collapsed_regions( exons ):
        """Get non-overlapping exon boundaries from which to get coverage
        note exons must be previously sorted
        """
        regions = []
        
        curr_start = exons[0][0]
        curr_stop = exons[0][1]
        
        # for converting between array position and exon position
        region_lengths = []
        for start, stop in exons[1:]:
            # if this exon does not overlap any previous exons
            if start > curr_stop + 1:
                # add the current collapsed exon and move the curr_start
                regions.append( (curr_start, curr_stop) )
                region_lengths.append( curr_stop - curr_start + 1 )
                curr_start = start
            
            curr_stop = max( stop, curr_stop )
        
        # add the last collapsed exon
        regions.append( (curr_start, curr_stop) )
        region_lengths.append( curr_stop - curr_start + 1 )
        regions_cumsum = numpy.append(
            0, numpy.cumsum( region_lengths ) )
        
        return regions, regions_cumsum
    
    def get_cvrg_structs( all_cvrg, regions ):
        cvrg_arrays = []
        for start, stop in regions:
            cvrg_arrays.append( \
                all_cvrg[ start:stop+1 ] )
        
        # flatten coverage array for locus
        gene_cvrg = numpy.hstack( cvrg_arrays )
        if gene_cvrg.sum() == 0:
            raise ValueError, 'No coverage in this locus: ' + \
                '{0:d}-{0:d}'.format( regions[0][0], regions[-1][1] )
        
        # the positions along the coverage array sorted by tss cvrg values
        # at those positions
        sorted_pos_s = gene_cvrg.argsort()[::-1]
        # the fraction of coverage up to this sorted coverage position
        cum_cvrg_frac = numpy.cumsum(
            numpy.sort( gene_cvrg ) / gene_cvrg.sum() )[::-1]
        
        return sorted_pos_s, cum_cvrg_frac
    
    def get_exon_quantl( start, stop, region_starts, regions_cumsum, 
                         sorted_pos_s, cum_cvrg_frac ):
        """Get the quantile of coverage along the cum_cvrg_frac at which this
        exon will be included as a valid tss exon
        
        Algorithm:
        1) Convert the genomic start stop positions to positions relative to
           the coverage array for this gene
        2) Identify the index of the first base from this exon is the 
           sorted_pos_s array. The index of this position cooresponds to
           its quantile in the cum_cvrg_frac array
        """
        def get_array_bndrys( start, stop ):
            # get the index of the containing region
            region_i = region_starts.searchsorted( start, side='right' ) - 1
            region_start = region_starts[ region_i ]
            bases_into_region = start - region_start
            array_start = regions_cumsum[ region_i ] + bases_into_region
            array_stop = array_start + (stop - start)
            
            return array_start, array_stop
        
        
        # get position relative to the coverage array
        array_start, array_stop = get_array_bndrys( start, stop )
        # find the first value in the sorted positions array that is a position 
        # within this exon and get the index of that position for indexing 
        # into the cum_cvrg_frac which returns the quantile at which this exon
        # will be introduced as a valid tss exon
        index = numpy.where( (sorted_pos_s >= array_start) &
                             (sorted_pos_s <= array_stop) 
                             )[0][0]
        
        return cum_cvrg_frac[ index ]
    
    
    collapsed_regions, regions_cumsum = get_collapsed_regions( exons )
    sorted_pos_s, cum_cvrg_frac = get_cvrg_structs(
        all_cvrg, collapsed_regions )
    
    exons_w_quantls = []
    # get region exon starts in order to find array coords from genomic coords
    region_starts = numpy.array( [ i[0] for i in collapsed_regions ] )
    for start, stop in exons:
        # get the maximum quantile cutoff at which this exon would be
        # introduced as a valid tss exon
        quantl = get_exon_quantl( start, stop, region_starts, regions_cumsum,
                                  sorted_pos_s, cum_cvrg_frac )
        
        exons_w_quantls.append( ((start, stop), quantl) )
    
    return exons_w_quantls

def filter_exons_by_quantiles( exons, all_cvrg, thresh_frac=THRESHOLD_FRAC ):
    """
    0) smooth cage signal if requested
    1) order bp position coverage values
    2) determine the highest coverage pos in each exon (get_exons_w_quantls)
    3) keep exons that stil have signal above thresh_frac
    """
    try:
        exons_w_quantls = get_exons_w_quantls( exons, all_cvrg )
    except ValueError:
        # print 'Warning no coverage at a locus.'
        return []
    
    upper_quantile_exons = []
    for exon, quantl in exons_w_quantls:
        if quantl > thresh_frac:
            upper_quantile_exons.append( exon )
    
    return upper_quantile_exons

def get_unique_internal_splice_exons( exons, strand, is_tss ):
    # get index of internal and external exons
    if (strand == '+' and is_tss) or (strand == '-' and not is_tss):
        int_index = 1
    else:
        int_index = 0
    ext_index = 0 if int_index == 1 else 1
    
    unique_largest_exons = {}
    # get the largest exon at each internal splice site
    for exon in exons:
        int_coord = exon[int_index]
        ext_coord = exon[ext_index]
        if int_coord in unique_largest_exons:
            if (ext_coord < unique_largest_exons[int_coord] \
                    and int_index == 1) \
                    or (ext_coord > unique_largest_exons[int_coord] \
                            and int_index == 0):
                unique_largest_exons[ int_coord ] = ext_coord
        else:
            unique_largest_exons[ int_coord ] = ext_coord
    
    # get exon tuple structure back
    unique_exons = []
    for int_coord, ext_coord in unique_largest_exons.iteritems():
        unique_exons.append( (min(int_coord, ext_coord), \
                                  max(int_coord, ext_coord) ) )
    
    sorted_unique_exons = sorted( unique_exons )
    
    return sorted_unique_exons

def refine_exon_bndrys( exons, all_cvrg, key, is_tss, 
                        cvrg_fraction=CVRG_FRACTION ):
    """Refine distal exons according to distal reads coverage across each exon
    """
    chrm, strand = key
    refined_exons = []
    for start, stop in exons:
        from_upstream = (strand == '+' and not is_tss) \
            or (strand == '-' and is_tss)
        refined_exon = all_cvrg.get_region_w_fraction_cvrg( \
            key, start, stop, cvrg_fraction, from_upstream )
        if None == refined_exon:
            continue
        refined_exons.append( refined_exon )
    
    return refined_exons

def find_distal_exons( clustered_exons, all_cvrg, chrm, strand, \
                           is_tss, cvrg_fraction ):
    """find exons which have distal read coverage and refine their boundaries
    """
    all_exons = []
    for exon_cluster_id, exons in enumerate(clustered_exons):
        # uniquify by internal coordinate
        unique_exons = get_unique_internal_splice_exons( \
            exons, strand, is_tss )
        
        # refine external boundary using distal read signal
        refined_exons = refine_exon_bndrys( \
            unique_exons, all_cvrg, (chrm, strand), is_tss, cvrg_fraction )
        
        all_exons.extend( refined_exons )
    
    return all_exons

def find_all_distal_exons( clustered_exons, all_cvrg, is_tss=True,
                           cvrg_fraction=CVRG_FRACTION ):
    """wrapper for find_distal_exons
    """
    all_exons = {}
    # process each chrm, strand combination separately
    keys = sorted( set( clustered_exons ).intersection( all_cvrg ) )
    for (chrm, strand) in keys:
        exons = find_distal_exons(
            clustered_exons[(chrm, strand)], all_cvrg, \
            chrm, strand, is_tss, cvrg_fraction)
        all_exons[(chrm, strand)] = exons
    
    return all_exons

if __name__ == "__main__":
    print "Please use discover_tes_exons_from_polya_reads.py or " + \
        "discover_tss_exons_from_cage_data.py."
