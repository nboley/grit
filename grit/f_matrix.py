"""
Copyright (c) 2011-2015 Nathan Boley

This file is part of GRIT.

GRIT is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GRIT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GRIT.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
sys.setrecursionlimit(10000)

import numpy
from scipy.spatial import KDTree

MIN_NUM_MAPPABLE_BASES = 1
LET_READS_OVERLAP = True

DEBUG=False

from scipy.stats import beta

import config

import networkx as nx

import frag_len

from itertools import product, izip, chain
from collections import defaultdict

from grit.files.reads import ( iter_coverage_intervals_for_read, get_read_group,
                               CAGEReads, RAMPAGEReads, PolyAReads )

class NoObservableTranscriptsError(Exception):
    pass

################################################################################
#
#
# Code for converting transcripts into 
#
#

def find_nonoverlapping_contig_indices( contigs, contig_boundaries ):
    """Convert transcripts into lists of non-overlapping contig indices.
    
    """
    non_overlapping_contig_indices = []
    for contig in contigs:
        assert contig[0] in contig_boundaries \
            and contig[1]+1 in contig_boundaries

        start_i = contig_boundaries.searchsorted( contig[0] )
        assert contig[0] == contig_boundaries[start_i]

        stop_i = start_i + 1
        while contig[1] > contig_boundaries[stop_i]-1:
            stop_i += 1            
        assert contig[1] == contig_boundaries[stop_i]-1

        non_overlapping_contig_indices.extend( xrange(start_i,stop_i) )
    
    return non_overlapping_contig_indices

def build_nonoverlapping_indices( transcripts, exon_boundaries ):
    # build the transcript composed of pseudo ( non overlapping ) exons. 
    # This means splitting the overlapping parts into new 'exons'
    for transcript in transcripts:
        yield find_nonoverlapping_contig_indices( 
            transcript.exons, exon_boundaries )
    
    return

################################################################################

def simulate_reads_from_exons( n, fl_dist, \
                               read_len, max_num_unmappable_bases, \
                               exon_lengths ):
    """Set read_len == None to get the full fragment.

    """
    
    import random
    exon_lengths = numpy.array( exon_lengths )

    exon_lens_cum = list( exon_lengths.cumsum() )
    exon_lens_cum.insert( 0, 0 )
    exon_lens_cum = numpy.array( exon_lens_cum )
    # 1111, 1112, 1122, 1222, 2222
    bin_counts = {}
    
    def simulate_read( ):
        # until we shouldnt choose another read
        while True:
            # choose a fragmnent length
            rv = random.random()
            fl_index = fl_dist.fl_density_cumsum.searchsorted( rv, side='left' )
            assert fl_index < len(fl_dist.fl_density_cumsum )
            fl = fl_index + fl_dist.fl_min
                        
            # make sure the fl is possible given the reads
            if fl > exon_lens_cum[-1]:
                continue
            if read_len != None \
                and not LET_READS_OVERLAP \
                and fl < 2*read_len:
                continue
            
            def get_bin( position ):
                # find the start exon for the first read
                return exon_lens_cum.searchsorted( position, side='right' )-1
            
            def get_bin_pairs( read_start, read_len ):
                ## The following diagram makes is much easier to 
                ## trace the subsequent code.
                #     5    4    Exon lens 5 and 4 
                # XXXXX|XXXX    Bndry at |
                #    RR         Read len 2
                #    34    ( 0 indexed ) 
                #     R R       Read Len 2
                #     4 5  ( 0 Indexed )
                #       RR      Read Len 2
                #       56 ( 0 Indexed )
                start_exon = get_bin( read_start )
                
                # find what bin these correspond to. For basically, we just need
                # to check if this crosses a junction.
                loop = 0
                while read_start + read_len > exon_lens_cum[start_exon + loop]:
                    loop += 1
                return tuple(range( start_exon, start_exon + loop ))
            
            # choose the start location
            start_location = random.randrange( 0, exon_lens_cum[-1] - fl + 1 )
            end_location = start_location + fl - 1
            
            # if we want to full fragment
            if read_len == None:
                bin1 = get_bin( start_location )
                bin2 = get_bin( end_location )
                assert bin1 != None and bin2 != None
                return tuple(xrange(bin1, bin2+1))
            # if we want just the pairs
            else:
                pair_1 = get_bin_pairs( start_location, read_len )
                if pair_1 == None: continue
                pair_2 = get_bin_pairs( end_location - read_len + 1, read_len )
                if pair_2 == None: continue
                return ( pair_1, pair_2 )            
    
    for loop in xrange( n ):
        bin = simulate_read()
        try:
            bin_counts[ bin ] += 1
        except KeyError:
            bin_counts[ bin ] = 1

    return bin_counts


def calc_long_read_expected_cnt_given_bin( bin, exon_lens, fl_dist ):
    density = 0
    exon_lens = [ exon_lens[exon] for exon in bin ]
    if len( exon_lens ) == 1: 
        # if no fragments fit within this, return 0
        if fl_dist.fl_min > exon_lens[0]:
            return 0.0
        
        # if the fragment is longer than the exon, then it can't have come 
        # from this exon, so set the maximum fragment length to the exon 
        # minimum of the exon length and the global max frag length
        upper_fl_bnd = min( fl_dist.fl_max, exon_lens[0] )
        
        # this is all just algebra on 'faster way'
        density += ( exon_lens[0] + 1 )*fl_dist.fl_density_cumsum[ 
            upper_fl_bnd - fl_dist.fl_min ]
        density -= fl_dist.fl_density_weighted_cumsum[ 
            upper_fl_bnd - fl_dist.fl_min ]
        
        return density
    
    gap = sum( exon_lens[1:-1] )
    for x in xrange( 0, exon_lens[0] ):
        # below, we use an inner loop to look at all possible fragment lengths.
        # of course, this isnt even close to necessary because we have 
        # precalcualted the cumsum. So, we really only need the minimum 
        # fragment length, which is just max( exon_lens[0] - x + gap + 0, 
        # fl_dist.fl_min ) and the maximum fragment length, 
        # min( exon_lens[0] - x + gap + exon_lens[-1], fl_dist.fl_max )
        # slow ( but correct ) old code
        """
        for y in xrange( 0, exon_lens[-1] + 1 ):
            f_size = exon_lens[0] - x + gap + y
            # skip fragments that are too long or too short
            if f_size < fl_dist.fl_min: continue
            if f_size > fl_dist.fl_max: continue
            density_1 += fl_dist.fl_density[ f_size - fl_dist.fl_min ]
        """

        min_fl = max( exon_lens[0] - x + gap + 1, fl_dist.fl_min )
        if min_fl > fl_dist.fl_max: continue
        max_fl = min( exon_lens[0] - x + gap + exon_lens[-1], fl_dist.fl_max )
        if max_fl < fl_dist.fl_min: continue
        assert min_fl <= max_fl
        density += fl_dist.fl_density_cumsum[ max_fl - fl_dist.fl_min ]
        if min_fl > fl_dist.fl_min:
            density -= fl_dist.fl_density_cumsum[ min_fl - fl_dist.fl_min - 1 ]
    
    return density

def find_short_read_bins_from_fragment_bin(
        bin, exon_lens, read_len, min_num_mappable_bases ):
    # one of the end positions of a fragment must be
    # in the short read's bin, so the question is how many
    # junctions each end are able to cross. The general rule is 
    # that the total length of the exons in the short read bin 
    # must be greater than the read length, but the length of
    # the internal exons ( ie, in 1,2,3 the internal exon is 2 )
    # must be shorter than the read - 2*min_num_mappable_bases ( 
    # to account for the mappability )
    def find_possible_end_bins( ordered_exons ):
        for end_index in xrange( 1, len( ordered_exons )+1 ):
            # get the new bin as the first end_index exons
            new_bin  = ordered_exons[0:end_index]
            # find the length of the internal exons. If they are too long,
            # then any new exons will make them even longer, so we are done
            internal_exons_len = sum( exon_lens[i] for i in new_bin[1:-1] )
            if internal_exons_len - 2*min_num_mappable_bases > read_len:
                return
            
            # calculate all of the exon lens. If they are too short, then
            # continue ( to add more exons ). We can't just add the beginning
            # and final exons to internal exons because of the single exon case
            # ( although we can get around this, but, premature optimization...)
            if sum( exon_lens[i] for i in new_bin ) < read_len:
                continue

            yield tuple(new_bin)
        
        return
    
    # first, do the 5' end
    starts = []
    for end in find_possible_end_bins( bin ):
        starts.append( end )
    
    # next, the 3' end, using the reversed bin
    stops = []
    for end in find_possible_end_bins( bin[::-1] ):
        stops.append( end[::-1] )
    
    # join all of the starts and stops
    paired_read_bins = set()
    for start in starts:
        for stop in stops:
            paired_read_bins.add( (start, stop) )
    
    single_end_read_bins = set( starts )
    single_end_read_bins.update( stops )
    
    return single_end_read_bins, paired_read_bins

def find_possible_read_bins_for_transcript( 
    trans_exon_indices, exon_lens, fl_dist, read_len, min_num_mappable_bases=1):
    """Find all possible bins. 

    """
    full_fragment_bins = set()
    paired_read_bins = {}
    single_end_read_bins = {}
    
    transcript_exon_lens = numpy.array(
        [exon_lens[i] for i in trans_exon_indices])
    transcript_exon_lens_cumsum = transcript_exon_lens.cumsum()
    transcript_exon_lens_cumsum = numpy.insert(transcript_exon_lens_cumsum,0,0)
    for index_1 in xrange(len(transcript_exon_lens)):
        for index_2 in xrange(index_1, len(transcript_exon_lens)):
            # make sure the bin is allowed by the fl dist
            if index_2 - index_1 > 2:
                min_frag_len = \
                    transcript_exon_lens_cumsum[index_2-1+1]   \
                    - transcript_exon_lens_cumsum[index_1+1] \
                    + 2*min_num_mappable_bases
            else:
                min_frag_len = 0
            
            max_frag_len = \
                transcript_exon_lens_cumsum[index_2+1] \
                - transcript_exon_lens_cumsum[index_1]
            
            if max_frag_len < fl_dist.fl_min or min_frag_len > fl_dist.fl_max:
                continue
            
            bin = tuple(trans_exon_indices[index_1:index_2+1])
        
            full_fragment_bins.add( bin )
            single, paired = find_short_read_bins_from_fragment_bin( 
                bin, exon_lens, read_len, min_num_mappable_bases )
        
            paired_read_bins[bin] = paired
            single_end_read_bins[bin] =  single
    
    return full_fragment_bins, paired_read_bins, single_end_read_bins

def iter_possible_fragment_bins_for_splice_graph( 
    graph, segment_bnds, fl_dist ):
    """Find all possible fragment bins. 

    """
    def seg_len(i):
        return segment_bnds[i+1] - segment_bnds[i] 
    
    def min_path_len(path):
        return sum(seg_len(i) for i in path[1:-1]) + min(2, len(path))

    def max_path_len(path):
        return sum(seg_len(i) for i in path)

    paths = [[i,] for i in graph.nodes()]
    while len(paths) > 0:
        curr_path = paths.pop()
        if max_path_len(curr_path) >= fl_dist.fl_min:
            yield curr_path
        for child in graph.successors(curr_path[-1]):
            child_path = curr_path + [child,]
            if min_path_len(child_path) <= fl_dist.fl_max:
                paths.append( child_path )
    
    return

def estimate_num_paired_reads_from_bin( 
        bin, transcript, exon_lens,
        fl_dist, read_len, min_num_mappable_bases=1 ):
    """

    """
    assert min_num_mappable_bases > 0, \
        "It doesn't make sense to map a read into a segment with 0 bases"
    #global DEBUG
    #if transcript == (20, 21):
    #    print transcript, bin, fl_dist
    #    DEBUG = True
    
    # calculate the exon lens for the first and second reads
    fr_exon_lens = [ exon_lens[i] for i in bin[0]  ]
    sr_exon_lens = [ exon_lens[i] for i in bin[1]  ]
    # find the length of the internal exons ( those not in either bin )
    pre_sr_exon_lens = sum( exon_lens[i] for i in transcript \
                            if i < bin[1][0] and i >= bin[0][0] )

    if DEBUG:
        print "Bin", bin
        print "FR Exon Lens", fr_exon_lens
        print "SR Exon Lens", sr_exon_lens
        print "Pre FR Exon Lens", pre_sr_exon_lens
    
    # loop through the possible start indices for the first bin to be satisifed
    # basically, after removing the spliced introns, we need at least min_num_ma
    # in the first exon, and in the last exon. This puts a lower bound on the 
    # read start at position ( relative to the transcript ) of 0 ( the read can 
    # start at the first position in exon 1 ) *or* last_exon_start + 1-read_len 
    min_start = max( 0, sum( fr_exon_lens[:-1] ) 
                         + min_num_mappable_bases - read_len )
    # make sure that there is enough room for the second read to start in 
    # sr_exon[0] given the fragment length constraints
    min_start = max( min_start, pre_sr_exon_lens + read_len - fl_dist.fl_max )

    max_start = fr_exon_lens[0] - min_num_mappable_bases
    # make sure that there are enough bases in the last exon for the read to fit
    max_start = min( max_start, sum( fr_exon_lens ) - read_len )
    
    if DEBUG:
        print "Start Bnds", min_start, max_start
    
    # find the range of stop indices.
    # first, the minimum stop is always at least the first base of the first 
    # exon in the second read bin
    min_stop = pre_sr_exon_lens + read_len
    # second, we know that at least min_num_mappable_bases are in the last exon
    # of the second read, because we assume that this read is possible. 
    min_stop = max( min_stop, pre_sr_exon_lens \
                              + sum( sr_exon_lens[:-1] ) \
                              + min_num_mappable_bases )

    # max stop can never be greater than the transcript length
    max_stop = pre_sr_exon_lens + sum( sr_exon_lens )
    # all but min_num_mappable_bases basepair is in the first exon of second rd
    # again, we assume that at least min_num_mappable_bases bases are in the
    # last exon because the read is possible. However, it could be that this 
    # extends *past* the last exon, which we account for on the line after.x
    max_stop = min( max_stop, \
                    pre_sr_exon_lens + sr_exon_lens[0] \
                    - min_num_mappable_bases + read_len )
    # all basepairs of the last exon in the second read are in the fragment. 
    # Make sure the previous calculation doesnt extend past the last exon.
    max_stop = min( max_stop, 
                    min_stop - min_num_mappable_bases + sr_exon_lens[-1] )
    
    # ensure that the min_stop is at least 1 read length from the min start
    min_stop = max( min_stop, min_start + read_len )

    if DEBUG:
        print "Stop Bnds", min_stop, max_stop
    def do():
        density = 0.0
        for start_pos in xrange( min_start, max_start+1 ):
            min_fl = max( min_stop - start_pos, fl_dist.fl_min )
            max_fl = min( max_stop - start_pos, fl_dist.fl_max )
            if min_fl > max_fl: continue
            density += fl_dist.fl_density_cumsum[ max_fl - fl_dist.fl_min ]
            if min_fl > fl_dist.fl_min:
                density -= fl_dist.fl_density_cumsum[min_fl - fl_dist.fl_min-1]
        return density
    
    density = do()
    if DEBUG:
        print "Density", density
        print
    
    #DEBUG = False
    return float( density )

def calc_expected_cnts( exon_boundaries, transcripts, fl_dist, 
                        r1_len, r2_len,
                        max_num_unmappable_bases=MIN_NUM_MAPPABLE_BASES,
                        max_memory_usage=3.5 ):
    assert r1_len == r2_len, "Paired reads must have the same lengths"
    read_len = r1_len
    
    # store all counts, and count vectors. Indexed by ( 
    # read_group, read_len, bin )
    cached_f_mat_entries = {}
    f_mat_entries = {}
    
    nonoverlapping_exon_lens = \
        numpy.array([ stop - start for start, stop in 
                      izip(exon_boundaries[:-1], exon_boundaries[1:])])
    
    # for each candidate trasncript
    for transcript_index, nonoverlapping_indices in enumerate(transcripts):
        if (len(exon_boundaries)*len(transcripts)*8.)/(1024**3) > max_memory_usage:
            raise MemoryError, \
                "Building the design matrix has exceeded the maximum allowed memory "
        nonoverlapping_indices = tuple(nonoverlapping_indices)
        f_mat_entries[nonoverlapping_indices] = {}
        # find all of the possible read bins for transcript given this 
        # fl_dist and read length
        ( full, pair, single 
            ) = find_possible_read_bins_for_transcript( 
            nonoverlapping_indices, nonoverlapping_exon_lens, fl_dist, 
            read_len, min_num_mappable_bases=1 )
        
        # add the expected counts for paired reads
        for full_bin, paired_bins in pair.iteritems():
            for bin in paired_bins:
                # we can only re-use cached full_bin/bin combos
                # because it's possible for the middle of a fragment
                # to skip a region in one transcript, and be spliced
                # out in another transcript
                key = (full_bin, bin)
                if key in cached_f_mat_entries:
                    pseudo_cnt = cached_f_mat_entries[ key ]
                else:
                    pseudo_cnt = estimate_num_paired_reads_from_bin(
                        bin, full_bin, 
                        nonoverlapping_exon_lens, fl_dist,
                        read_len, max_num_unmappable_bases )

                    cached_f_mat_entries[ key ] = pseudo_cnt
                if pseudo_cnt > 0:
                    f_mat_entries[nonoverlapping_indices][bin] = pseudo_cnt
    
    return f_mat_entries

def convert_f_matrices_into_arrays( f_mats, normalize=True ):
    expected_cnts = []
    observed_cnts = []
    zero_entries = []
    
    for key, (expected, observed) in f_mats.iteritems():
        expected_cnts.append( expected )
        observed_cnts.append( observed )
    
    # normalize the expected counts to be fraction of the highest
    # read depth isoform
    expected_cnts = numpy.array( expected_cnts )
    
    zero_entries =  ( expected_cnts.sum(0) == 0 ).nonzero()[0]
    
    if zero_entries.shape[0] > 0:
        expected_cnts = expected_cnts[:, expected_cnts.sum(0) > 0 ]
    
    if normalize:
        assert expected_cnts.sum(0).min() > 0
        expected_cnts = expected_cnts/expected_cnts.sum(0)
    
    observed_cnts = numpy.array( observed_cnts )
    
    assert expected_cnts.sum(0).min() > 0
    assert  expected_cnts.sum(1).min() > 0
    
    return expected_cnts, observed_cnts, zero_entries.tolist()

def build_observed_cnts( binned_reads, fl_dists ):
    rv = {}
    for ( read_len, read_group, bin ), value in binned_reads.iteritems():
        rv[ ( read_len, read_group, bin ) ] = value
    
    return rv

def build_expected_and_observed_arrays( 
        expected_cnts, observed_cnts, normalize=True ):
    expected_mat = []
    observed_mat = []
    unobservable_transcripts = set()
    
    # find and sort all transcripts
    transcripts = set()
    for bin, transcript_cnts in sorted(expected_cnts.iteritems()):
        for transcript in transcript_cnts:
            transcripts.add(transcript)
    transcripts = sorted(transcripts)

    # turn the cnts into a more structured format
    for bin, transcript_cnts in sorted(expected_cnts.iteritems()):
        expected_cnts = [transcript_cnts[t] for t in transcripts]
        expected_mat.append( expected_cnts )
        try:
            observed_mat.append( observed_cnts[bin] )
        except KeyError:
            observed_mat.append( 0 )
    
    if len( expected_mat ) == 0:
        raise ValueError, "No expected reads."
    
    expected_mat = numpy.array( expected_mat, dtype=numpy.double )
    if normalize:
        nonzero_entries = expected_mat.sum(0).nonzero()[0]
        unobservable_transcripts = set(range(expected_mat.shape[1])) \
            - set(nonzero_entries.tolist())
        expected_mat = expected_mat/(expected_mat.sum(0)+1e-12)
    
    observed_mat = numpy.array( observed_mat, dtype=numpy.int )
    return expected_mat, observed_mat, unobservable_transcripts

def build_expected_and_observed_rnaseq_counts(gene, reads, fl_dists):
    # find the set of non-overlapping exons, and convert the transcripts to 
    # lists of these non-overlapping indices. All of the f_matrix code uses
    # this representation.     
    exon_boundaries = numpy.array(gene.find_nonoverlapping_boundaries())
    transcripts_non_overlapping_exon_indices =list(build_nonoverlapping_indices(
            gene.transcripts, exon_boundaries ))
    
    binned_reads = bin_rnaseq_reads( 
        reads, gene.chrm, gene.strand, exon_boundaries)
    observed_cnts = build_observed_cnts( binned_reads, fl_dists )    
    read_groups_and_read_lens =  set( (RG, read_len) for RG, read_len, bin 
                                        in binned_reads.iterkeys() )
    """
    # TODO XXX
    print exon_boundaries
    print observed_cnts.values()
    """
    expected_cnts = defaultdict(lambda: defaultdict(float))
    for (rg, (r1_len,r2_len)), (fl_dist, marginal_frac) in fl_dists.iteritems():
        for transcript, read_bins_and_vals in calc_expected_cnts( 
                exon_boundaries, transcripts_non_overlapping_exon_indices, 
                fl_dist, r1_len, r2_len).iteritems():
            for read_bin, expected_bin_cnt in read_bins_and_vals.iteritems():
                assert r1_len == r2_len
                expected_cnts[(r1_len, rg, read_bin)][transcript] += (
                    marginal_frac*expected_bin_cnt )

    return dict(expected_cnts), observed_cnts

def build_expected_and_observed_transcript_bndry_counts( 
        gene, reads, bndry_type=None ):
    if bndry_type == None:
        # try to infer the boundary type from the reads type
        if reads.type == 'CAGE': bndry_type = "five_prime"
        elif reads.type == 'PolyA': bndry_type = "three_prime"
        elif reads.type == 'RAMPAGE': bndry_type = "five_prime"
        else: assert False, "Unsupported boundary read type (%s)" % reads.type
    assert bndry_type in ["five_prime", "three_prime"]
    
    cvg_array = numpy.zeros(gene.stop-gene.start+1)
    cvg_array += reads.build_read_coverage_array( 
        gene.chrm, gene.strand, gene.start, gene.stop )
    
    # find all the bndry peak regions
    peaks = list()
    for transcript in gene.transcripts:
        if bndry_type == 'five_prime':
            peaks.append( transcript.find_promoter() )
        elif bndry_type == 'three_prime':
            peaks.append( transcript.find_polya_region() )
        else: assert False
    
    peak_boundaries = set()
    for start, stop in peaks:
        peak_boundaries.add( start )
        peak_boundaries.add( stop + 1 )
    peak_boundaries = numpy.array( sorted( peak_boundaries ) )
    pseudo_peaks = zip(peak_boundaries[:-1], peak_boundaries[1:])
    
    # build the design matrix. XXX FIXME
    expected_cnts = defaultdict( lambda: [0.]*len(gene.transcripts) )
    for transcript_i, peak in enumerate(peaks):
        nonoverlapping_indices = \
            find_nonoverlapping_contig_indices( 
                [peak,], peak_boundaries )
        # calculate the count probabilities, adding a fudge to deal with 0
        # frequency bins
        tag_cnt = cvg_array[
            peak[0]-gene.start:peak[1]-gene.start].sum() \
            + 1e-6*len(nonoverlapping_indices)
        for i in nonoverlapping_indices:
            ps_peak = pseudo_peaks[i]
            ps_tag_cnt = cvg_array[
                ps_peak[0]-gene.start:ps_peak[1]-gene.start].sum()
            expected_cnts[ ps_peak ][transcript_i] \
                = (ps_tag_cnt+1e-6)/tag_cnt
        
    # count the reads in each non-overlaping peak
    observed_cnts = {}
    for (start, stop) in expected_cnts.keys():
        observed_cnts[ (start, stop) ] \
            = int(round(cvg_array[start-gene.start:stop-gene.start].sum()))
    
    return expected_cnts, observed_cnts

def cluster_bins(expected_rnaseq_cnts):
    if config.DEBUG_VERBOSE:
        config.log_statement( "Normalizing bin frequencies" )
    
    clustered_bins = defaultdict(list)
    for bin, transcripts_and_cnts in expected_rnaseq_cnts.items():
        row = numpy.array([x[1] 
                           for x in sorted(transcripts_and_cnts.iteritems())])
        key = tuple((100000*row/row.sum()).round().tolist())
        clustered_bins[key].append(bin)
    return clustered_bins.values()


def cluster_rows(expected_rnaseq_array, observed_rnaseq_array):
    if config.DEBUG_VERBOSE:
        config.log_statement( "Normalizing bin frequencies" )
    
    clustered_rows = defaultdict(list)
    for i, row in enumerate(expected_rnaseq_array):
        key = tuple((100000*row/row.sum()).round().tolist())
        clustered_rows[key].append(i)
    clusters = clustered_rows.values()
    
    """ This worked quickly, but could use way too much memory
    norm_rows = []
    for i, row in enumerate(expected_rnaseq_array):
        norm_rows.append(row/row.sum())
        
    if config.DEBUG_VERBOSE:
        config.log_statement( "Building KDTree to cluster bins" )
    tree = KDTree(numpy.vstack(norm_rows), leafsize=30)
    
    edges = set()
    points = set(xrange(len(norm_rows)))
    while len(points) > 0:
        i = points.pop()
        neighbors = tree.query_ball_point(norm_rows[i], 1e-6, 1)
        for j in neighbors:
            edges.add((i, j))
            points.discard(j)
        if config.DEBUG_VERBOSE and len(points)%10 == 0:
            config.log_statement("%i rows remain to be clustered" % len(points))

    graph = nx.Graph()
    graph.add_nodes_from(xrange(expected_rnaseq_array.shape[0]))
    graph.add_edges_from(edges)
    clusters = nx.connected_components(graph)
    print len(clusters)
    print sorted(clusters)
    """
    
    """
    Use the tree to find pairs. We dont use this because we only care about distance.
    #for i, row in enumerate(norm_rows):
    #    config.log_statement("%i/%i-%s" % (i, len(norm_rows), tree.query_ball_point(row, 1e-6, 1)))
    pairs = tree.query_pairs(1e-6, 1)

    if config.DEBUG_VERBOSE:
        config.log_statement( "Finding clusters")    
    graph = nx.Graph()
    graph.add_nodes_from(xrange(expected_rnaseq_array.shape[0]))
    graph.add_edges_from(pairs)
    clusters = nx.connected_components(graph)
    """
    
    """
    Brute force code for finding the matching points. This works, but is slow.
    norm_rows = []
    edges = []
    for i, row in enumerate(expected_rnaseq_array):
        if config.DEBUG_VERBOSE and i%10 == 0:
            config.log_statement( "Clustering bin %i/%s in RNAseq array" % (
                    i, expected_rnaseq_array.shape) )
        norm_row = row/row.sum()
        matching_indices = [ j for j, x in enumerate(norm_rows)
                             if float(numpy.abs(x-norm_row).sum()) < 1e-6  ]
        for j in matching_indices:
            edges.append(( i, j))
        norm_rows.append( norm_row  )
    
    graph = nx.Graph()
    graph.add_nodes_from(xrange(expected_rnaseq_array.shape[0]))
    graph.add_edges_from(edges)
    clusters = nx.connected_components(graph)
    print clusters
    assert False
    """
    
    new_expected_array = numpy.zeros( 
        (len(clusters), expected_rnaseq_array.shape[1]) )
    new_observed_array = numpy.zeros( len(clusters), dtype=int )
    cluster_mapping = {}
    for i, node in enumerate(clusters):
        cluster_mapping[i] = node
        new_expected_array[i,:] = expected_rnaseq_array[node,].sum(0)
        new_observed_array[i] = observed_rnaseq_array[node].sum()

    return new_expected_array, new_observed_array, cluster_mapping

def find_nonoverlapping_exons_covered_by_segment(exon_bndrys, start, stop):
    """Return the pseudo bins that a given segment has at least one basepair in.

    """
    bin_1 = exon_bndrys.searchsorted(start, side='right')-1
    # if the start falls before all bins
    if bin_1 == -1: return ()

    bin_2 = exon_bndrys.searchsorted(stop, side='right')-1
    # if the stop falls after all bins
    if bin_2 == len( exon_bndrys ) - 1: return ()
    
    if DEBUG:
        assert bin_1 == -1 or start >= exon_bndrys[ bin_1  ]
        assert stop < exon_bndrys[ bin_2 + 1  ]

    return tuple(xrange( bin_1, bin_2+1 ))
 

def bin_rnaseq_reads( reads, chrm, strand, exon_boundaries, include_read_type=True ):
    """Bin reads into non-overlapping exons.

    exon_boundaries should be a numpy array that contains
    pseudo exon starts.
    """
    if not reads.reads_are_stranded: strand = '.'
    
    # first get the paired reads
    gene_start = int(exon_boundaries[0])
    gene_stop = int(exon_boundaries[-1])
    paired_reads = list( reads.iter_paired_reads(
            chrm, strand, gene_start, gene_stop+1) )
    
    # find the unique subset of contiguous read sub-locations
    read_locs = set()
    for r in chain(*paired_reads):
        for start, stop in iter_coverage_intervals_for_read( r ):
            read_locs.add( (start, stop) )
    
    # build a mapping from contiguous regions into the non-overlapping exons (
    # ie, exon segments ) that they overlap
    read_locs_into_bins = {}
    for start, stop in read_locs:
        read_locs_into_bins[(start, stop)] = \
            find_nonoverlapping_exons_covered_by_segment( 
                exon_boundaries, start, stop )

    def build_bin_for_read( read ):
        bin = set()
        for start, stop in iter_coverage_intervals_for_read( read ):
            bin.update( read_locs_into_bins[(start, stop)] )
        return tuple(sorted(bin))
    
    # finally, aggregate the bins
    binned_reads = defaultdict( int )
    for r1, r2 in paired_reads:
        if r1.rlen == 0: 
            rlen = sum( x[1] for x in r1.cigar if x[0] == 0 )
        else: 
            rlen = r1.rlen
            if rlen != r2.rlen:
                if config.DEBUG_VERBOSE:
                    config.log_statement(
                        "WARNING: read lengths are not the same for %s and %s" % (
                            r1.qname, r2.qname),
                        log=True, display=False)
                    config.log_statement(
                        str(r1), log=True, display=False)
                    config.log_statement(
                        str(r2), log=True, display=False)                    
                continue
        
        rg = get_read_group( r1, r2 )
        bin1 = build_bin_for_read( r1 )
        bin2 = build_bin_for_read( r2 )
        # skip any reads that don't completely overlap the gene
        if bin1 == () or bin2== () or any(x==() for x in chain(bin1, bin2)): continue
        assert len(bin1) > 0
        assert len(bin2) > 0
        if include_read_type: key = ( rlen, rg, tuple(sorted((bin1,bin2))))
        else: key = tuple(sorted((bin1,bin2)))
        binned_reads[key] += 1
    
    return dict(binned_reads)

def bin_single_end_rnaseq_reads(reads, chrm, strand, exon_boundaries):
    # first get the paired reads
    gene_start = int(exon_boundaries[0])
    gene_stop = int(exon_boundaries[-1])
    paired_reads = list( reads.iter_paired_reads(
            chrm, strand, gene_start, gene_stop) )
    bins = defaultdict(int)
    for r in chain(*paired_reads):
        bin = []
        for start, stop in iter_coverage_intervals_for_read( r ):
            bin.extend(
                find_nonoverlapping_exons_covered_by_segment( 
                    exon_boundaries, start, stop ))
        bins[(r.inferred_length, tuple(bin))] += 1
    return bins

def build_element_expected_and_observed(
        rnaseq_reads, exon_boundaries, gene):
    binned_reads = bin_single_end_rnaseq_reads(
        rnaseq_reads, gene.chrm, gene.strand, exon_boundaries)
    num_reads = float(sum(binned_reads.values()))
    fragment_cnts = find_expected_single_bin_freqs(
        exon_boundaries, binned_reads.keys())
    return fragment_cnts, binned_reads
    #### calculate the read length normalized bin counts
    # estimate the marginal fraction of each read length
    read_len_freqs = defaultdict(int)
    for (read_len, bin), cnt in binned_reads.iteritems():
        read_len_freqs[read_len] += cnt/num_reads
    # normalize and group the bin type fragment counts
    rd_len_grpd_fragment_cnts = defaultdict(int)
    for (read_len, bin), cnt in fragment_cnts.iteritems():
        rd_len_grpd_fragment_cnts[bin] += cnt*read_len_freqs[read_len]
    num_frags = sum(rd_len_grpd_fragment_cnts.values())
    for bin, cnt in list(rd_len_grpd_fragment_cnts.iteritems()):
        rd_len_grpd_fragment_cnts[bin] = cnt
    
    # group binned reads by read length
    rd_len_grpd_binned_reads = defaultdict(int)
    for (read_len, bin), cnt in binned_reads.iteritems():
        rd_len_grpd_binned_reads[bin] += cnt
    return rd_len_grpd_fragment_cnts, rd_len_grpd_binned_reads

def find_expected_single_bin_freqs(exon_boundaries, read_lens_and_bins):
    def bin_len(i):
        return exon_boundaries[i+1] - exon_boundaries[i]
    
    bin_counts = {}
    for read_len, bin_indices in read_lens_and_bins:
        assert sum(bin_len(i) for i in bin_indices) >= read_len
        if len(bin_indices) == 1: 
            num_frags = bin_len(bin_indices[0]) - read_len
        elif len(bin_indices) == 2:
            num_frags = min(bin_len(bin_indices[0]), 
                            bin_len(bin_indices[1]), 
                            read_len)
        else:
            interior_len = sum(bin_len(i) for i in bin_indices[1:-1])
            num_frags = min(bin_len(bin_indices[0]), 
                            bin_len(bin_indices[-1]), 
                            read_len-interior_len)
        
        bin_counts[(read_len, bin_indices)] = num_frags
    return bin_counts

class DesignMatrix(object):
    def filter_design_matrix(self):        
        return
    
    def _build_rnaseq_arrays(self, gene, rnaseq_reads, fl_dists):
        # bin the rnaseq reads
        expected_rnaseq_cnts, observed_rnaseq_cnts = \
            build_expected_and_observed_rnaseq_counts( 
                gene, rnaseq_reads, fl_dists )
        clustered_bins = cluster_bins(expected_rnaseq_cnts)
        for cluster in clustered_bins:
            print cluster
        print
        
        # if no transcripts are observable given the fl dist, then return nothing
        if len( expected_rnaseq_cnts ) == 0:
            self.array_types.append('RNASeq')
            self.obs_cnt_arrays.append(None)
            self.expected_freq_arrays.append(None)
            return 
        
        # build the expected and observed counts, and convert them to frequencies
        ( expected_rnaseq_array, observed_rnaseq_array, unobservable_rnaseq_trans ) = \
              build_expected_and_observed_arrays( 
                expected_rnaseq_cnts, observed_rnaseq_cnts, normalize=True ) 

        del expected_rnaseq_cnts, observed_rnaseq_cnts
        
        if config.DEBUG_VERBOSE:
            config.log_statement( "Clustering bins in RNAseq array" )
        expected_rnaseq_array, observed_rnaseq_array, clusters = cluster_rows(
            expected_rnaseq_array, observed_rnaseq_array)
        
        self.array_types.append('RNASeq')
        self.obs_cnt_arrays.append(observed_rnaseq_array)
        self.expected_freq_arrays.append(expected_rnaseq_array)
        self.unobservable_transcripts.update(unobservable_rnaseq_trans)
    
    def _build_gene_bnd_arrays(self, gene, reads, reads_type):
        # bin the CAGE data
        expected_cnts, observed_cnts = \
            build_expected_and_observed_transcript_bndry_counts( 
            gene, reads )
        expected_array, observed_array, unobservable_trans = \
            build_expected_and_observed_arrays( 
            expected_cnts, observed_cnts, normalize=True )
        del expected_cnts, observed_cnts

        self.array_types.append(reads_type)
        self.obs_cnt_arrays.append(observed_array)
        self.expected_freq_arrays.append(expected_array)
        self.unobservable_transcripts.update(unobservable_trans)
        return
    
    def transcript_indices(self):
        """Sorted list transcript indices that the expected array was built for.
        """
        # it doesn't matter which design matric we use, because they 
        # al have the same number of transcripts
        for array in self.expected_freq_arrays:
            if array == None: continue
            num_transcripts = array.shape[1]
            break
        
        indices = set(xrange(num_transcripts)) - self.filtered_transcripts
        return numpy.array(sorted(indices))
    
    def expected_and_observed(self, bam_cnts=None):
        """Build expected and observed arrays. 
        
        If bam cnts is provided, then add the out-of-gene bins
        """
        # find the transcripts that we want to build the array for
        indices = self.transcript_indices()
        # add in the out of gene counts
        if bam_cnts != None:
            indices = [-1,] + indices.tolist()
            indices = numpy.array(indices)+1
        
        if self._expected_and_observed != None and \
                self._cached_bam_cnts == bam_cnts and \
                sorted(self._cached_indices) == sorted(indices):
            return self._expected_and_observed
        
        # stack all of the arrays, and filter out transcripts to skip
        exp_arrays_to_stack = []
        obs_arrays_to_stack = []
        for i, (expected, observed) in enumerate(izip(
                self.expected_freq_arrays, self.obs_cnt_arrays)):
            if expected == None:
                assert observed == None
                continue
            if bam_cnts != None: 
                #obs_arrays_to_stack.append( bam_cnts[i]-sum(observed) )
                observed = numpy.hstack((bam_cnts[i]-sum(observed), observed))
                expected = numpy.vstack(
                    (numpy.zeros(expected.shape[1]), expected))
                #print "zeros", numpy.zeros((expected.shape[0], 1)).shape, expected.shape,
                expected = numpy.hstack(
                    (numpy.zeros((expected.shape[0], 1)), expected))
                #print expected.shape,
                expected[0,0] = 1
            
            exp_arrays_to_stack.append(expected)
            obs_arrays_to_stack.append(observed)

        # stack all of the data type arrays
        expected = numpy.vstack(exp_arrays_to_stack)[:,indices]
        observed = numpy.hstack(obs_arrays_to_stack)
        
        # find which bins have 0 expected reads
        bins_to_keep = (expected.sum(1) > 1e-6)
        self._expected_and_observed = (
            expected[bins_to_keep,], observed[bins_to_keep])
        self._cached_bam_cnts = bam_cnts
        self._cached_indices = indices
        return self._expected_and_observed

    def find_transcripts_to_filter(self,expected,observed,max_num_transcripts):
        # cluster bins
        expected, observed = cluster_rows(expected, observed)
        
        num_transcripts = expected.shape[1]
        low_expression_ts = set(self.unobservable_transcripts)
        if num_transcripts <= max_num_transcripts: 
            return low_expression_ts
        

        # transcripts to remove
        test = (observed+1)/expected.max(1)
        for index in numpy.arange(observed.shape[0])[test.argsort()]:
            if num_transcripts - len(low_expression_ts) <= max_num_transcripts:
                break
            new_low_expression = set(expected[index,].nonzero()[0])
            # if this would remove every transcript, skip it
            if len( low_expression_ts.union(new_low_expression) ) == num_transcripts:
                continue
            low_expression_ts.update( new_low_expression )

        return low_expression_ts
    
    def __init__(self, gene, fl_dists,
                 rnaseq_reads, five_p_reads, three_p_reads,
                 max_num_transcripts=None):
        assert fl_dists != None
        self.array_types = []

        self.obs_cnt_arrays = []
        self.expected_freq_arrays = []
        self.unobservable_transcripts = set()

        self._cached_bam_cnts = None
        self._cached_indices = None
        
        self.filtered_transcripts = None
        self.max_num_transcripts = max_num_transcripts
        
        # initialize sum currently unset variables
        self.num_rnaseq_reads, self.num_fp_reads, self.num_tp_reads = (
            None, None, None)
        self._expected_and_observed = None
        
        if len( gene.transcripts ) == 0:
            raise ValueError, "No transcripts"
        
        if five_p_reads != None:
            if config.DEBUG_VERBOSE:
                config.log_statement( "Building TSS arrays for %s" % gene.id )
            self._build_gene_bnd_arrays(gene, five_p_reads, 'five_p_reads')
            self.num_fp_reads = sum(self.obs_cnt_arrays[-1])
        else:
            self.expected_freq_arrays.append(None)
            self.obs_cnt_arrays.append(None)
            self.num_fp_reads = None
        
        if config.DEBUG_VERBOSE:
            config.log_statement( "Building RNAseq arrays for %s" % gene.id )
        self._build_rnaseq_arrays(gene, rnaseq_reads, fl_dists)
        if self.obs_cnt_arrays[-1] != None:
            self.num_rnaseq_reads = sum(self.obs_cnt_arrays[-1])
        
        if three_p_reads != None:
            if config.DEBUG_VERBOSE:
                config.log_statement( "Building TES arrays for %s" % gene.id )
            self._build_gene_bnd_arrays(gene, three_p_reads, 'three_p_reads')
            self.num_tp_reads = sum(self.obs_cnt_arrays[-1])
        else:
            self.expected_freq_arrays.append(None)
            self.obs_cnt_arrays.append(None)
            self.num_tp_reads = None

        if all( mat == None for mat in self.obs_cnt_arrays ):
            raise NoObservableTranscriptsError, "No observable transcripts"
        
        # initialize the filtered_transcripts to the unobservable transcripts
        self.filtered_transcripts = set(list(self.unobservable_transcripts))
        
        # update the set to satisfy the max_num_transcripts restriction
        if max_num_transcripts != None:
            if config.DEBUG_VERBOSE:
                config.log_statement("Filtering design matrix for %s" % gene.id)

            expected, observed = self.expected_and_observed()
            
            self.filtered_transcripts =  self.find_transcripts_to_filter(
                expected, observed, max_num_transcripts)
        
        return


def tests( ):
    exon_lens = [100,1,100]
    transcript = range( len(exon_lens) )
    fl_dist = frag_len.build_uniform_density( 100, 100 )
    read_len = 50
    
    full, paired, single = find_possible_read_bins_for_transcript( 
        transcript, exon_lens, fl_dist, read_len )

    for bin in sorted(full):
        print "Fr Bin:", bin
    
    bins = sorted( full )
    bin_counts = simulate_reads_from_exons(10000, fl_dist, None, 0, exon_lens)
    total_cnt = float(sum(bin_counts.values()))
    simulated_freqs = dict((bin, cnt/total_cnt) 
                           for bin, cnt in bin_counts.iteritems())
    analytical_cnts = dict( 
        ( bin, calc_long_read_expected_cnt_given_bin( 
                bin, exon_lens, fl_dist ) ) 
        for bin in bins )
    
    total_cnt = sum( analytical_cnts.values() )
    analytical_freqs = dict( ( bin, float(cnt)/total_cnt ) \
                             for bin, cnt in analytical_cnts.iteritems() )
    print str("Frag Bin").ljust( 30 ), "Est %/cnt \tObs %/cnt"
    for bin in bins:
        sim_freq = simulated_freqs[bin] if bin in simulated_freqs else 0.000
        sim_cnt = bin_counts[bin] if bin in bin_counts else 0
        print "%s %.3f/%i  \t%.3f/%i" % (
            str(bin).ljust( 30 ), 
            analytical_freqs[bin], analytical_cnts[bin],
            sim_freq, sim_cnt )

    print "SHORT READ SIMS:"    
    
    #fl_dist = frag_len.build_uniform_density( 105, 400 )
    #exon_lens = [ 500, 500, 50, 5, 5, 500, 500 ]
    #read_len = 100
    max_num_unmappable_bases = 1
    transcript = range( len(exon_lens) )
    
    print "Transcript:", transcript
    print "Exon Lens:", exon_lens
    print fl_dist

    exon_boundaries = numpy.array(exon_lens).cumsum().tolist()
    exon_boundaries.insert(0, 0)
    bin_counts = simulate_reads_from_exons( 
        1000000, fl_dist, read_len, 0, exon_lens  )
    total_cnt = float(sum(bin_counts.values()))
    simulated_freqs = dict( ( bin, cnt/total_cnt ) 
                            for bin, cnt in bin_counts.iteritems() )
    analytical_cnts = calc_expected_cnts(
        exon_boundaries, [tuple(transcript),], fl_dist,
        read_len, read_len, max_num_unmappable_bases )[tuple(transcript)]
    total_cnt = sum( analytical_cnts.values() )
    analytical_freqs = dict( ( bin, float(cnt)/total_cnt ) 
                             for bin, cnt in analytical_cnts.iteritems() )
    print str("Frag Bin").ljust( 40 ), "Est %/cnt  \tObs %/cnt"
    for bin in sorted(analytical_freqs.keys()):
        sim_freq = simulated_freqs[bin] if bin in simulated_freqs else 0.000
        sim_cnt = bin_counts[bin] if bin in bin_counts else 0
        print "%s %.3f/%i\t%.3f/%i" % (
            str(bin).ljust( 40 ), 
            analytical_freqs[bin], analytical_cnts[bin],
            sim_freq, sim_cnt )
    
    assert False

    global foo
    def foo():
        for loop in xrange(100):
            analytical_cnts = dict( ( bin, estimate_num_paired_reads_from_bin( \
                    bin, transcript, exon_lens, fl_dist, \
                    read_len, max_num_unmappable_bases ) ) \
                            for bin in bins )

    #import cProfile    
    #cProfile.run("foo()")
    
    return
        
if __name__ =='__main__':
    import grit
    import grit.config as config
    tests()
