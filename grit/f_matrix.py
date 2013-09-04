# Copyright (c) 2011-2012 Nathan Boley

import sys

import numpy
MAX_NUM_UNMAPPABLE_BASES = 0
RUN_INLINE_TESTS = False
LET_READS_OVERLAP = True
PROMOTER_SIZE = 50

DEBUG=False

from igraph import Graph

import frag_len

from itertools import product, izip, chain
from collections import defaultdict

from files.reads import ( iter_coverage_regions_for_read, get_read_group,
                          CAGEReads, RAMPAGEReads, PolyAReads )


################################################################################
#
#
# Code for converting transcripts into 
#
#

def find_nonoverlapping_boundaries( transcripts ):
    boundaries = set()
    for transcript in transcripts:
        for exon in transcript.exons:
            boundaries.add( exon[0] )
            boundaries.add( exon[1]+1 )
    
    return numpy.array( sorted( boundaries ) )

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
            if fl >= exon_lens_cum[-1]:
                continue
            if read_len != None \
                and not LET_READS_OVERLAP \
                and fl < 2*read_len:
                continue
            
            def get_bin( position ):
                # find the start exon for the first read
                exon = exon_lens_cum.searchsorted( position )
                
                if position > 0:
                    exon -= 1
                    
                return exon
            
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
            start_location = random.randrange( 0, exon_lens_cum[-1] - fl )
            end_location = start_location + fl
            
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
                pair_2 = get_bin_pairs( end_location - read_len, read_len )
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
        # from this exon.
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
        
        min_fl = max( exon_lens[0] - x + gap + 0, fl_dist.fl_min )
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

def find_possible_read_bins( trans_exon_indices, exon_lens, fl_dist, \
                             read_len, min_num_mappable_bases=1 ):
    """Find all possible bins. 

    """
    full_fragment_bins = set()
    paired_read_bins = set()
    single_end_read_bins = set()
    
    transcript_exon_lens = numpy.array([ exon_lens[i] for i in trans_exon_indices ])
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
        
            paired_read_bins.update( paired )
            single_end_read_bins.update( single )
    
    return full_fragment_bins, paired_read_bins, single_end_read_bins

def estimate_num_paired_reads_from_bin( 
        bin, transcript, exon_lens,
        fl_dist, read_len, min_num_mappable_bases=1 ):
    """

    """
    # calculate the exon lens for the first and second reads
    fr_exon_lens = [ exon_lens[i] for i in bin[0]  ]
    sr_exon_lens = [ exon_lens[i] for i in bin[1]  ]
    # find the length of the internal exons ( those not in either bin )
    pre_sr_exon_lens = sum( exon_lens[i] for i in transcript \
                            if i < bin[1][0] and i >= bin[0][0] )

    if DEBUG:
        print fr_exon_lens
        print sr_exon_lens
        print pre_sr_exon_lens
    
    # loop through the possible start indices for the first bin to be satisifed
    # basically, after removing the spliced introns, we need at least min_num_ma
    # in the first exon, and in the last exon. This puts a lower bound on the 
    # read start at position ( relative to the transcript ) of 0 ( the read can 
    # start at the first position in exon 1 ) *or* last_exon_start + 1-read_len 
    min_start = max( 0, sum( fr_exon_lens[:-1] ) \
                         + min_num_mappable_bases - read_len )
    # make sure that there is enough room for the second read to start in 
    # sr_exon[0] given the fragment length constraints
    min_start = max( min_start, pre_sr_exon_lens + read_len - fl_dist.fl_max )

    max_start = fr_exon_lens[0] - min_num_mappable_bases
    # make sure that there are enough bases in the last exon for the read to fit
    max_start = min( max_start, sum( fr_exon_lens ) - read_len )
    
    if DEBUG:
        print min_start, max_start
    
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
        print min_stop, max_stop
    
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
    
    return float( density )

def calc_expected_cnts( exon_boundaries, transcripts, fl_dists_and_read_lens, \
                        max_num_unmappable_bases=MAX_NUM_UNMAPPABLE_BASES,
                        max_memory_usage=3.5 ):
    # store all counts, and count vectors. Indexed by ( 
    # read_group, read_len, bin )
    cached_f_mat_entries = {}
    f_mat_entries = {}
    
    nonoverlapping_exon_lens = \
        numpy.array([ stop - start for start, stop in 
                      izip(exon_boundaries[:-1], exon_boundaries[1:])])
    
    # for each candidate trasncript
    n_bins = 0
    for fl_dist, read_len in fl_dists_and_read_lens:
        for transcript_index, nonoverlapping_indices in enumerate(transcripts):
            if (n_bins*len(transcripts)*8.)/(1024**3) > max_memory_usage:
                raise MemoryError, "Building the design matrix has exceeded the maximum allowed memory "
            
            # find all of the possible read bins for transcript given this 
            # fl_dist and read length
            full, pair, single =  \
                find_possible_read_bins( nonoverlapping_indices, \
                                         nonoverlapping_exon_lens, fl_dist, \
                                         read_len, min_num_mappable_bases=1 )
            
            # add the expected counts for paired reads
            for bin in pair:
                key = (read_len, hash(fl_dist), bin)
                if key in cached_f_mat_entries:
                    pseudo_cnt = cached_f_mat_entries[ key ]
                else:
                    pseudo_cnt = estimate_num_paired_reads_from_bin(
                        bin, nonoverlapping_indices, 
                        nonoverlapping_exon_lens, fl_dist,
                        read_len, max_num_unmappable_bases )

                    cached_f_mat_entries[ key ] = pseudo_cnt
                
                if not f_mat_entries.has_key( key ):
                    n_bins += 1
                    f_mat_entries[key]= numpy.zeros( 
                        len(transcripts), dtype=float )
                    
                f_mat_entries[key][ transcript_index ] = pseudo_cnt
        
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
        rv[ ( read_len, hash(fl_dists[read_group]), bin ) ] = value
    
    return rv

def build_expected_and_observed_arrays( 
        expected_cnts, observed_cnts, normalize=True ):
    expected_mat = []
    observed_mat = []
    unobservable_transcripts = set()
    
    for key, val in sorted(expected_cnts.iteritems()):
        # skip bins with 0 expected reads
        if sum( val) == 0:
            continue
        
        expected_mat.append( val )
        try:
            observed_mat.append( observed_cnts[key] )
        except KeyError:
            observed_mat.append( 0 )
    
    if len( expected_mat ) == 0:
        raise ValueError, "No expected reads."
    
    expected_mat = numpy.array( expected_mat, dtype=numpy.double )
    
    if normalize:
        nonzero_entries = expected_mat.sum(0).nonzero()[0]
        unobservable_transcripts = set(range(expected_mat.shape[1])) \
            - set(nonzero_entries.tolist())
        observed_mat = numpy.array( observed_mat, dtype=numpy.int )
        expected_mat = expected_mat[:,nonzero_entries]
        expected_mat = expected_mat/expected_mat.sum(0)
    
    return expected_mat, observed_mat, unobservable_transcripts

def build_expected_and_observed_rnaseq_counts( gene, reads, fl_dists ):    
    # find the set of non-overlapping exons, and convert the transcripts to 
    # lists of these non-overlapping indices. All of the f_matrix code uses
    # this representation.     
    exon_boundaries = find_nonoverlapping_boundaries(gene.transcripts)
    transcripts_non_overlapping_exon_indices = \
        list(build_nonoverlapping_indices( 
                gene.transcripts, exon_boundaries ))
    
    binned_reads = bin_reads( 
        reads, gene.chrm, gene.strand, exon_boundaries)
        
    observed_cnts = build_observed_cnts( binned_reads, fl_dists )    
    read_groups_and_read_lens =  { (RG, read_len) for RG, read_len, bin 
                                   in binned_reads.iterkeys() }
    
    fl_dists_and_read_lens = [ (fl_dists[RG], read_len) for read_len, RG  
                               in read_groups_and_read_lens ]
    
    expected_cnts = calc_expected_cnts( 
        exon_boundaries, transcripts_non_overlapping_exon_indices, 
        fl_dists_and_read_lens)
    
    return expected_cnts, observed_cnts

def build_expected_and_observed_transcript_bndry_counts( 
        gene, reads, bndry_type=None ):
    if bndry_type == None:
        # try to infer the boundary type from the reads type
        if isinstance(reads, CAGEReads): bndry_type = "CAGE"
        elif isinstance(reads, RAMPAGEReads): bndry_type = "CAGE"
        elif isinstance(reads, PolyAReads): bndry_type = "POLYA"
        else: assert False, "Can't infer boundary read type (%s)" % type(reads)
    assert bndry_type in ["CAGE", "POLYA"]
    
    cvg_array = numpy.zeros(gene.stop-gene.start+1)
    cvg_array += reads.build_read_coverage_array( 
        gene.chrm, gene.strand, gene.start, gene.stop )
    
    # find all the bndry peak regions
    peaks = list()
    for transcript in gene.transcripts:
        if bndry_type == 'CAGE':
            peaks.append( transcript.find_promoter() )
        elif bndry_type == 'POLYA':
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

def cluster_rows(expected_rnaseq_array, observed_rnaseq_array):
    corr_coefs = numpy.corrcoef(expected_rnaseq_array)
    if isinstance(corr_coefs, numpy.ndarray):
        edges = set()
        row, col = (corr_coefs == 1).nonzero()
        for row_i, col_i in zip( row, col ):
            edges.add((int(row_i), int(col_i)))

        graph = Graph( expected_rnaseq_array.shape[0] )
        graph.add_edges(edges) #[ (int(e[0]), int(e[1])) for e in edges])
        clusters = graph.clusters()
        
        new_expected_array = numpy.zeros( 
            (len(clusters), expected_rnaseq_array.shape[1]) )
        new_observed_array = numpy.zeros( len(clusters), dtype=int )
        for i, node in enumerate(graph.clusters()):
            new_expected_array[i,:] = expected_rnaseq_array[node,].sum(0)
            new_observed_array[i] = observed_rnaseq_array[node].sum()

        expected_rnaseq_array = new_expected_array
        observed_rnaseq_array = new_observed_array
    
    return expected_rnaseq_array, observed_rnaseq_array

def build_design_matrices( gene, rnaseq_reads, fl_dists, all_promoter_reads=[], 
                           max_num_transcripts=None  ):
    if len( gene.transcripts ) == 0:
        return numpy.zeros(0), numpy.zeros(0), []
        
    # bin the rnaseq reads
    expected_rnaseq_cnts, observed_rnaseq_cnts = \
        build_expected_and_observed_rnaseq_counts( 
            gene, rnaseq_reads, fl_dists )
    if len( expected_rnaseq_cnts ) == 0:
        return numpy.zeros(0), numpy.zeros(0), []
        
    expected_rnaseq_array, observed_rnaseq_array, unobservable_rnaseq_trans = \
        build_expected_and_observed_arrays( 
            expected_rnaseq_cnts, observed_rnaseq_cnts, True )
    del expected_rnaseq_cnts, observed_rnaseq_cnts
    
    expected_rnaseq_array, observed_rnaseq_array = cluster_rows(
        expected_rnaseq_array, observed_rnaseq_array)
    
    num_transcripts = expected_rnaseq_array.shape[1]
    if max_num_transcripts != None and num_transcripts > max_num_transcripts:
        low_expression_ts = set()
        test = (observed_rnaseq_array+1)/expected_rnaseq_array.max(1)
        for index in numpy.arange(observed_rnaseq_array.shape[0])[test.argsort()]:
            if num_transcripts - len(low_expression_ts) <= max_num_transcripts:
                break
            new_low_expression = set(expected_rnaseq_array[index,].nonzero()[0])
            if len( low_expression_ts.union(new_low_expression) ) == num_transcripts:
                break
            low_expression_ts.update( new_low_expression )
        
        high_expression_ts = numpy.array(
            list(set(range(num_transcripts)) - low_expression_ts))
        expected_rnaseq_array = expected_rnaseq_array[:,high_expression_ts]
        
        new_unobservable_transcripts = [] + list(unobservable_rnaseq_trans)
        for low_exp_t_i in low_expression_ts:
            new_unobservable_transcripts.append( low_exp_t_i + sum(
                    x <= low_exp_t_i for x in unobservable_rnaseq_trans ))
        unobservable_rnaseq_trans = set( new_unobservable_transcripts )

        expected_rnaseq_array, observed_rnaseq_array = cluster_rows(
            expected_rnaseq_array, observed_rnaseq_array)
    
    # rename the arrays, in case there is no 
    expected_array = expected_rnaseq_array
    observed_array = observed_rnaseq_array
    unobservable_trans = unobservable_rnaseq_trans
    
    # deal with the, optional, promoter reads
    for promoter_reads in all_promoter_reads:
        # bin the CAGE data
        expected_promoter_cnts, observed_promoter_cnts = \
            build_expected_and_observed_transcript_bndry_counts( 
            gene, promoter_reads )
        expected_prom_array, observed_prom_array, unobservable_prom_trans = \
            build_expected_and_observed_arrays( 
            expected_promoter_cnts, observed_promoter_cnts, False )
        del expected_promoter_cnts, observed_promoter_cnts

        # combine the arrays
        observed_array = numpy.delete( observed_array, 
                      numpy.array(list(unobservable_prom_trans)) )
        observed_prom_array = numpy.delete( observed_prom_array, 
                      numpy.array(list(unobservable_trans)) )
        observed_array = numpy.hstack((observed_prom_array, observed_array))

        expected_array = numpy.delete( expected_array, 
                      numpy.array(list(unobservable_prom_trans)), axis=1 )
        expected_prom_array = numpy.delete( expected_prom_array, 
                      numpy.array(list(unobservable_trans)), axis=1 )   

        expected_array = numpy.vstack((expected_prom_array, expected_array))
        unobservable_trans = unobservable_trans.union(
            unobservable_prom_trans)
    
    return expected_array, observed_array, unobservable_trans

def find_nonoverlapping_exons_covered_by_segment(exon_bndrys, start, stop):
    """Return the pseudo bins that a given segment has at least one basepair in.

    """
    # BUG XXX
    start += 1
    stop +=1 
    
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
 

def bin_reads( reads, chrm, strand, exon_boundaries ):
    """Bin reads into non-overlapping exons.

    exon_boundaries should be a numpy array that contains
    pseudo exon starts.
    """
    # first get the paired reads
    gene_start = int(exon_boundaries[0])
    gene_stop = int(exon_boundaries[-1])
    paired_reads = list( reads.iter_paired_reads(
            chrm, strand, gene_start, gene_stop) )
    
    # find the unique subset of contiguous read sub-locations
    read_locs = set()
    for r in chain(*paired_reads):
        for chrm, strand, start, stop in iter_coverage_regions_for_read( 
                r, reads, reads.reverse_read_strand,reads.pairs_are_opp_strand):
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
        for chrm, strand, start, stop in iter_coverage_regions_for_read(
                read, reads, 
                reads.reverse_read_strand, reads.pairs_are_opp_strand):
            bin.update( read_locs_into_bins[(start, stop)] )
        return tuple(sorted(bin))
    
    # finally, aggregate the bins
    binned_reads = defaultdict( int )
    for r1, r2 in paired_reads:
        if r1.rlen != r2.rlen:
            print >> sys.stderr, "WARNING: read lengths are then same"
            continue
        
        rlen = r1.rlen
        rg = get_read_group( r1, r2 )
        bin1 = build_bin_for_read( r1 )
        bin2 = build_bin_for_read( r2 )
        binned_reads[( rlen, rg, tuple(sorted((bin1,bin2))))] += 1
    
    return dict(binned_reads)


def tests( ):
    exon_lens = [ 500, 500, 5, 5, 5, 100, 200, 1000]
    transcript = range( len(exon_lens) )
    fl_dist = frag_len.build_uniform_density( 110, 600 )
    read_len = 100
    
    full, paired, single = find_possible_read_bins( \
        transcript, exon_lens, fl_dist, read_len )

    for bin in sorted(full):
        print "Fr Bin:", bin
    
    bins = sorted( full )
    bin_counts = simulate_reads_from_exons(10000, fl_dist, None, 0, exon_lens)
    total_cnt = float(sum(bin_counts.values()))
    simulated_freqs = dict((bin, cnt/total_cnt) 
                           for bin, cnt in bin_counts.iteritems())
    
    analytical_cnts = dict( ( bin, calc_long_read_expected_cnt_given_bin( \
                                       bin, exon_lens, fl_dist ) ) \
                                for bin in bins )
    
    total_cnt = sum( analytical_cnts.values() )
    analytical_freqs = dict( ( bin, float(cnt)/total_cnt ) \
                             for bin, cnt in analytical_cnts.iteritems() )

    for bin in bins:
        sim = simulated_freqs[bin] if bin in simulated_freqs else 0.000
        print str(bin).ljust( 30 ), round(analytical_freqs[bin],3 ),round(sim,3)

    print "SHORT READ SIMS:"    

    fl_dist = frag_len.build_uniform_density( 105, 400 )
    exon_lens = [ 500, 500, 50, 5, 5, 500, 500 ]
    read_len = 100
    max_num_unmappable_bases = 1
    transcript = range( len(exon_lens) )
    
    print transcript
    print exon_lens
    print fl_dist.fl_min, fl_dist.fl_max
    
    full, paired, single = find_possible_read_bins( \
        transcript, exon_lens, fl_dist, read_len )    
    bins = sorted( paired )
    
    bin_counts = simulate_reads_from_exons( 
        10000, fl_dist, read_len, 0, exon_lens  )
    total_cnt = float(sum(bin_counts.values()))
    simulated_freqs = dict( ( bin, cnt/total_cnt ) \
                            for bin, cnt in bin_counts.iteritems() )
    
    analytical_cnts = dict( ( bin, estimate_num_paired_reads_from_bin( \
                    bin, transcript, exon_lens, fl_dist, \
                    read_len, max_num_unmappable_bases ) ) \
                            for bin in bins )
            
    total_cnt = sum( analytical_cnts.values() )
    analytical_freqs = dict( ( bin, float(cnt)/total_cnt ) \
                             for bin, cnt in analytical_cnts.iteritems() )
    
    for bin in bins:
        sim = simulated_freqs[bin] if bin in simulated_freqs else 0.000
        print str(bin).ljust( 40 ), \
              str(round(analytical_freqs[bin],4 )).ljust(8), \
              str(round(sim,4 )).ljust(8)
    

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
    tests()
    
