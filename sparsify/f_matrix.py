# Copyright (c) 2011-2012 Nathan Boley

import numpy
MAX_NUM_UNMAPPABLE_BASES = 0
RUN_INLINE_TESTS = False
LET_READS_OVERLAP = True

DEBUG=False

import frag_len

from itertools import product

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
                while read_start + read_len > exon_lens_cum[ start_exon + loop ]:
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

def find_possible_read_bins( transcript, exon_lens, fl_dist, \
                             read_len, min_num_mappable_bases=1 ):
    """Find all possible bins. 

    """
    full_fragment_bins = set()
    paired_read_bins = set()
    single_end_read_bins = set()
    
    exons = numpy.array( list( transcript ) )
    transcript_exon_lens = numpy.array([ exon_lens[exon] for exon in transcript  ])
    transcript_exon_lens_cumsum = transcript_exon_lens.cumsum()
    transcript_exon_lens_cumsum = numpy.insert( transcript_exon_lens_cumsum, 0, 0 )
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
            
            bin = tuple(exons[index_1:index_2+1])
        
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
    # basically, after removing the spliced introns, we need at least min_num_ma..
    # in the first exon, and in the last exon. This puts a lower bound on the read start
    # at position ( relative to the transcript ) of 0 ( the read can start at the first 
    # position in exon 1 ) *or* last_exon_start + 1 - read_len 
    min_start = max( 0, sum( fr_exon_lens[:-1] ) + min_num_mappable_bases - read_len )
    # make sure that there is enough room for the second read to start in sr_exon[0]
    # given the fragment length constraints
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
    # second, we know that at least min_num_mappable_bases are in the last exon of 
    # second read, because we assume that this read is possible. 
    min_stop = max( min_stop, pre_sr_exon_lens \
                              + sum( sr_exon_lens[:-1] ) \
                              + min_num_mappable_bases )

    # max stop can never be greater than the transcript length
    max_stop = pre_sr_exon_lens + sum( sr_exon_lens )
    # all but min_num_mappable_bases basepair is in the first exon of second read
    # again, we assume that at least min_num_mappable_bases bases are in the
    # last exon because the read is possible. However, it could be that this extends 
    # *past* the last exon, which we account for on the line after.x
    max_stop = min( max_stop, \
                    pre_sr_exon_lens + sr_exon_lens[0] \
                    - min_num_mappable_bases + read_len )
    # all basepairs of the last exon in the second read are in the fragment. Make sure
    # the previous calculation doesnt extend past the last exon.
    max_stop = min( max_stop, min_stop - min_num_mappable_bases + sr_exon_lens[-1] )
    
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
                density -= fl_dist.fl_density_cumsum[ min_fl - fl_dist.fl_min - 1 ]
        return density
    
    density = do()
    
    return float( density )

cached_f_mat_entries = {}
def build_f_matrix( transcripts, binned_reads, gene, fl_dists, \
                    max_num_unmappable_bases=MAX_NUM_UNMAPPABLE_BASES ):    
    num_trans = len( transcripts )

    # store all counts, and count vectors. Indexed by ( 
    # read_group, read_len, bin )
    f_mat_entries = {}

    # for each candidate trasncript
    for transcript_index, transcript in enumerate( transcripts ):
        # build the transcript composed of pseudo ( non overlapping ) exons. 
        # This means splitting the overlapping parts into new 'exons'
        nonoverlapping_transcript = \
            transcript.build_nonoverlapping_transcript( gene )
        nonoverlapping_exon_lens = numpy.array( gene.nonoverlapping_exon_lens )

        # loop through each read type - seperated by fl dist and read length
        for read_len, fl_group in binned_reads.read_groups:
            fl_dist = fl_dists[ fl_group ]
            
            # get the number of reads for this read type. Note that this isn't 
            # precisely correct - we should be scaling by the experiment cnts, 
            # but it should be pretty close
            num_reads = binned_reads.read_group_cnts[ ( read_len, fl_group ) ]
            
            # find all of the possible read bins for transcript given this 
            # fl_dist and read length
            full, pair, single =  \
                find_possible_read_bins( nonoverlapping_transcript, \
                                         nonoverlapping_exon_lens, fl_dist, \
                                         read_len, min_num_mappable_bases=1 )
            
            # add the expected counts for paired reads
            pseudo_cnts = []
            for bin in pair:
                key = (read_len, fl_group, bin)
                if key in cached_f_mat_entries:
                    pseudo_cnt = cached_f_mat_entries[ key ]
                else:
                    pseudo_cnt = estimate_num_paired_reads_from_bin(
                        bin, nonoverlapping_transcript, 
                        nonoverlapping_exon_lens, fl_dist,
                        read_len, max_num_unmappable_bases )
                    
                    cached_f_mat_entries[ key ] = pseudo_cnt

                # scale the count by the number of reads. If there are twice as 
                # many reads in experiment 1 as 2, then the expected ratio of 
                # reads is 2:1 even if the expression levels are identical
                pseudo_cnts.append( pseudo_cnt*num_reads )
            
            total_pseudo_cnts = float(sum(pseudo_cnts))
            
            for bin, pseudo_cnt in zip(pair, pseudo_cnts):
                if not f_mat_entries.has_key( (read_len, fl_group, bin) ):
                    f_mat_entries[(read_len, fl_group, bin)]= [0]*num_trans
                
                f_mat_entries[(read_len, fl_group, bin)][ transcript_index ] \
                    =  pseudo_cnt
            

    # associate the read counts with each key    
    new_f_mat = {}
    for key in f_mat_entries:
        read_len, read_group, bin = key
        try:
            read_cnts = binned_reads.binned_reads[ (read_len, read_group, bin) ]
        except KeyError:
            read_cnts = 0
        
        # dont allow unexplainable reads
        if sum(f_mat_entries[key]) == 0:
            continue
        
        new_f_mat[key] = [ f_mat_entries[key], read_cnts ]
    
    return new_f_mat

def convert_f_matrices_into_arrays( f_mats ):
    expected_cnts = []
    observed_cnts = []
        
    for key, (expected, observed) in f_mats.iteritems():
        expected_cnts.append( expected )
        observed_cnts.append( observed )
    
    # normalize the expected counts to be fraction of the highest
    # read depth isoform
    expected_cnts = numpy.array( expected_cnts )
    expected_cnts = expected_cnts/expected_cnts.sum(0).max()
    
    return expected_cnts, numpy.array( observed_cnts )


def tests( ):
    from transcripts import Transcript
    
    exon_lens = [ 500, 500, 5, 5, 5, 100, 200, 1000]
    transcript = Transcript((0,),(len(exon_lens)-1,), range(1,len(exon_lens)-1))
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
    transcript = Transcript((0,), (len(exon_lens)-1,),range(1,len(exon_lens)-1))
    
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
    
