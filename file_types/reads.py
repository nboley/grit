import sys, os
from itertools import chain
from collections import defaultdict

DEBUG = False

def clean_chr_name( chrm ):
    if chrm.startswith( "chr" ):
        chrm = chrm[3:]
    # convert the dmel chrm to M, to be consistent with ucsc
    if chrm == 'dmel_mitochondrion_genome':
        chrm = "M"
    return chrm

def get_chrm( read, bam_obj ):
    chrm = bam_obj.getrname( read.tid )
    chrm = clean_chr_name( chrm )
    return chrm

def get_strand( read, reverse_read_strand, pairs_are_opp_strand ):
    # make sure that the read is on the correct strand
    if pairs_are_opp_strand  \
            and ( ( read.is_read1 and not read.is_reverse ) \
                      or ( read.is_read2 and read.is_reverse ) ):
        strand = '+'
    elif not pairs_are_opp_strand and not read.is_reverse:
        strand = '+'
    else:
        strand = '-'

    if reverse_read_strand:
        if strand == '+':
            strand = '-'
        else:
            strand = '+'

    return strand

def get_read_group( r1, r2 ):        
    return 'mean'
    r1_read_group = [ val for key, val in r1.tags if key == 'RG' ]
    r1_read_group = r1_read_group[0] if len( r1_read_group ) == 1 else 'mean'
    r2_read_group = [ val for key, val in r2.tags if key == 'RG' ]
    r2_read_group = r2_read_group[0] if len( r2_read_group ) == 1 else 'mean'
    if r1_read_group == r2_read_group:
        return r1_read_group
    else: 
        print "WARNING: Read groups do not match."
        return None


def read_pairs_are_on_same_strand( bam_obj, num_reads_to_check=25000 ):
    # keep track of which fractiona re on the sam strand
    paired_cnts = {'no_mate': 0, 'same_strand': 1e-4, 'diff_strand': 1e-4}
    
    num_good_reads = 0
    num_observed_reads = 0
    for read in bam_obj:
        num_observed_reads += 1
        if read.is_paired and read.mate_is_unmapped:
            continue

        if not read.is_paired:
            paired_cnts['no_mate'] += 1        
        elif read.is_reverse != read.mate_is_reverse:
            paired_cnts['diff_strand'] += 1
        else:
            paired_cnts['same_strand'] += 1
            
        # keep collecting reads until we observe enough
        num_good_reads += 1
        if num_good_reads > num_reads_to_check:
            break
    
    # if the reads are single ended, then return True ( 
    #    because it doesnt really matter )
    if paired_cnts['no_mate'] == num_reads_to_check:
        return True
    
    if float(paired_cnts['same_strand'])/paired_cnts['diff_strand'] > 10:
        return True
    elif float(paired_cnts['diff_strand'])/paired_cnts['same_strand'] > 10:
        return False
    else:
        print >> sys.stderr, "Paired Cnts:", paired_cnts, "Num Reads", num_observed_reads
        raise ValueError, "Reads appear to be a mix of reads on the same and different strands."

def iter_coverage_regions_for_read( 
    read, bam_obj, reverse_read_strand, pairs_are_opp_strand ):
    """Find the regions covered by this read

    """
    strand = get_strand( read, reverse_read_strand, pairs_are_opp_strand )

    # get the chromosome, correcting for alternate chrm names
    chrm = get_chrm( read, bam_obj )

    # we loop through each contig in the cigar string to deal
    # with junctions reads.
    # add 1 to the start because bam files are 0 based
    rv = []
    start = read.pos
    for contig_type, length in read.cigar:
        # if this is a match, add it 
        if contig_type == 0:
            yield ( chrm, strand, start, start + length - 1 )
            start += length
        # skip reference insertions
        elif contig_type == 1:
            pass
        # move past refernce deletions
        elif contig_type == 2:
            start += length
        # skip past skipped regions
        elif contig_type == 3:
            start += length
        # skip past soft clipped regions, because the
        # actual aligned sequence doesnt start until we've moved
        # past the clipped region
        elif contig_type == 4:
            start += length
        # hard clipped regions are not present int he aligned 
        # sequence, so do nothing
        elif contig_type == 5:
            pass
        else:
            print >> sys.stderr, "Unrecognized cigar format:", read.cigar

    return

def clean_qname( qname ):
    if qname[-2] == '/' and qname[-1] in '12':
        return qname[:-2]
    return qname

try:
    import pysam

    class Reads( pysam.Samfile ):
        """Subclass the samfile object to include a method that returns reads 
           and their pairs.


        """
        def iter_paired_reads(self, chrm, strand, start, stop, 
                              min_read_len=0, ignore_partial_alignments=True):
            chrm = clean_chr_name( chrm )
            if chrm not in self.references:
                chrm = 'chr' + chrm
            
            # whether or not the gene is on the positive strand
            gene_strnd_is_rev = ( strand == '-' )
            #chrm = clean_chr_name( chrm )

            # get all of the first pairs
            def iter_pair1_reads():
                for read in self.fetch( chrm, start, stop  ):
                    if read.alen < min_read_len: continue
                    if ignore_partial_alignments and read.alen != read.qlen:
                        continue
                    if get_strand( read, False, False ) == strand:
                        yield read

            # index the pair 2 reads
            reads_pair2 = {}
            for read in self.fetch( chrm, start, stop ):
                if read.alen < min_read_len: continue
                if ignore_partial_alignments and read.alen != read.qlen:
                    continue
                if get_strand( read, False, False ) != strand:
                    if clean_qname(read.qname) in reads_pair2:
                        if DEBUG:
                            print "Multiple reads on the same strand"
                    reads_pair2[ clean_qname(read.qname) ] = read

            # iterate through the read pairs
            for read1 in iter_pair1_reads():
                try:
                    read2 = reads_pair2[ clean_qname(read1.qname) ]
                # if there is no mate, skip this read
                except KeyError:
                    if DEBUG:
                        print "No mate: ", read1.pos, read1.aend, clean_qname(read1.qname)
                    continue
                
                assert ( read1.qlen == read1.aend - read1.pos ) \
                       or ( len( read1.cigar ) > 1 )
                assert ( read2.qlen == read2.aend - read2.pos ) \
                       or ( len( read2.cigar ) > 1 )

                if read1.qlen != read2.qlen:
                    if DEBUG:
                        print( "ERROR: unequal read lengths %i and %i\n", \
                                   read1.qlen, read2.qlen )
                    continue

                yield read1, read2

            return

except ImportError:
    pass


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
 

def bin_reads( reads, chrm, strand, exon_boundaries, 
               reverse_read_strand, pairs_are_opp_strand ):
    """Bin reads into non-overlapping exons.

    exon_boundaries should be a numpy array that contains
    pseudo exon starts.
    """
    assert isinstance( reads, Reads )
    
    # first get the paired reads
    gene_start = int(exon_boundaries[0])
    gene_stop = int(exon_boundaries[-1])
    paired_reads = list( reads.iter_paired_reads(
            chrm, strand, gene_start, gene_stop) )
    
    # find the unique subset of contiguous read sub-locations
    read_locs = set()
    for r in chain(*paired_reads):
        for chrm, strand, start, stop in iter_coverage_regions_for_read( 
                r, reads, reverse_read_strand, pairs_are_opp_strand):
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
                read, reads, reverse_read_strand, pairs_are_opp_strand):
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
