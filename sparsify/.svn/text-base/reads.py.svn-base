import pysam
import numpy
import os

from itertools import product
from collections import namedtuple, defaultdict

import gene_models

PRINT_READ_PAIRS = False
DEBUG_VERBOSE = False
MINIMUM_INTRON_SIZE = 20
VERBOSE = False
DEBUG = True

GenomicInterval = namedtuple('GenomicInterval', ['chr', 'strand', 'start', 'stop'])
# Intron = namedtuple('GenomicInterval', ['chr', 'strand', 'start', 'stop'])

def is_junction_read( read ):
    return ( read.qlen != read.aend - read.pos )

def get_bndrys_from_jn_read( read  ):
    """Given a junction read, return the end and start of the exons it comes from.

    """
    if len( read.cigar ) != 3:
        if DEBUG_VERBOSE:
            print "Weird Cigar: ", read.cigar
        return None

    if read.cigar[1][1] < MINIMUM_INTRON_SIZE:
        if DEBUG_VERBOSE:
            print "Weird Cigar: ", read.cigar
        return None
    
    # I SHOULD HAVE TO add 1 to the start to make it 1 based
    # according to the pysam docs, but this doesn't seem to actually
    # be true. Dont know where the bug is.
    start = read.pos + read.cigar[0][1] # +1
    # I add 1 to this to make it 1 based ( and it works in this case )
    stop = read.aend - read.cigar[2][1] + 1
    try:
        assert( stop == read.pos + read.cigar[0][1] + read.cigar[1][1] + 1 )
    except:
        print stop, read.pos, read.cigar
        raise
    
    return (start, stop)


class Reads( pysam.Samfile ):
    """Subclass the samfile object to include a method that returns reads 
       and their pairs.


    """
    def iter_paired_reads( self, gene_bndry ):
        # whether or not the gene is on the positive strand
        gene_strnd_is_rev = ( gene_bndry.strand == '-' )
        
        # get all of the first pairs
        def iter_pair1_reads():
            for read in self.fetch( \
                gene_bndry.chr, gene_bndry.start, gene_bndry.stop  ):
                if read.is_read1 \
                    and read.is_reverse == gene_strnd_is_rev:
                    yield read
        
        # index the pair 2 reads
        reads_pair2 = {}
        for read in self.fetch( gene_bndry.chr, gene_bndry.start, gene_bndry.stop  ):
            if not read.is_read1 \
                and read.is_reverse != gene_strnd_is_rev:
                reads_pair2[ read.qname ] = read
        
        # iterate through the read pairs
        for read1 in iter_pair1_reads():
            try:
                read2 = reads_pair2[ read1.qname ]
            # if there is no mate, skip this read
            except KeyError:
                if DEBUG_VERBOSE:
                    print "No mate: ", read1.pos, read1.aend, read1.qname
                continue
            
            assert ( read1.qlen == read1.aend - read1.pos ) \
                   or ( len( read1.cigar ) > 1 )
            assert ( read2.qlen == read2.aend - read2.pos ) \
                   or ( len( read2.cigar ) > 1 )
            
            if read1.qlen != read2.qlen:
                print( "ERROR: unequal read lengths %i and %i\n", read1.qlen, read2.qlen )
                continue
            
            yield read1, read2
        
        return

def build_nonoverlapping_bin_observations( f_mat ):
    """Build single bin counts from f_mat dictionary.
    
    """
    def iter_single_bins( pseudo_bin ):
        for item in pseudo_bin[0]:
            yield item
        for jn in zip( pseudo_bin[0][:-1], pseudo_bin[0][1:] ):
            yield jn
        for item in pseudo_bin[1]:
            yield item
        for jn in zip( pseudo_bin[1][:-1], pseudo_bin[1][1:] ):
            yield jn
        
    single_bins = defaultdict( lambda: [ numpy.zeros(len(f_mat.values()[0][0])), 0 ] )
    for pseudo_bin, ( exp, obs ) in f_mat.iteritems():
        for single_bin in iter_single_bins( pseudo_bin[2] ):
            single_bins[ single_bin ][0] += exp
            single_bins[ single_bin ][1] += obs

    """
    for key, val in sorted(single_bins.iteritems()):
        print str(key).ljust( 20), val
    """

    return single_bins

## Associate reads with bins
class BinnedReads( object ):
    def _find_nonoverlapping_exons_covered_by_segment( self, start, stop ):
        """Return the pseudo bins that a given segment has at least one basepair in.

        """
        bin_1 = self.nonoverlapping_exon_bndrys.searchsorted( start, side='right' ) - 1
        # if the start falls before all bins
        if bin_1 == -1: return ()

        bin_2 = self.nonoverlapping_exon_bndrys.searchsorted( stop, side='right' ) - 1
        # if the stop falls after all bins
        if bin_2 == len( self.nonoverlapping_exon_bndrys ) - 1: return ()
        
        if DEBUG:
            assert bin_1 == -1 or start >= self.nonoverlapping_exon_bndrys[ bin_1  ]
            assert stop < self.nonoverlapping_exon_bndrys[ bin_2 + 1  ]
        
        return tuple(xrange( bin_1, bin_2+1 ))
    
    @staticmethod
    def _iter_contiguous_segments_in_read( read ):
        """Iter through contiguous regions in a read. For a simple read, this
           is simply the entire segment. For a read with one junction, there 
           are 2 segments. Etc.
        """
        contiguous_segments = []
        # make this one based
        start = read.pos+1
        for region_type, region_size in read.cigar:
            # If this is a 'match' region
            if region_type == 0:
                # subtract one to make the interval closed closed
                contiguous_segments.append( ( start, start-1+region_size ) )
                start += region_size
            # if this is a skipped region ( ie, an intron ), then we need to
            # increment the start.
            elif region_type == 3:
                start += region_size
            # if it's another region type, skip this read
            else:
                return []
                
        return contiguous_segments

    def _find_nonoverlapping_exons_covered_by_read( self, read ):
        # fast path for non-junction reads
        if len( read.cigar ) == 1:
            return self._find_nonoverlapping_exons_covered_by_segment( read.pos, read.aend )
        # if this is a junction read
        else:
            all_covered_nonoverlapping_exons = []
            contiguous_segments = list( self._iter_contiguous_segments_in_read( read ) )
            for index, segment in enumerate( contiguous_segments ):
                covered_nonoverlapping_exons = \
                    self._find_nonoverlapping_exons_covered_by_segment( *segment )
                
                if len( covered_nonoverlapping_exons ) == 0: return ()
                
                # if this is an internal splice, make sure that the end
                # corresponds to an exon boundary
                if index < len(contiguous_segments) - 1:
                    # we add one to the segment, because the current segment gives us the last
                    # covered base in the read, and we want the start ( inclusive ) of the 
                    # next segment.
                    if segment[1]+1 != self.nonoverlapping_exon_bndrys[ \
                            covered_nonoverlapping_exons[-1] + 1 ]:
                        return ()
                # is this isn't the first splice, make sure that the start corresponds
                # to a bin end.
                if index > 0:
                    # we add one to the segment, because the current segment gives us the last
                    # covered base in the read, and we want to the start ( inclusive ) of the 
                    # next segment
                    if segment[0] != self.nonoverlapping_exon_bndrys[ 
                            covered_nonoverlapping_exons[0] ]:
                        return ()
                    
                all_covered_nonoverlapping_exons.extend( covered_nonoverlapping_exons )
                
            return tuple( all_covered_nonoverlapping_exons )
    
    @staticmethod
    def get_read_group( r1, r2 ):        
        r1_read_group = [ val for key, val in r1.tags if key == 'RG' ]
        r1_read_group = r1_read_group[0] if len( r1_read_group ) == 1 else 'mean'
        r2_read_group = [ val for key, val in r2.tags if key == 'RG' ]
        r2_read_group = r2_read_group[0] if len( r2_read_group ) == 1 else 'mean'
        if r1_read_group == r2_read_group:
            return r1_read_group
        else: 
            print "WARNING: Read groups do not match."
            print r1
            print r2
            return None
    
    def bin_reads( self, read_group_mappings ):
        # actually bin the reads
        binned_reads = defaultdict(int)
        for r1, r2 in self.reads.iter_paired_reads( self.gene.boundaries ):
            # find the read type
            read_group = self.get_read_group( r1, r2 )
            if read_group == None: continue
            if read_group_mappings != None:
                read_group = read_group_mappings[ read_group ]
            
            read_len = r1.rlen
            if r1.rlen != r2.rlen:
                raise NotImplementedError, "Read lens of paired reads must be the same."
            
            pseudo_bin = ( \
                self._find_nonoverlapping_exons_covered_by_read( r1 ), \
                self._find_nonoverlapping_exons_covered_by_read( r2 ) \
            )
                        
            # if one of the ends didnt fall in any bins
            if () in pseudo_bin:
                continue
            
            # because read 2 is always before read 1 in the stranded protool,
            # and the bins are always in genome direction, we have to switch the bins.
            if self.gene.strand == '-':
                pseudo_bin = ( pseudo_bin[1], pseudo_bin[0] )
            
            if pseudo_bin[0][0] > pseudo_bin[1][0]:
                # print "WARNING: a read has been observed with the pairs in the wrong order ( ie, in a negative stranded gene, the first read was before the second read in genomic coordinates )"
                continue
            
            binned_reads[ ( read_len, read_group, pseudo_bin ) ] += 1
        
        if DEBUG_VERBOSE:
            for key, val in binned_reads.iteritems():
                print str(key).ljust( 100 ), val
        
        return binned_reads
    
    def find_connected_exons( self ):
        """For now we just look at junctions.

        """
        # get all of the possible junctions from the exon boundaries
        # note that there is a non-overlapping criteria, so this is more 
        # complicated than it looks.
        junctions = list( self.gene.iter_possible_junctions() )
        # build a set of the junctions we can map to for faster lookups
        jns_bins_map = defaultdict( list )
        for i, j in junctions:
            jns_bins_map[ self.gene.exon_bndrys[i][1], self.gene.exon_bndrys[j][0] ].append( (i,j) )
        
        def find_possible_junctions( read ):
            loc = get_bndrys_from_jn_read( read )
            try:
                return jns_bins_map[loc]
            # if this doesn't correspond to a real junction, or
            # the read isn't actually a junction ( like the multi
            # junction reads ) then return an empty list
            except KeyError:
                if DEBUG_VERBOSE:
                    print "Key Error: ", loc
                return []

        connected_exons = set()
        for r1, r2 in self.reads.iter_paired_reads( self.gene.boundaries ):
            if is_junction_read( r1 ):
                connected_exons.update( find_possible_junctions( r1 ) )
            if is_junction_read( r2 ):
                connected_exons.update( find_possible_junctions( r2 ) )
        
        return connected_exons
    
    def find_connected_exons_read_coverage( self ):
        connected_exons_coverage = defaultdict( int )
        for ( rl, cluster_grp, bin ), cnt in self.binned_reads.iteritems():
            for read_bins in bin:
                for exon in read_bins:
                    connected_exons_coverage[ exon ] += ( float(cnt)/len(read_bins) )
        
        return connected_exons_coverage
    
    def __init__( self, gene, reads, read_group_mappings=None ):
        assert isinstance( gene, gene_models.GeneModel )
        assert isinstance( reads, Reads )
        
        self.sourcefile = os.path.basename( reads.filename )
        self.reads = reads
        self.gene = gene
                                      
        self.nonoverlapping_exon_bndrys = numpy.array( gene.nonoverlapping_exon_boundaries )
        
        self.binned_reads = self.bin_reads( read_group_mappings )
        self.total_num_reads = float( sum( self.binned_reads.values() ) )
        
        self.connected_exon_pairs = sorted( self.find_connected_exons() )

        read_groups = defaultdict( int )
        for (read_len, read_group, possible_bins), cnt \
                in self.binned_reads.iteritems():
            read_groups[ ( read_len, read_group ) ] += cnt
        
        self.read_groups = sorted( read_groups.keys() )
        self.read_group_cnts = read_groups
                
        return
