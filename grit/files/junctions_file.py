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
import numpy

from itertools import product, izip
import re

from reads import GenomicIntervals, GenomicInterval
import itertools
from collections import defaultdict, namedtuple

VERBOSE = False

CONSENSUS_PLUS = 'GTAG'
CONSENSUS_MINUS = 'CTAC'

def get_jn_type( chrm, upstrm_intron_pos, dnstrm_intron_pos, 
                 fasta, jn_strand="UNKNOWN" ):
    # get first 2 bases from 5' and 3' ends of intron to determine 
    # splice junction strand by comparing to consensus seqs
    # subtract one from start since fasta is 0-based closed-open
    intron_seq = \
        fasta.fetch( 'chr'+chrm , upstrm_intron_pos-1, upstrm_intron_pos+1) + \
        fasta.fetch( 'chr'+chrm , dnstrm_intron_pos-2, dnstrm_intron_pos )
    
    # if junction matchs consensus set strand
    # else return non-canonical
    if intron_seq.upper() == CONSENSUS_PLUS:
        canonical_strand = '+'
    elif intron_seq.upper() == CONSENSUS_MINUS:
        canonical_strand = '-'
    else:
        if jn_strand == "UNKNOWN":
            return 'non-canonical', '.'
        else:
            return 'non-canonical'
    
    # if we don't know what the strand should be, then use the canonical seq
    # to infer the strand of the jn
    if jn_strand == "UNKNOWN":
        return 'canonical', canonical_strand
    # otherwise, if we know the jn's strand
    else:
        if jn_strand == canonical_strand:
            return 'canonical'
        return 'canonical_wrong_strand'
    
    assert False

_junction_named_tuple_slots = [
    "region", "type", "cnt", "uniq_cnt", "source_read_offset", "source_id" ]
    
_JnNamedTuple = namedtuple( "Junction", _junction_named_tuple_slots )
                 
class Junction( _JnNamedTuple ):
    valid_jn_types = set(( "infer", "canonical", 
                           "canonical_wrong_strand", "non_canonical" ))
    
    def __new__( self, region,
                 jn_type=None, cnt=None, uniq_cnt=None,
                 source_read_offset=None, source_id=None ):
        # do type checking
        #if not isinstance( region, GenomicInterval ):
        #    raise ValueError, "regions must be of type GenomicInterval"
        if region.strand not in ("+", "-" ):
            raise ValueError, "Unrecognized strand '%s'" % strand
        if jn_type != None and jn_type not in self.valid_jn_types:
            raise ValueError, "Unrecognized jn type '%s'" % jn_type
        
        if cnt != None: cnt = int( cnt )
        if uniq_cnt != None: uniq_cnt = int( uniq_cnt )
        
        return _JnNamedTuple.__new__( 
            Junction, region,
            jn_type, cnt, uniq_cnt, source_read_offset, source_id )
    
    chrm = property( lambda self: self.region.chr )
    strand = property( lambda self: self.region.strand )
    start = property( lambda self: self.region.start )
    stop = property( lambda self: self.region.stop )

    def build_gff_line( self, group_id=None, fasta_obj=None ):
        if self.type == None and fasta_obj != None:
            intron_type = get_jn_type( 
                self.chr, self.start, self.stop, fasta_obj, self.strand )
        else:
            intron_type = self.type
        
        group_id_str = str(group_id) if group_id != None else ""

        if self.source_read_offset != None:
            group_id_str += ' source_read_offset "{0}";'.format( 
                self.source_read_offset )
        
        if self.uniq_cnt != None:
            group_id_str += ' uniq_cnt "{0}";'.format( self.uniq_cnt )
        
        if intron_type != None:
            group_id_str += ' type "{0}";'.format( intron_type )
        
        count = self.cnt if self.cnt != None else 0
        
        return create_gff_line( 
            self.region, group_id_str, score=count, feature='intron' )


def read_spans_single_intron( read ):
    # quickly check if read could spans a single intron
    if len( read.cigar ) != 3:
        return False
    
    # check that cigar regions are alternating matching(exon) and \
    # skipping(intron) regions
    for i, cigar_region in enumerate(read.cigar):
        if (i % 2 == 0 and cigar_region[0] != 0) or \
                (i % 2 == 1 and cigar_region[0] != 3):
            return False
    
    # if the read is not primary then do not include it in the counts
    if read.is_secondary:
        return False
    
    return True

def extract_junctions_in_contig( reads, chrm, strand, 
                                 reverse_strand, pairs_are_opp_strand ):
    # store how many times we have observed each read
    query_name_cnts = defaultdict( int )
    jn_reads = []
    for read in reads.fetch(chrm):
        # increment the number of times we've seen this read
        query_name_cnts[ read.qname ] += 1
        if read_spans_single_intron( read ):
            jn_reads.append( read )

    for read in jn_reads:
        if strand != get_strand( read, reverse_strand, pairs_are_opp_strand ):
            continue
        
        # add one to left_intron since bam files are 0-based
        upstrm_intron_pos = read.pos + read.cigar[0][1] + 1
        dnstrm_intron_pos = upstrm_intron_pos + read.cigar[1][1] - 1
        
        # increment count of junction reads at this read position
        # for this intron or initialize it to 1
        all_junctions[ (upstrm_intron_pos, dnstrm_intron_pos) ] += 1
    
    return sorted( all_junctions.iteritems() )

def iter_junctions( reads_fn, fasta_fn, stranded, reverse_strand, pairs_are_opp_strand ):
    """Get all of the junctions represented in a reads object
    """
    
    # build reads object
    reads = Samfile( reads_fn, "rb" )
    fasta_obj = None if fasta_fn == None else Fastafile( fasta_fn )
    
    # store the number of reads across the junction for each relative 
    # read position
    all_junctions = defaultdict(int)
    unique_junctions = defaultdict(int)
    
    for ref in reads.references:
        find_junctions_in_chrm( ref, all_junctions, unique_junctions )
    
    jns = []
    read_offsets = []
    cnts = []
    uniq_cnts = []
    for (chrm, strand, start, stop, read_offset ), cnt \
            in sorted(all_junctions.iteritems()):

        uniq_cnt = unique_junctions[ (chrm, strand, start, stop, read_offset ) ]
        jn = Junction( GenomicInterval(chrm, strand, start, stop), 
                       source_read_offset=read_offset, 
                       cnt=cnt, uniq_cnt=uniq_cnt )
        yield jn
        
    return

