
# Copyright (c) 2011-2012 Nathan Boley

import sys
import numpy

from itertools import product, izip
import re

import itertools
from collections import defaultdict, namedtuple

from reads import get_strand

VERBOSE = False

CONSENSUS_PLUS = 'GTAG'
CONSENSUS_MINUS = 'CTAC'

MIN_INTRON_LEN = 20

def get_jn_type( chrm, upstrm_intron_pos, dnstrm_intron_pos, 
                 fasta, jn_strand="UNKNOWN" ):
    # get first 2 bases from 5' and 3' ends of intron to determine 
    # splice junction strand by comparing to consensus seqs
    # subtract one from start since fasta is 0-based closed-open
    intron_seq = \
        fasta.fetch( chrm , upstrm_intron_pos-1, upstrm_intron_pos+1) + \
        fasta.fetch( chrm , dnstrm_intron_pos-2, dnstrm_intron_pos )
    
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

def extract_junctions_in_region( reads, chrm, strand, start=None, end=None, 
                                 allow_introns_to_span_start=False,
                                 allow_introns_to_span_end=False):
    reads.reload()
    all_junctions = defaultdict(int)
    for i, read in enumerate(reads.fetch(chrm, start, end)):
        # increment the number of times we've seen this read
        if not read_spans_single_intron( read ):
            continue

        if strand != get_strand( read, reads.reverse_read_strand, 
                                 reads.pairs_are_opp_strand ):
            continue
        
        # find the introns
        gaps = [ (i, size) for i, (code, size) in enumerate(read.cigar)
                 if code == 3 and size > MIN_INTRON_LEN]
        
        # skip reads without exactly 1 substatnial gap
        if len( gaps ) != 1: continue        
        intron_index, intron_len = gaps[0]
        
        # Find the start base of the intron. We need to add all of the 
        # match and skip bases, and delete the deletions
        n_pre_intron_bases = 0
        for code, size in read.cigar[:intron_index]:
            if code == 0: n_pre_intron_bases += size
            elif code == 1: n_pre_intron_bases -= size
            elif code == 2: n_pre_intron_bases += size
            elif code == 3: n_pre_intron_bases += size
            elif code == 4: n_pre_intron_bases += size
        
        # add one to left_intron since bam files are 0-based
        upstrm_intron_pos = read.pos + n_pre_intron_bases + 1
        dnstrm_intron_pos = upstrm_intron_pos + read.cigar[intron_index][1] - 1
        
        assert upstrm_intron_pos - dnstrm_intron_pos + 1 < MIN_INTRON_LEN
        
        # Filter out reads that aren't fully in the region
        if start != None:
            if dnstrm_intron_pos < start: continue
            if False == allow_introns_to_span_start and upstrm_intron_pos < start:
                continue
        
        if end != None:
            if upstrm_intron_pos > end: continue
            if False == allow_introns_to_span_end and dnstrm_intron_pos > end:
                continue
        
        # increment count of junction reads at this read position
        # for this intron or initialize it to 1
        all_junctions[ (upstrm_intron_pos, dnstrm_intron_pos) ] += 1
    
    return sorted( all_junctions.iteritems() )

def extract_junctions_in_contig( reads, chrm, strand ):
    return extract_junctions_in_region( 
        reads, chrm, strand, start=None, end=None )
