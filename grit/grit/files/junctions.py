# Copyright (c) 2011-2012 Nathan Boley

import sys
import numpy
import time

from itertools import product, izip
import re
import itertools
from collections import defaultdict, namedtuple

import multiprocessing

from reads import get_strand, get_contigs_and_lens

from grit import config

CONSENSUS_PLUS = 'GTAG'
CONSENSUS_MINUS = 'CTAC'
                                  

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

def iter_jns_in_read( read ):
    """Iter junctions in read.

    Returns 0-based closed-closed intron coordiantes. 
    """

    # quickly check if read could span a single intron
    if len( read.cigar ) < 3:
        return
    
    # find all of the intron indices
    intron_indices = [i for i, (contig_type, length) in enumerate(read.cigar)
                      if contig_type == 3]
    
    # return if any of the junctions are too short
    if any(read.cigar[i][1] < config.MIN_INTRON_SIZE for i in intron_indices):
        return

    # forbid jucntion calling in reads with ref insertions or deletions
    if any(code in (1,2) for code, length in read.cigar):
        return
    
    # only accept reads with exactly 1 junction
    if len(intron_indices) != 1: 
        return
    
    intron_index = intron_indices[0]
    # iterate thorough all of the junctions
    for intron_index in intron_indices:
        # if the intron is at the beggining or end of the read, that doesn't 
        # make any sense, so skip this read
        if intron_index == 0: return
        if intron_index == (len(read.cigar)-1): return
        
        # skip introns that aren't flanked by a reference match
        if read.cigar[intron_index-1][0] != 0: continue
        if read.cigar[intron_index+1][0] != 0: continue

        # skip introns whose reference match is too short
        if read.cigar[intron_index-1][1] < config.MIN_INTRON_FLANKING_SIZE: 
            continue
        if read.cigar[intron_index+1][1] < config.MIN_INTRON_FLANKING_SIZE: 
            continue

        # Find the start base of the intron. We need to add all of the 
        # match and skip bases, and delete the deletions
        n_pre_intron_bases = 0
        for code, size in read.cigar[:intron_index]:
            if code == 0: n_pre_intron_bases += size
            elif code == 1: pass #n_pre_intron_bases += size
            elif code == 2: n_pre_intron_bases += size
            elif code == 3: n_pre_intron_bases += size
            #elif code == 4: n_pre_intron_bases += size
        
        upstrm_intron_pos = read.pos + n_pre_intron_bases
        dnstrm_intron_pos = upstrm_intron_pos + read.cigar[intron_index][1] - 1
        
        yield (upstrm_intron_pos, dnstrm_intron_pos)
        
    return

def extract_junctions_in_region( reads, chrm, strand, start=None, end=None, 
                                 allow_introns_to_span_start=False,
                                 allow_introns_to_span_end=False,
                                 only_unique=False ):
    reads.reload()
    all_junctions = defaultdict(lambda: defaultdict(int))
    for i, read in enumerate(reads.iter_reads(chrm, strand, start, end)):
        # check for uniqueness, if possible
        try: 
            if only_unique and int(read.opt('NH')) > 1: continue
        except KeyError: 
            pass
        
        for upstrm_intron_pos, dnstrm_intron_pos in iter_jns_in_read( read ):
            assert ( upstrm_intron_pos - dnstrm_intron_pos + 1 
                     < config.MIN_INTRON_SIZE )

            # Filter out reads that aren't fully in the region
            if start != None:
                if dnstrm_intron_pos < start: continue
                if not allow_introns_to_span_start and upstrm_intron_pos<start:
                    continue

            if end != None:
                if upstrm_intron_pos > end: continue
                if not allow_introns_to_span_end and dnstrm_intron_pos>end:
                    continue

            # increment count of junction reads at this read position
            # for this intron or initialize it to 1
            all_junctions[(upstrm_intron_pos, dnstrm_intron_pos)][read.pos] += 1
    
    rv = []
    for jn, cnts in sorted(all_junctions.iteritems()):
        cnts = numpy.array(cnts.values(), dtype=float)
        ps = cnts/cnts.sum()
        entropy = max(0, -float((ps*numpy.log2(ps)).sum()))
        rv.append((jn, int(cnts.sum()), entropy))
    
    return rv

def extract_junctions_in_contig( reads, chrm, strand ):
    return extract_junctions_in_region( 
        reads, chrm, strand, start=None, end=None )

def load_junctions_worker(all_jns, all_jns_lock, 
                          segments_queue, segments_queue_lock, reads):
    jns = defaultdict(list)
    while len(segments_queue) > 0:
        with segments_queue_lock:
            if len(segments_queue) == 0: break
            chrm, strand, start, stop = segments_queue.pop()
        if config.VERBOSE: 
            config.log_statement("Finding jns in '%s:%s:%i:%i'" % 
                          (chrm, strand, start, stop))
        jns[(chrm, strand)].extend(
            extract_junctions_in_region(
                reads, chrm, strand, start, stop, True))
    
    # finally, block until we can offload the remaining junctions
    with all_jns_lock:
        for key, region_jns in jns.iteritems():
            if key not in all_jns: all_jns_key = []
            else: all_jns_key = all_jns[key]
            all_jns_key.extend( region_jns )
            all_jns[key] = all_jns_key
    del jns
    if config.VERBOSE: config.log_statement( "" )
    return

def load_junctions_in_bam( reads, regions=None, nthreads=1):
    if regions == None:
        regions = []
        for contig, contig_len in zip(*get_contigs_and_lens([reads,])):
            for strand in '+-':
                regions.append( (contig, strand, 0, contig_len) )
    
    if nthreads == 1:
        jns = defaultdict(list)
        for chrm, strand, region_start, region_stop in regions:
            jns[(chrm, strand)].extend( extract_junctions_in_region( 
                    reads, chrm, strand, region_start, region_stop ) )
        return jns
    else:
        from multiprocessing import Process, Manager
        manager = Manager()
        all_jns = manager.dict()
        all_jns_lock = multiprocessing.Lock()
        
        segments_queue = manager.list()
        segments_queue_lock = multiprocessing.Lock()
        
        for chrm, strand, region_start, region_stop in regions:
            # add all the regions to search for junctions in
            seg_len = min(5000, int((region_stop - region_start + 1)/nthreads))
            pos = region_start
            while pos < region_stop:
                segments_queue.append( (chrm, strand, pos, pos+seg_len) )
                pos += seg_len
            # make sure the last region doesnt exten past the stop
            segments_queue[-1] = (
                chrm, strand, segments_queue[-1][2], region_stop)
        
        ps = []
        for i in xrange(nthreads):
            p = Process(target=load_junctions_worker,
                        args=( all_jns, all_jns_lock, 
                               segments_queue, segments_queue_lock, reads))
            
            p.start()
            ps.append( p )

        if config.VERBOSE:
            config.log_statement( "Waiting on jn finding children" )
        while len(segments_queue) > 0:
            if config.VERBOSE:
                config.log_statement( 
                    "Waiting on jn finding children (%i in queue)" 
                    % len(segments_queue) )
            time.sleep( 0.5 )

        if config.VERBOSE:
            config.log_statement("Waiting on jn finding children (0 in queue)")
        for p in ps: p.join()
        #while any( not p.is_alive() for p in ps ):

        if config.VERBOSE:
            config.log_statement("Merging junctions from threads")
        junctions = {}
        for key in all_jns.keys():
            junctions[key] = sorted(all_jns[key])
        return junctions
    assert False
