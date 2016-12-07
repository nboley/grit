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

import os
import time
import numpy
import scipy
import math

from collections import defaultdict, namedtuple
from itertools import chain

import multiprocessing
import queue

from scipy.stats import beta, binom

import networkx as nx

from copy import copy

ReadCounts = namedtuple('ReadCounts', ['Promoters', 'RNASeq', 'Polya'])

from .frag_len import build_fl_dists_from_fls_dict

from .transcript import Transcript, Gene

from .files.junctions import filter_jns

from . import config

from .files.reads import MergedReads, RNAseqReads, CAGEReads, \
    RAMPAGEReads, PolyAReads, \
    fix_chrm_name_for_ucsc, get_contigs_and_lens, calc_frag_len_from_read_data, \
    iter_paired_reads, extract_jns_and_reads_in_region, TooManyReadsError
from . import files.junctions

from .files.bed import create_bed_line

class Bin(object):
    start = None
    stop = None
    
    def reverse_strand(self, contig_len):
        kwargs = copy(self.__dict__)
        kwargs['start'] = contig_len-1-self.stop
        kwargs['stop'] = contig_len-1-self.start
        if 'left_labels' in kwargs:
            (kwargs['left_labels'], kwargs['right_labels'] 
             ) = (kwargs['right_labels'], kwargs['left_labels'])
        return type(self)(**kwargs)
    
    def length(self):
        return self.stop - self.start + 1
    
    def mean_cov( self, cov_array ):
        return numpy.median(cov_array[self.start:self.stop+1])
        return cov_array[self.start:self.stop].mean()
    
    #def reverse_strand(self, contig_len):
    #    return Bin(contig_len-1-self.stop, contig_len-1-self.start, 
    #               self.right_label, self.left_label, self.type)

    #def reverse_coords(self, contig_len):
    #    return Bin(contig_len-1-self.stop, contig_len-1-self.start, 
    #               self.left_label, self.right_label, self.type)
        
    _bndry_color_mapping = {
        'CONTIG_BNDRY': '0,0,0',
        'GENE_BNDRY': '0,0,0',
        
        'POLYA': '255,255,0',

        'CAGE_PEAK': '0,255,0',
        
        'D_JN': '173,255,47',
        'R_JN': '0,0,255',
        
        'ESTART': '0,0,0',
        'ESTOP': '0,0,0'
    }
    
    def find_bndry_color( self, bndry ):
        return self._bndry_color_mapping[ bndry ]
    
    def _find_colors( self, strand ):
        if self.type != None:
            if self.type =='GENE':
                return '0,0,0'
            if self.type =='CAGE_PEAK':
                return '0,0,0'
            if self.type =='EXON':
                return '0,0,0'
            if self.type =='EXON_EXT':
                return '0,0,255'
            if self.type =='RETAINED_INTRON':
                return '255,255,0'
            if self.type =='TES_EXON':
                return '255,0,0'
            if self.type =='TSS_EXON':
                return '0,255,0'
            if self.type =='SE_GENE':
                return '255,255,0'
            if self.type =='INTERGENIC_SPACE':
                return '254,254,34'
        
        if strand == '+':
            left_label, right_label = self.left_label, self.right_label
        else:
            assert strand == '-'
            left_label, right_label = self.right_label, self.left_label
        
        
        if left_label == 'D_JN' and right_label  == 'R_JN':
            return '108,108,108'
        if left_label == 'D_JN' and right_label  == 'D_JN':
            return '135,206,250'
        if left_label == 'R_JN' and right_label  == 'R_JN':
            return '135,206,250'
        if left_label == 'R_JN' and right_label  == 'D_JN':
            return '0,0,255'
        if left_label == 'R_JN' and right_label  == 'POLYA':
            return '255,0,0'
        if left_label == 'POLYA' and right_label  == 'POLYA':
            return ' 240,128,128'
        if left_label == 'D_JN' and right_label  == 'POLYA':
            return '240,128,128'
        if left_label == 'POLYA' and right_label  == 'D_JN':
            return '147,112,219'
        if left_label == 'POLYA' and right_label  == 'R_JN':
            return '159,153,87'
        if left_label == 'ESTART' and right_label  == 'ESTOP':
            return '159,153,87'
        
        return ( self.find_bndry_color(left_label), 
                 self.find_bndry_color(right_label) )

class SegmentBin(Bin):
    def __init__(self, start, stop, left_labels, right_labels,
                 type=None, 
                 fpkm_lb=None, fpkm=None, fpkm_ub=None,
                 cnt=None):
        self.start = start
        self.stop = stop
        assert stop - start >= 0

        self.left_labels = sorted(left_labels)
        self.right_labels = sorted(right_labels)

        self.type = type
        
        self.fpkm_lb = fpkm_lb
        self.fpkm = fpkm
        self.fpkm_ub = fpkm_ub
        
        self.cnt = cnt
    
    def set_expression(self, fpkm_lb, fpkm, fpkm_ub):
        assert not any (math.isnan(x) for x in (fpkm_lb, fpkm, fpkm_ub))
        self.fpkm_lb = fpkm_lb
        self.fpkm = fpkm
        self.fpkm_ub = fpkm_ub
        return self
    
    def __repr__( self ):
        type_str =  ( 
            "(%s-%s)" % ( ",".join(self.left_labels), 
                        ",".join(self.right_labels) ) 
            if self.type == None else self.type )
        loc_str = "%i-%i" % ( self.start, self.stop )
        rv = "%s:%s" % (type_str, loc_str)
        if self.fpkm != None:
            rv += ":%.2f-%.2f TPM" % (self.fpkm_lb, self.fpkm_ub)
        return rv

class TranscriptBoundaryBin(SegmentBin):
    def __init__(self, start, stop, left_labels, right_labels,
                 type, cnts):
        assert type in ('CAGE_PEAK', 'POLYA')
        SegmentBin.__init__(
            self, start=start, stop=stop, 
            left_labels=left_labels, right_labels=right_labels,
            type=type)
        self.cnts = cnts
    
    def set_tpm(self, total_num_reads):
        SegmentBin.set_expression(
            self, 
            1e6*beta.ppf(0.01, self.cnts.sum()+1e-6, total_num_reads+1e-6),
            1e6*beta.ppf(0.50, self.cnts.sum()+1e-6, total_num_reads+1e-6),
            1e6*beta.ppf(0.99, self.cnts.sum()+1e-6, total_num_reads+1e-6))
        return self

class  TranscriptElement( Bin ):
    def __init__( self, start, stop, type, fpkm):
        self.start = start
        self.stop = stop
        assert stop - start >= 0
        self.type = type
        self.fpkm = fpkm
    
    def length( self ):
        return self.stop - self.start + 1
        
    def __repr__( self ):
        type_str = self.type
        loc_str = "%i-%i" % ( self.start, self.stop )
        rv = "%s:%s" % (type_str, loc_str)
        if self.fpkm != None:
            rv += ":%.2f FPKM" % self.fpkm
        return rv

class GeneElements(object):
    def __init__( self, chrm, strand  ):
        self.chrm = chrm
        self.strand = strand
        
        self.regions = []
        self.element_segments = []
        self.elements = []

    def find_coverage(self, reads):
        cov = numpy.zeros(self.stop-self.start+1, dtype=float)
        for x in self.regions:
            seg_cov = reads.build_read_coverage_array( 
                self.chrm, self.strand, x.start, x.stop )
            cov[x.start-self.start:x.stop-self.start+1] = seg_cov
        #if gene.strand == '-': cov = cov[::-1]
        return cov
    
    def base_is_in_gene(self, pos):
        return all( r.start <= pos <= r.stop for r in self.regions )
    
    def reverse_strand( self, contig_len ):
        rev_gene = GeneElements( self.chrm, self.strand )
        for bin in reversed(self._regions):
            rev_gene._regions.append( bin.reverse_strand( contig_len ) )
        for bin in reversed(self._element_segments):
            rev_gene._element_segments.append( bin.reverse_strand( contig_len ) )
        for bin in reversed(self._exons):
            rev_gene._exons.append( bin.reverse_strand( contig_len ) )
        
        return rev_gene

    def shift( self, shift_amnt ):
        rev_gene = GeneElements( self.chrm, self.strand )
        for bin in self._regions:
            rev_gene._regions.append( bin.shift( shift_amnt ) )
        for bin in self._element_segments:
            rev_gene._element_segments.append( bin.shift( shift_amnt ) )
        for bin in self._exons:
            rev_gene._exons.append( bin.shift( shift_amnt ) )
        
        return rev_gene
        
    @property
    def start(self):
        return min( x.start for x in self.regions )

    @property
    def stop(self):
        return max( x.stop for x in self.regions )

    def __repr__( self ):
        loc_str = "GENE:%s:%s:%i-%i" % ( 
            self.chrm, self.strand, self.start, self.stop )
        return loc_str

    def write_elements_bed( self, ofp ):
        feature_mapping = { 
            'GENE': 'gene',
            'CAGE_PEAK': 'promoter',
            'SE_GENE': 'single_exon_gene',
            'TSS_EXON': 'tss_exon',
            'EXON': 'internal_exon',
            'TES_EXON': 'tes_exon',
            'INTRON': 'intron',
            'POLYA': 'polya',
            'INTERGENIC_SPACE': 'intergenic',
            'RETAINED_INTRON': 'retained_intron',
            'UNKNOWN': 'UNKNOWN'
        }

        color_mapping = { 
            'GENE': '200,200,200',
            'CAGE_PEAK': '153,255,000',
            'SE_GENE': '000,000,200',
            'TSS_EXON': '140,195,59',
            'EXON': '000,000,000',
            'TES_EXON': '255,51,255',
            'INTRON': '100,100,100',
            'POLYA': '255,0,0',
            'INTERGENIC_SPACE': '254,254,34',
            'RETAINED_INTRON': '255,255,153',
            'UNKNOWN': '0,0,0'
        }

        chrm = self.chrm
        if config.FIX_CHRM_NAMES_FOR_UCSC:
            chrm = fix_chrm_name_for_ucsc(chrm)

        # write the gene line
        bed_line = create_bed_line( chrm, self.strand, 
                                    self.start, self.stop+1, 
                                    feature_mapping['GENE'],
                                    score=1000,
                                    color=color_mapping['GENE'],
                                    use_thick_lines=True,
                                    blocks=[(x.start, x.stop) for x in self.regions])
        ofp.write( bed_line + "\n"  )
        try: max_min_fpkm = max(1e-1, max(bin.fpkm for bin in self.elements)) \
           if len(self.elements) > 0 else 0
        except:  max_min_fpkm = 1000
        
        for element in self.elements:
            region = ( chrm, self.strand, element.start, element.stop)

            blocks = []
            use_thick_lines=(element.type != 'INTRON')
            element_type = element.type
            if element_type == None: 
                element_type = 'UNKNOWN'
                continue
            
            try: fpkm = element.fpkm_ub
            except: fpkm = element.fpkm
            score = min(1000, int(1000*fpkm/max_min_fpkm))
            score = 1000
            grp_id = element_type + "_%s_%s_%i_%i" % region

            # also, add 1 to stop because beds are open-closed ( which means no net 
            # change for the stop coordinate )
            bed_line = create_bed_line( chrm, self.strand, 
                                        element.start, element.stop+1, 
                                        feature_mapping[element_type],
                                        score=score,
                                        color=color_mapping[element_type],
                                        use_thick_lines=use_thick_lines,
                                        blocks=blocks)
            ofp.write( bed_line + "\n"  )

        return

    def writeGff( self, ofp ):
        """
            chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
        """
        if self.strand == '-':
            writetable_bins = self.reverse_strand( contig_len )
        else:
            writetable_bins = self
        
        for bin in writetable_bins:
            if filter != None and bin.type != filter:
                continue
            chrm = elements.chrm
            if config.FIX_CHRM_NAMES_FOR_UCSC:
                chrm = fix_chrm_name_for_ucsc(self.chrm)
            # add 1 because gffs are 1-based
            region = GenomicInterval(chrm, self.strand, 
                                     bin.start+1, bin.stop+1)
            grp_id = "%s_%s_%i_%i" % region
            ofp.write( create_gff_line(region, grp_id) + "\n" )
        
        return


def cluster_intron_connected_segments( segments, introns ):
    if len(segments) == 0:
        return []
    segments = sorted(segments)
    segment_starts = numpy.array([x[0] for x in segments])
    segment_stops = numpy.array([x[1] for x in segments])

    edges = set()
    for start, stop in introns:
        # Skip junctions that dont fall into any segment
        if start-1 < segment_starts[0]: continue
        if stop+1 >= segment_stops[-1]: continue
        
        # find which bin the segments fall into. Note, that since the 
        # segments don't necessarily tile the genome, it's possible
        # for the returned bin to not actually contain the junction
        start_bin = segment_starts.searchsorted( start-1, side='right' )-1
        assert start_bin >= 0
        stop_bin = segment_starts.searchsorted( stop+1, side='right' )-1

        # since the read coverage is determined in part determined by 
        # the junctions, we should never see a junction that doesn't fall
        # into a segment
        try:
            assert ( segment_starts[start_bin] <= 
                     start-1 <= segment_stops[start_bin] ), str([
                         segment_starts[start_bin],
                     start-1, segment_stops[start_bin],
                         segment_starts[start_bin+1],
                         start-1, segment_stops[start_bin+1]])

            assert ( segment_starts[stop_bin] <= 
                     stop+1 <= segment_stops[stop_bin]), str([
                         segment_starts[stop_bin],
                         stop-1, segment_stops[stop_bin],
                         segment_starts[stop_bin+1],
                         stop-1, segment_stops[stop_bin+1]])
        except:
            raise
            continue
        #if start > segment_stops[start_bin]: continue
        #if stop > segment_stops[stop_bin]: continue
        # XXX - dont rememeber why I was doing this
        #assert stop_bin < len(segment_starts)-1, \
        #    str([stop_bin, len(segment_stops), segment_stops[stop_bin]])
        if start_bin != stop_bin:
            edges.add((int(min(start_bin, stop_bin)), 
                       int(max(start_bin, stop_bin))))
    
    genes_graph = nx.Graph()
    genes_graph.add_nodes_from(range(len(segment_starts)))
    genes_graph.add_edges_from(edges)
    
    segments = []
    for g in nx.connected_components(genes_graph):
        g = sorted(g)
        segment = []
        prev_i = g[0]
        segment.append( [segment_starts[prev_i], ])
        for i in g[1:]:
            # if we've skipped at least one node, then add
            # a new region onto this segment
            if i > prev_i + 1:
                segment[-1].append( segment_stops[prev_i] )
                segment.append( [segment_starts[i], ])
                prev_i = i
            # otherwise, we've progressed to an adjacent sergments
            # so just merge the adjacent intervals
            else:
                assert i == prev_i + 1
                prev_i += 1
        segment[-1].append( segment_stops[g[-1]] )
        
        segments.append(segment)
    
    return segments

def find_empty_regions( cov, thresh=1e-6, 
                        min_length=config.MAX_EMPTY_REGION_SIZE ):
    x = numpy.diff( numpy.asarray( cov >= thresh, dtype=int ) )
    stops = (numpy.nonzero(x==1)[0]).tolist()
    if cov[-1] < thresh: stops.append(len(x))
    starts = (numpy.nonzero(x==-1)[0] + 1).tolist()
    if cov[0] < thresh: starts.insert(0, 0)
    assert len(starts) == len(stops)
    return [ x for x in zip(starts, stops) if x[1]-x[0]+1 >= min_length ]


def merge_adjacent_intervals(
        intervals, max_merge_distance=None):
    if len(intervals) == 0: return []
    intervals.sort()
    merged_intervals = [list(intervals[0]),]
    prev_stop = merged_intervals[-1][1]
    for start, stop in intervals[1:]:
        if start - max_merge_distance - 1 <= prev_stop:
            merged_intervals[-1][1] = stop
        else:
            merged_intervals.append([start, stop])
        prev_stop = stop
    return merged_intervals

def find_transcribed_regions( cov, thresh=1e-6 ):
    empty_regions = find_empty_regions(cov, thresh, 0)
    if len(empty_regions) == 0: 
        return [[0, len(cov)-1],]

    transcribed_regions = []
    if empty_regions[0][0] == 0:
        transcribed_regions.append([empty_regions.pop(0)[1]+1,])
    else:
        transcribed_regions.append([0,])
    for start, stop in empty_regions:
        transcribed_regions[-1].append(start-1)
        transcribed_regions.append([stop+1,])
    if transcribed_regions[-1][0] == len(cov):
        transcribed_regions.pop()
    else:
        transcribed_regions[-1].append(len(cov)-1)

    return transcribed_regions
    
    """
    # code to try and merge low signal segments, but I think that 
    # this is the wrong appraoch. I should be merging introns that 
    # could have all come from a uniform distribution
    if len(transcribed_regions) == 0: return []
    
    # n distinct 
    seg_graph = nx.Graph()  
    seg_graph.add_nodes_from(xrange(len(transcribed_regions)))  
    y = cov[1:] - cov[:-1]
    y[y<0] = 0
    cnt_data = [ (numpy.count_nonzero(y[start:stop+1]) + 1, 
                  start, stop)
                 for start, stop in transcribed_regions]
    #if cnt_data[0][1] > 0:
    #    cnt_data.insert(0, (1, 0, cnt_data[0][1]-1))
    #if cnt_data[-1][-1] < len(cov):
    #    cnt_data.append((1, cnt_data[-1][-1]+1, len(cov)-1))
    
    for i, (nz, start, stop) in enumerate(cnt_data[1:]):
        prev_nz, p_start, p_stop = cnt_data[i]
        old_cnt = p_stop - p_start + 1
        new_cnt = stop - p_start + 1
        merged_p = float(prev_nz + nz)/(stop - p_start + 1)
        if ( start - p_stop < 10000
             and nz < binom.ppf(1-1e-6, p=merged_p, n=new_cnt)
             and prev_nz < binom.ppf(1-0.5/len(cnt_data), 
                                     p=merged_p, n=old_cnt) ):
            seg_graph.add_edge(i, i+1)

    merged_regions = []
    for regions in nx.connected_components(seg_graph):
        merged_regions.append((cnt_data[min(regions)][1], 
                               cnt_data[max(regions)][2]))
    return merged_regions
    """

def find_transcribed_regions_and_jns_in_segment(xxx_todo_changeme, 
        rnaseq_reads, promoter_reads, polya_reads,
        ref_elements, ref_elements_to_include):
    #reg = numpy.array([1,7,6,0,0,0,4,0])
    #print reg, find_empty_regions(reg, min_length=4)
    #print reg, find_transcribed_regions(reg, min_empty_region_length=4)
    (contig, r_start, r_stop) = xxx_todo_changeme
    reg_len = r_stop-r_start+1
    cov = { '+': numpy.zeros(reg_len, dtype=float), 
            '-': numpy.zeros(reg_len, dtype=float) }
    jn_reads = {'+': defaultdict(int), '-': defaultdict(int)}
    num_unique_reads = [0.0, 0.0, 0.0]
    fragment_lengths =  defaultdict(lambda: defaultdict(int))
    
    for reads_i,reads in enumerate((promoter_reads, rnaseq_reads, polya_reads)):
        if reads == None: continue

        ( p1_rds, p2_rds, r_plus_jns, r_minus_jns, inner_cov, reads_n_uniq
          ) = extract_jns_and_reads_in_region(
              (contig, '.', r_start, r_stop), reads)
        cov['+'] += inner_cov['+']
        cov['-'] += inner_cov['-']
        num_unique_reads[reads_i] += reads_n_uniq
        
        for jn, cnt in r_plus_jns.items(): 
            jn_reads['+'][jn] += cnt 
        for jn, cnt in r_minus_jns.items(): 
            jn_reads['-'][jn] += cnt 
        

        # update the fragment length dist
        if reads == rnaseq_reads:
            for qname, r1_data in p1_rds.items():
                # if there are multiple mappings, or the read isn't paired,
                # don't use this for fragment length estimation
                if len(r1_data) != 1 or r1_data[0].map_prb != 1.0: 
                    continue
                try: r2_data = p2_rds[qname]
                except KeyError: continue
                if len(r2_data) != 1 or r2_data[0].map_prb != 1.0: 
                    continue
                # estimate the fragment length, and update the read lengths
                assert r1_data[0].map_prb == 1.0, str(r1_data)
                assert r2_data[0].map_prb == 1.0, str(r2_data)
                assert r1_data[0].read_grp == r2_data[0].read_grp
                fl_key = ( r1_data[0].read_grp, 
                           (r1_data[0].read_len, r2_data[0].read_len))
                # we add 1 because we know that this read is unique
                frag_len = calc_frag_len_from_read_data(r1_data[0], r2_data[0])
                fragment_lengths[fl_key][frag_len] += 1.0
    
    # add pseudo coverage for annotated jns. This is so that the clustering
    # algorithm knows which gene segments to join if a jn falls outside of 
    # a region with observed transcription
    # we also add pseudo coverage for other elements to provide connectivity
    if ref_elements != None and len(ref_elements_to_include) > 0:
        ref_jns = []
        for strand in "+-":
            for element_type, (start, stop) in ref_elements.iter_elements(
                    contig, strand, 
                    r_start-config.MIN_INTRON_SIZE-1, 
                    r_stop+config.MIN_INTRON_SIZE+1):
                if element_type == 'intron':
                    ref_jns.append((strand, start-r_start, stop-r_start))
                    jn_reads[strand][(start,stop)] += 0
                # add in exons
                elif element_type in ref_elements_to_include:
                    cov[strand][
                        max(0,start-r_start):max(0, stop-r_start+1)] += 1
        # add in the junction coverage
        for strand, start, stop in ref_jns:
            cov[strand][max(0, start-1-config.MIN_INTRON_SIZE):
                        max(0, start-1+1)] += 1
            cov[strand][stop:stop+config.MIN_INTRON_SIZE+1] += 1
    
    transcribed_regions = {'+': [], '-': []}
    for strand, counts in cov.items():
        transcribed_regions[strand].extend(
            sorted(find_transcribed_regions(counts)))
    
    jn_reads['+'] = sorted(jn_reads['+'].items())
    jn_reads['-'] = sorted(jn_reads['-'].items())
    
    return ( 
        transcribed_regions, jn_reads, 
        ReadCounts(*num_unique_reads), fragment_lengths )

def split_genome_into_segments(contig_lens, region_to_use, 
                               min_segment_length=5000):
    """Return non-overlapping segments that cover the genome.

    The segments are closed-closed, and strand specific.
    """
    if region_to_use != None:
        r_chrm, (r_start, r_stop) = region_to_use
    else:
        r_chrm, r_start, r_stop = None, 0, 1000000000000
    total_length = sum(contig_lens.values())
    segment_length = max(min_segment_length, 
                         int(total_length/float(config.NTHREADS*1000)))
    segments = []
    # sort by shorter contigs, so that the short contigs (e.g. mitochondrial)
    # whcih usually take longer to finish are started first
    for contig, contig_length in sorted(
            iter(contig_lens.items()), key=lambda x:x[1]):
        if region_to_use != None and r_chrm != contig: 
            continue
        for start in range(r_start, min(r_stop,contig_length), segment_length):
            segments.append(
                (contig, start, 
                 min(contig_length, start+segment_length-1)))
    return segments

class GlobalGeneSegmentData(object):
    def __init__(self, contig_lens):
        self.manager = multiprocessing.Manager()
        self.num_unique_reads = ReadCounts(multiprocessing.Value('d', 0.0), 
                                           multiprocessing.Value('d', 0.0), 
                                           multiprocessing.Value('d', 0.0))
        self.transcribed_regions = {}
        self.jns = {}
        for strand in "+-":
            for contig in list(contig_lens.keys()):
                self.transcribed_regions[(contig, strand)] = self.manager.list()
                self.jns[(contig, strand)] = self.manager.list()
        self.frag_lens = self.manager.dict()
        self.lock = multiprocessing.Lock()

    def update_all_data(self, frag_lens, transcribed_regions, jns, rd_cnts):
        with self.lock:
            for key, cnt in frag_lens.items():
                if key not in self.frag_lens: self.frag_lens[key] = cnt
                else: self.frag_lens[key] += cnt
            for key, regions in transcribed_regions.items():
                self.transcribed_regions[key].extend(regions)
            for key, vals in jns.items():
                self.jns[key].extend( vals )
            for i, val in enumerate(rd_cnts):
                self.num_unique_reads[i].value += val
        return

    def shutdown(self):
        self.manager.shutdown()
        
def find_segments_and_jns_worker(
        segments, global_gene_data,
        rnaseq_reads, promoter_reads, polya_reads,
        ref_elements, ref_elements_to_include ):
    rnaseq_reads = rnaseq_reads.reload()
    if promoter_reads != None: 
        promoter_reads = promoter_reads.reload()
    if polya_reads != None: 
        polya_reads = polya_reads.reload()

    local_frag_lens = defaultdict(int)
    local_transcribed_regions = defaultdict(list)
    local_jns = defaultdict(list)
    local_rd_cnts = [0.0, 0.0, 0.0]
    
    # just use this to keep track of where we are in the queue
    length_of_segments = segments.qsize()
    while True:
        try: 
            config.log_statement("Waiting for segment")
            segment = segments.get(timeout=1.0)
        except queue.Empty: 
            continue
        if segment == 'FINISHED': 
            config.log_statement("")
            break
        config.log_statement("Finding genes and jns in %s" % str(segment) )
        try:
            ( r_transcribed_regions, r_jns, r_n_unique_reads, r_frag_lens,
                ) = find_transcribed_regions_and_jns_in_segment(
                    segment, rnaseq_reads, promoter_reads, polya_reads, 
                    ref_elements, ref_elements_to_include) 
        except TooManyReadsError:
            seg1 = list(segment)
            seg1[2] = segment[1] + (segment[2]-segment[1])/2
            seg2 = list(segment)
            seg2[1] = seg1[2]
            segments.put(seg1)
            segments.put(seg2)
            config.log_statement("")
            continue

        for (rd_key, rls), fls in r_frag_lens.items():
            for fl, cnt in fls.items():
                local_frag_lens[(rd_key, rls, fl)] += cnt
        local_transcribed_regions[(segment[0], '+')].extend([
            (start+segment[1], stop+segment[1])
            for start, stop in r_transcribed_regions['+']])
        local_transcribed_regions[(segment[0], '-')].extend([
            (start+segment[1], stop+segment[1])
            for start, stop in r_transcribed_regions['-']])

        local_jns[(segment[0], '+')].extend(r_jns['+'])
        local_jns[(segment[0], '-')].extend(r_jns['-'])

        for i, val in enumerate(r_n_unique_reads):
            local_rd_cnts[i] += val

        if sum(local_rd_cnts) > 1e5:
            global_gene_data.update_all_data(
                local_frag_lens,
                local_transcribed_regions,
                local_jns,
                local_rd_cnts)
            
            local_frag_lens = defaultdict(int)
            local_transcribed_regions = defaultdict(list)
            local_jns = defaultdict(list)
            local_rd_cnts = [0.0, 0.0, 0.0]

    global_gene_data.update_all_data(
        local_frag_lens,
        local_transcribed_regions,
        local_jns,
        local_rd_cnts)
    
    return

def load_gene_bndry_bins( genes, contig, strand, contig_len ):  
    if config.VERBOSE:
        config.log_statement( 
            "Loading gene boundaries from annotated genes in %s:%s" % (  
                contig, strand) )  
  
    regions_graph = nx.Graph()
    for gene in genes:
        if gene.chrm != contig: continue  
        if gene.strand != strand: continue
        regions = [tuple(x) for x in gene.find_transcribed_regions()]
        regions_graph.add_nodes_from(regions)
        regions_graph.add_edges_from(zip(regions[:-1], regions[1:]))

    # group overlapping regions
    all_regions = sorted(regions_graph.nodes())
    if len(all_regions) == 0: return []  

    grpd_regions = [[],]
    curr_start, curr_stop = all_regions[0]
    for x in all_regions:
        if x[0] < curr_stop:
            curr_stop = max(x[1], curr_stop)
            grpd_regions[-1].append(x)
        else:
            curr_start, curr_stop = x
            grpd_regions.append([x,])
    # add edges for overlapping regions
    for grp in grpd_regions:
        regions_graph.add_edges_from(zip(grp[:-1], grp[1:]))

    # build gene objects with the intervals  
    gene_bndry_bins = []  
    for regions_cluster in nx.connected_components(regions_graph):
        gene_bin = GeneElements( contig, strand )
        regions = sorted(files.gtf.flatten(regions_cluster))
        for start, stop in regions:
            gene_bin.regions.append(
                SegmentBin(start, stop, ["ESTART",], ["ESTOP",], "GENE"))
        gene_bndry_bins.append( gene_bin )  

    # XXX TODO expand gene boundaries
    # actually, it's probably better just to go through discovery
    
    return gene_bndry_bins

def find_all_gene_segments( rnaseq_reads, promoter_reads, polya_reads,
                            ref_genes, ref_elements_to_include,
                            region_to_use=None ):
    config.log_statement("Finding gene segments")

    contig_lens = dict(list(zip(*get_contigs_and_lens( 
        [ reads for reads in [rnaseq_reads, promoter_reads, polya_reads]
          if reads != None ] ))))

    config.log_statement("Spawning gene segment finding children")    
    segments_queue = multiprocessing.Queue()
    global_gene_data = GlobalGeneSegmentData(contig_lens)
    
    ref_element_types_to_include = set()
    if ref_elements_to_include.junctions: 
        ref_element_types_to_include.add('intron')
    if ref_elements_to_include.TSS: 
        ref_element_types_to_include.add('tss_exon')
    if ref_elements_to_include.TES: 
        ref_element_types_to_include.add('tes_exon')
    if ref_elements_to_include.promoters: 
        ref_element_types_to_include.add('promoter')
    if ref_elements_to_include.polya_sites: 
        ref_element_types_to_include.add('polya')
    if ref_elements_to_include.exons: 
        ref_element_types_to_include.add('exon')
    # to give full gene connectivity
    if ref_elements_to_include.genes:
        ref_element_types_to_include.add('intron')
        ref_element_types_to_include.add('exon')
    
    pids = []
    for i in range(config.NTHREADS):
        pid = os.fork()
        if pid == 0:
            find_segments_and_jns_worker(
                segments_queue, 
                global_gene_data,
                rnaseq_reads, promoter_reads, polya_reads,
                ref_genes, ref_element_types_to_include)
            os._exit(0)
        pids.append(pid)

    config.log_statement("Populating gene segment queue")        
    segments = split_genome_into_segments(contig_lens, region_to_use)
    for segment in segments: 
        segments_queue.put(segment)
    for i in range(config.NTHREADS): segments_queue.put('FINISHED')
    
    while segments_queue.qsize() > 2*config.NTHREADS:
        config.log_statement(
            "Waiting on gene segment finding children (%i/%i segments remain)" 
            %(segments_queue.qsize(), len(segments)))        
        time.sleep(0.5)
    
    for i, pid in enumerate(pids):
        config.log_statement(
            "Waiting on gene segment finding children (%i/%i children remain)" 
            %(len(pids)-i, len(pids)))
        os.waitpid(pid, 0) 
            
    config.log_statement("Merging gene segments")
    merged_transcribed_regions = {}
    for key, intervals in global_gene_data.transcribed_regions.items():
        merged_transcribed_regions[
            key] = merge_adjacent_intervals(
                intervals, config.MAX_EMPTY_REGION_SIZE)
    transcribed_regions = merged_transcribed_regions
    
    config.log_statement("Filtering junctions")    
    filtered_jns = defaultdict(dict)
    for contig in list(contig_lens.keys()):
        plus_jns = defaultdict(int)
        for jn, cnt in global_gene_data.jns[(contig, '+')]: plus_jns[jn] += cnt
        minus_jns = defaultdict(int)
        for jn, cnt in global_gene_data.jns[(contig, '-')]: minus_jns[jn] += cnt
        filtered_jns[(contig, '+')] = filter_jns(plus_jns, minus_jns)
        filtered_jns[(contig, '-')] = filter_jns(minus_jns, plus_jns)

    config.log_statement("Building FL dist")        
    fl_dists = build_fl_dists_from_fls_dict(dict(global_gene_data.frag_lens))
        
    if ref_elements_to_include.junctions:
        for gene in ref_genes:
            for jn in gene.extract_elements()['intron']:
                if jn not in filtered_jns[(gene.chrm, gene.strand)]:
                    filtered_jns[(gene.chrm, gene.strand)][jn] = 0
        
    config.log_statement("Clustering gene segments")    
    # build bins for all of the genes and junctions, converting them to 1-based
    # in the process
    new_genes = []
    new_introns = []
    for contig, contig_len in contig_lens.items():
        for strand in '+-':
            key = (contig, strand)
            jns = [ (start, stop, cnt) 
                    for (start, stop), cnt 
                    in sorted(filtered_jns[key].items()) ]
            for start, stop, cnt in jns:
                new_introns.append(
                    SegmentBin(start, stop, ["D_JN",], ["R_JN",], "INTRON"))
            intervals = cluster_intron_connected_segments(
                transcribed_regions[key], 
                [(start, stop) for start, stop, cnt in jns ] )
            # add the intergenic space, since there could be interior genes
            for segments in intervals: 
                new_gene = GeneElements( contig, strand )
                for start, stop in segments:
                    new_gene.regions.append( 
                        SegmentBin(start, stop, ["ESTART",],["ESTOP",],"GENE"))
                if new_gene.stop-new_gene.start+1 < config.MIN_GENE_LENGTH: 
                    continue
                new_genes.append(new_gene)

    try: 
        num_unique_reads = ReadCounts(*[
            float(x.value) for x in global_gene_data.num_unique_reads])
    except AttributeError:
        num_unique_reads = ReadCounts(*global_gene_data.num_unique_reads)

    global_gene_data.shutdown()

    config.log_statement("")    
        
    return new_genes, fl_dists, num_unique_reads 
