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

import sys, os
import time
import math
import traceback

import shutil

import numpy
from scipy.stats import beta, binom

from collections import defaultdict, namedtuple
from itertools import chain, izip
from bisect import bisect
from copy import copy

import multiprocessing
import Queue

import networkx as nx

from genes import (
    GeneElements, SegmentBin, TranscriptBoundaryBin, 
    find_all_gene_segments, TranscriptElement
)

from files.reads import MergedReads, RNAseqReads, CAGEReads, \
    RAMPAGEReads, PolyAReads, \
    fix_chrm_name_for_ucsc, get_contigs_and_lens, calc_frag_len_from_read_data, \
    iter_paired_reads, extract_jns_and_reads_in_region
import files.junctions
from files.bed import create_bed_line
from files.gtf import parse_gtf_line, load_gtf

from peaks import call_peaks, build_control_in_gene

from elements import find_jn_connected_exons

from frag_len import FlDist, find_fls_from_annotation

from transcript import Transcript, Gene
import f_matrix     
import frequency_estimation
frequency_estimation.LHD_ABS_TOL = 1e-1
frequency_estimation.PARAM_ABS_TOL = 1e-3

import config

class ThreadSafeFile( file ):
    def __init__( self, *args ):
        args = list( args )
        args.insert( 0, self )
        file.__init__( *args )
        self.lock = multiprocessing.Lock()

    def write( self, line ):
        with self.lock:
            file.write( self, line )
            self.flush()


def filter_exon(exon, wig, num_start_bases_to_skip=0, num_stop_bases_to_skip=0):
    '''Find all the exons that are sufficiently homogenous and expressed.
    
    '''
    start = exon.start + num_start_bases_to_skip
    end = exon.stop - num_stop_bases_to_skip
    if start >= end - 10: return False
    vals = wig[start:end+1]
    n_div = max( 1, int(len(vals)/config.MAX_EMPTY_REGION_SIZE) )
    div_len = len(vals)/n_div
    for i in xrange(n_div):
        seg = vals[i*div_len:(i+1)*div_len]
        if seg.mean() < config.MIN_EXON_AVG_CVG:
            return True

    return False

def filter_exons( exons, rnaseq_cov, 
                  num_start_bases_to_skip=0, 
                  num_stop_bases_to_skip=0 ):
    for exon in exons:
        if not filter_exon( exon, rnaseq_cov, 
                            num_start_bases_to_skip, 
                            num_stop_bases_to_skip ):
            yield exon
    
    return        


def reverse_strand(bins_iter, contig_len):
    rev_bins = []
    for bin in reversed(bins_iter):
        rev_bins.append( bin.reverse_strand( contig_len ) )
    return rev_bins

class SpliceGraph(nx.DiGraph):
    pass

def find_cage_peak_bins_in_gene( gene, cage_reads, rnaseq_reads ):
    rnaseq_cov = gene.find_coverage( rnaseq_reads )
    print rnaseq_cov
    cage_cov = gene.find_coverage( cage_reads)
    print cage_cov
    assert False
    # threshold the CAGE data. We assume that the CAGE data is a mixture of 
    # reads taken from actually capped transcripts, and random transcribed 
    # regions, or RNA seq covered regions. We zero out any bases where we
    # can't reject the null hypothesis that the observed CAGE reads all derive 
    # from the background, at alpha = 0.001. 
    rnaseq_cov = numpy.array( rnaseq_cov+1-1e-6, dtype=int)
    max_val = rnaseq_cov.max()
    thresholds = config.TOTAL_MAPPED_READS*beta.ppf( 
        config.CAGE_FILTER_ALPHA, 
        numpy.arange(max_val+1)+1, 
        numpy.zeros(max_val+1)+(config.TOTAL_MAPPED_READS+1) 
    )
    max_scores = thresholds[ rnaseq_cov ]
    cage_cov[ cage_cov < max_scores ] = 0    
    
    raw_peaks = find_peaks( cage_cov, window_len=config.CAGE_PEAK_WIN_SIZE, 
                            min_score=config.MIN_NUM_CAGE_TAGS,
                            max_score_frac=config.MAX_CAGE_FRAC,
                            max_num_peaks=100)
    
    cage_peak_bins = [] #Bins( gene.chrm, gene.strand )
    for pk_start, pk_stop in raw_peaks:
        cnt = float(cage_cov[pk_start:pk_stop+1].sum())
        bin = SegmentBin(
            gene.start+pk_start, gene.start+pk_stop, 
            "CAGE_PEAK_START", "CAGE_PEAK_STOP", "CAGE_PEAK",
            1e6*beta.ppf(0.01, cnt+1e-6, cage_reads.num_reads+1e-6),
            1e6*beta.ppf(0.50, cnt+1e-6, cage_reads.num_reads+1e-6),
            1e6*beta.ppf(0.99, cnt+1e-6, cage_reads.num_reads+1e-6) )
        cage_peak_bins.append(bin)
    
    return cage_peak_bins

def find_polya_peak_bins_in_gene( gene, polya_reads, rnaseq_reads ):
    polya_cov = gene.find_coverage(polya_reads)
    
    # threshold the polya data. We assume that the polya data is a mixture of 
    # reads taken from actually capped transcripts, and random transcribed 
    # regions, or RNA seq covered regions. We zero out any bases where we
    # can't reject the null hypothesis that the observed polya reads all derive 
    # from the background, at alpha = 0.001. 
    """
    rnaseq_cov = gene.find_coverage( rnaseq_reads )
    rnaseq_cov = numpy.array( rnaseq_cov+1-1e-6, dtype=int)
    max_val = rnaseq_cov.max()
    thresholds = TOTAL_MAPPED_READS*beta.ppf( 
        0.1, 
        numpy.arange(max_val+1)+1, 
        numpy.zeros(max_val+1)+(TOTAL_MAPPED_READS+1) 
    )
    max_scores = thresholds[ rnaseq_cov ]
    polya_cov[ polya_cov < max_scores ] = 0    
    """
    
    raw_peaks = find_peaks( polya_cov, window_len=30, 
                            min_score=config.MIN_NUM_POLYA_TAGS,
                            max_score_frac=0.05,
                            max_num_peaks=100)
    polya_sites = [] #Bins( gene.chrm, gene.strand )
    if len( raw_peaks ) == 0:
        return polya_sites
    
    for pk_start, pk_stop in raw_peaks:
        cnt = float(polya_cov[pk_start:pk_stop+1].sum())
        bin = SegmentBin(
            gene.start+pk_start, gene.start+pk_stop, 
            "POLYA_PEAK_START", "POLYA_PEAK_STOP", "POLYA",
            1e6*beta.ppf(0.01, cnt+1e-6, polya_reads.num_reads+1e-6),
            1e6*beta.ppf(0.50, cnt+1e-6, polya_reads.num_reads+1e-6),
            1e6*beta.ppf(0.99, cnt+1e-6, polya_reads.num_reads+1e-6) )
        polya_sites.append(bin)
    
    return polya_sites

def find_peaks( cov, window_len, min_score, max_score_frac, max_num_peaks ):    
    def overlaps_prev_peak( new_loc ):
        for start, stop in peaks:
            if not( new_loc > stop or new_loc + window_len < start ):
                return True
        return False
    
    # merge the peaks
    def grow_peak( start, stop, grow_size=
                   max(1, window_len/4), min_grow_ratio=config.MAX_CAGE_FRAC ):
        # grow a peak at most max_num_peaks times
        max_mean_signal = cov[start:stop+1].mean()
        for i in xrange(max_num_peaks):
            curr_signal = cov[start:stop+1].sum()
            if curr_signal < min_score:
                return ( start, stop )
            
            downstream_sig = float(cov[max(0, start-grow_size):start].sum())/grow_size
            upstream_sig = float(cov[stop+1:stop+1+grow_size].sum())/grow_size
            
            # if neither passes the threshold, then return the current peak
            if max(upstream_sig, downstream_sig) \
                    < min_grow_ratio*curr_signal/float(stop-start+1): 
                return (start, stop)
            
            # if the expansion isn't greater than the min ratio, then return
            if max(upstream_sig,downstream_sig) < \
                    config.MAX_CAGE_FRAC*max_mean_signal:
                return (start, stop)
            
            # otherwise, we know one does
            if upstream_sig > downstream_sig:
                stop += grow_size
            else:
                start = max(0, start - grow_size )
        
        if config.VERBOSE:
            config.log_statement( 
                "Warning: reached max peak iteration at %i-%i ( signal %.2f )"
                    % (start, stop, cov[start:stop+1].sum() ) )
        return (start, stop )
    
    peaks = []
    peak_scores = []
    cumsum_cvg_array = (
        numpy.append(0, numpy.cumsum( cov )) )
    scores = cumsum_cvg_array[window_len:] - cumsum_cvg_array[:-window_len]
    indices = numpy.argsort( scores )
    min_score = max( min_score, config.MAX_CAGE_FRAC*scores[ indices[-1] ] )
    for index in reversed(indices):
        if not overlaps_prev_peak( index ):
            score = scores[ index ]
            new_peak = grow_peak( index, index + window_len )
            # if we are below the minimum score, then we are done
            if score < min_score:
                break
            
            # if we have observed peaks, and the ratio between the highest
            # and the lowest is sufficeintly high, we are done
            if len( peak_scores ) > 0:
                if float(score)/peak_scores[0] < max_score_frac:
                    break
                        
            peaks.append( new_peak ) 
            peak_scores.append( score )
    
    if len( peaks ) == 0:
        return []
    
    # merge cage peaks together
    def merge_peaks( peaks_and_scores ):
        peaks_and_scores = sorted( list(x) for x in peaks_and_scores )
        peak, score = peaks_and_scores.pop()
        new_peaks = [peak,]
        new_scores = [score,]
        while len(peaks_and_scores) >  0:
            last_peak = new_peaks[-1]
            peak, score = peaks_and_scores.pop()
            new_peak = (min(peak[0], last_peak[0]), max(peak[1], last_peak[1]))
            if (new_peak[1] - new_peak[0]) <= 1.5*( 
                    last_peak[1] - last_peak[0] + peak[1] - peak[0] ):
                new_peaks[-1] = new_peak
                new_scores[-1] += score
            else:
                new_peaks.append( peak )
                new_scores.append( score )
        
        return zip( new_peaks, new_scores )
    
    peaks_and_scores = sorted( zip(peaks, peak_scores) )
    
    for i in xrange( 99 ):
        if i == 100: assert False
        old_len = len( peaks_and_scores )
        peaks_and_scores = merge_peaks( peaks_and_scores )
        if len( peaks_and_scores ) == old_len: break
    
        
    new_peaks_and_scores = []
    for peak, score in peaks_and_scores:
        peak_scores = cov[peak[0]:peak[1]+1]
        max_score = peak_scores.max()
        good_indices = (peak_scores >= max_score*config.MAX_CAGE_FRAC).nonzero()[0]
        new_peak = [
                peak[0] + int(good_indices.min()), 
                peak[0] + int(good_indices.max())  ]
        new_score = float(cov[new_peak[0]:new_peak[1]+1].sum())
        new_peaks_and_scores.append( (new_peak, new_score) )
    
    peaks_and_scores = sorted( new_peaks_and_scores )
    max_score = max( s for p, s in peaks_and_scores )
    return [ pk for pk, score in peaks_and_scores \
                 if score >= config.MAX_CAGE_FRAC*max_score
                 and score > min_score ]

"""
def estimate_peak_expression(peaks, peak_cov, distal_reads):
    for peak in peaks:
        cnt = float(peak_cov[peak.start:peak.stop+1].sum())
        peak.tpm = cnt/distal_reads.num_reads
        peak.tpm_ub = 1e6*beta.ppf(0.99, cnt+1, distal_reads.num_reads)/distal_reads.num_reads
        peak.tpm_lb = 1e6*beta.ppf(0.01, cnt+1, distal_reads.num_reads)/distal_reads.num_reads
    
    return peaks
"""

def iter_retained_intron_connected_exons(
        left_exons, right_exons, retained_introns,
        max_num_exons=None):
    graph = nx.DiGraph()
    left_exons = sorted((x.start, x.stop) for x in left_exons)
    right_exons = sorted((x.start, x.stop) for x in right_exons)

    graph.add_nodes_from(chain(left_exons, right_exons))
    retained_jns = [ (x.start+1, x.stop-1) for x in retained_introns ]
    edges = find_jn_connected_exons(set(chain(left_exons, right_exons)), 
                                    retained_jns, '+')
    graph.add_edges_from( (start, stop) for jn, start, stop in edges )    
    cntr = 0
    for left in left_exons:
        for right in right_exons:
            if left[1] > right[0]: continue
            for x in nx.all_simple_paths(
                    graph, left, right, max_num_exons-cntr+1):
                cntr += 1
                if max_num_exons != None and cntr > max_num_exons:
                    raise ValueError, "Too many retained introns"
                yield (x[0][0], x[-1][1])
    return

def merge_tss_exons(tss_exons):
    grpd_exons = defaultdict(list)
    for exon in tss_exons:
        grpd_exons[exon.stop].append(exon)
    merged_tss_exons = []
    for stop, exons in grpd_exons.iteritems():
        exons.sort(key=lambda x:x.start)
        curr_start = exons[0].start
        score = exons[0].score
        for exon in exons[1:]:
            if exon.start - curr_start < config.TSS_EXON_MERGE_DISTANCE:
                score += exon.score
            else:
                merged_tss_exons.append( Bin(curr_start, stop,
                                             exon.left_label, 
                                             exon.right_label,
                                             exon.type, score) )
                curr_start = exon.start
                score = exon.score
            
        merged_tss_exons.append( Bin(curr_start, stop,
                                     exon.left_label, 
                                     exon.right_label,
                                     exon.type, score) )
    return merged_tss_exons

def merge_tes_exons(tes_exons):
    grpd_exons = defaultdict(list)
    for exon in tes_exons:
        grpd_exons[exon.start].append(exon)
    merged_tes_exons = []
    for start, exons in grpd_exons.iteritems():
        exons.sort(key=lambda x:x.stop)
        curr_stop = exons[0].stop
        score = exons[0].score
        for exon in exons[1:]:
            if exon.stop - curr_stop < config.TES_EXON_MERGE_DISTANCE:
                curr_stop = exon.stop
                score += exon.score
            else:
                merged_tes_exons.append( Bin(start, curr_stop,
                                            exon.left_label, 
                                            exon.right_label,
                                            exon.type, score) )
                curr_stop = exon.stop
                score = exon.score

        merged_tes_exons.append( Bin(start, curr_stop,
                                     exon.left_label, 
                                     exon.right_label,
                                     exon.type, score) )
    return merged_tes_exons


def find_transcribed_fragments_covering_region(
        segment_graph, segment_id,
        max_frag_len, use_genome_coords=False):
    if use_genome_coords:
        raise NotImplemented, "This code path probably doesn't work anymore, it hasn't been updated since this method was changed to use the splice graph directly"
    
    def build_neighbor_paths(side):
        assert side in ('BEFORE', 'AFTER')
        complete_paths = []
        partial_paths = [([segment_id,], 0),]
        while len(partial_paths) > 0:
            curr_path, curr_path_len = partial_paths.pop()
            if side == 'BEFORE':
                neighbors = list(
                    segment_graph.predecessors(curr_path[0]))
            else:
                neighbors = list(
                    segment_graph.successors(curr_path[-1]))
            if len(neighbors) == 0: 
                complete_paths.append((curr_path, curr_path_len))
            else:
                for child in neighbors:
                    if segment_graph.node[child]['type'] in ('TSS', 'TES'):
                        complete_paths.append((curr_path, curr_path_len))
                        continue
                    assert segment_graph.node[child]['type'] == 'segment'
                    
                    if side == 'BEFORE':
                        new_path = [child,] + curr_path
                    else:
                        new_path = curr_path + [child,]
                    
                    new_path_len = ( 
                        segment_graph.node[child]['bin'].length()+curr_path_len)
                    if new_path_len >= max_frag_len:
                        complete_paths.append((new_path, new_path_len))
                    else:
                        partial_paths.append((new_path, new_path_len))
        
        if side == 'BEFORE': return [x[:-1] for x, x_len in complete_paths]
        else: return [x[1:] for x, x_len in complete_paths]
    
    def build_genome_segments_from_path(segment_indexes):        
        merged_intervals = merge_adjacent_intervals(
            zip(segment_indexes, segment_indexes), 
            max_merge_distance=0)
        coords = []
        for start, stop in merged_intervals:
            coords.append([genes_graph.node[start]['bin'].start,
                           genes_graph.node[stop+1]['bin'].stop,])
        return coords
        
    complete_before_paths = build_neighbor_paths('BEFORE')
    complete_after_paths = build_neighbor_paths('AFTER')    
    transcripts = []
    for bp in complete_before_paths:
        for ap in complete_after_paths:
            segments = bp + [segment_id,] + ap
            #assert segments == sorted(segments), str(segments)
            if use_genome_coords:
                before_trim_len = max(
                    0, sum(genes_graph.node[i]['node'].length() 
                           for i in bp) - max_frag_len)
                after_trim_len = max(
                    0, sum(genes_graph.node[i]['node'].length() 
                           for i in ap) - max_frag_len)
                segments = build_genome_segments_from_path(segments)
                segments[0][0] = segments[0][0] + before_trim_len
                segments[-1][1] = segments[-1][1] - after_trim_len
            transcripts.append( sorted(segments) )
    return transcripts

def extract_jns_and_paired_reads_in_gene(gene, reads):
    pair1_reads = defaultdict(list)
    pair2_reads = defaultdict(list)
    plus_jns = defaultdict(int)
    minus_jns = defaultdict(int)
    
    for region in gene.regions:
        ( r_pair1_reads, r_pair2_reads, r_plus_jns, r_minus_jns, 
          ) = extract_jns_and_reads_in_region(
            (gene.chrm, gene.strand, region.start, region.stop), reads)
        for jn, cnt in r_plus_jns.iteritems(): 
            plus_jns[jn] += cnt 
        for jn, cnt in r_minus_jns.iteritems(): 
            minus_jns[jn] += cnt 
        for qname, read_mappings in r_pair1_reads.iteritems(): 
            pair1_reads[qname].extend(read_mappings)
        for qname, read_mappings in r_pair2_reads.iteritems(): 
            pair2_reads[qname].extend(read_mappings)

    paired_reads = list(iter_paired_reads(pair1_reads, pair2_reads))
    jns, opp_strand_jns = (
        (plus_jns, minus_jns) if gene.strand == '+' else (minus_jns, plus_jns)) 
    return paired_reads, jns, opp_strand_jns

def find_widest_path(splice_graph):
    return None, 0
    try: 
        tss_max = max(data['bin'].fpkm_lb
                      for node, data in splice_graph.nodes(data=True)
                      if data['type'] == 'TSS')
    except:
        tss_max = 0
    try:
        tes_max = max(data['bin'].fpkm_lb
                      for node, data in splice_graph.nodes(data=True)
                      if data['type'] == 'TES')
    except:
        tes_max = 0
    return None, max(tss_max, tss_max)

    splice_graph = copy(splice_graph)
    assert nx.is_directed_acyclic_graph(splice_graph)
    config.log_statement("Finding widest path")
    curr_paths = [ [[node,], data['bin'].fpkm, data['bin'].fpkm] 
                  for node, data in splice_graph.nodes(data=True)
                  if data['type'] == 'TSS']
    max_path = None
    max_min_fpkm = min( x[1] for x in curr_paths )/10. if len(curr_paths) > 0 else 0
    while len(curr_paths) > 0:
        # find a path to check
        curr_path, curr_fpkm, prev_fpkm = curr_paths.pop()
        # if the paths fpkm is less than the maximum,
        # then there is no need to consider anything else
        if curr_fpkm < max_min_fpkm: continue
        successors = splice_graph.successors(curr_path[-1])

        # if there are no successors, then still allow this path
        # to contribute in case we cant build a full path
        # (e.g. there are no TES sites)
        if len(successors) == 0:
            if max_min_fpkm < curr_fpkm:
                max_path = curr_path
                max_min_fpkm = curr_fpkm

        config.log_statement("%i paths in queue (%e, %e, %i, %i)" % (
            len(curr_paths), max_min_fpkm, prev_fpkm, 
            len(successors), len([x for x, data in splice_graph.nodes(data=True)
                                  if data['bin'].fpkm > max_min_fpkm])))
            
        # otherwise, build new paths and update the fpkms            
        for successor in successors:
            bin = splice_graph.node[successor]['bin']
            new_path = curr_path + [successor,]
            new_fpkm = min(bin.fpkm, curr_fpkm)
            if new_fpkm < max_min_fpkm: 
                splice_graph.remove_node(successor)
                continue
            if bin.type == 'POLYA':
                if new_fpkm > max_min_fpkm:
                    max_path = new_path
                    max_min_fpkm = new_fpkm
            else:
                curr_paths.append((new_path, new_fpkm, bin.fpkm))
        curr_paths.sort(key=lambda x:len(x[0]))
    return max_path, max_min_fpkm

def build_splice_graph_and_binned_reads_in_gene( 
        gene, 
        rnaseq_reads, tss_reads, tes_reads, ref_elements ):
    assert isinstance( gene, GeneElements )
    config.log_statement( 
        "Extracting reads and jns in Chrm %s Strand %s Pos %i-%i" %
        (gene.chrm, gene.strand, gene.start, gene.stop) )
    
    # initialize the cage peaks with the reference provided set
    tss_regions = [ Bin(pk_start, pk_stop, 
                     "CAGE_PEAK_START", "CAGE_PEAK_STOP", "CAGE_PEAK")
                 for pk_start, pk_stop in ref_elements['promoter'] ]
    
    # initialize the polya peaks with the reference provided set
    tes_regions = [ Bin( pk_start, pk_stop, 
                      "POLYA_PEAK_START", "POLYA_PEAK_STOP", "POLYA")
                 for pk_start, pk_stop in ref_elements['polya'] ]
    
    # build and pair rnaseq reads, and extract junctions
    (paired_rnaseq_reads, observed_jns, opp_strand_jns 
     ) = extract_jns_and_paired_reads_in_gene(gene, rnaseq_reads)
        
    jns = files.junctions.filter_jns(
        observed_jns, opp_strand_jns, set(ref_elements['introns']))
    # add in connectivity junctions
    for distal_reads in (tss_reads, tes_reads):
        if distal_reads == None: continue
        for jn, cnt, entropy in files.junctions.load_junctions_in_bam(
              distal_reads, 
              [ (gene.chrm, gene.strand, r.start, r.stop) for r in gene.regions]
              )[(gene.chrm, gene.strand)]:
            jns[jn] += 0
    
    # add in reference junctions
    for jn in ref_elements['introns']: jns[jn] += observed_jns[jn]
        
    config.log_statement( 
        "Building exon segments in Chrm %s Strand %s Pos %i-%i" %
        (gene.chrm, gene.strand, gene.start, gene.stop) )
    
    # build the pseudo exon set
    segment_bnds = set([gene.start, gene.stop+1])
    segment_bnd_labels = defaultdict(set)
    segment_bnd_labels[gene.start].add("GENE_BNDRY")
    segment_bnd_labels[gene.stop+1].add("GENE_BNDRY")
    # add in empty regions - we will use these to filter bins
    # that fall outside of the gene
    for r1, r2 in izip(gene.regions[:-1], gene.regions[1:]):
        assert r1.stop+2 < r2.start-1+2
        segment_bnds.add(r1.stop+1)
        segment_bnd_labels[r1.stop+1].add( 'EMPTY_START' )
        segment_bnds.add(r2.start-1+1)
        segment_bnd_labels[r2.start-1+1].add( 'EMPTY_STOP' )
    
    for (start, stop) in jns:
        segment_bnds.add(start)
        segment_bnd_labels[start].add('D_JN' if gene.strand == '+' else 'R_JN')
        segment_bnds.add(stop+1)
        segment_bnd_labels[stop+1].add('R_JN' if gene.strand == '+' else 'D_JN')
    
    if tss_reads != None:
        control_cov = build_control_in_gene(
            gene, paired_rnaseq_reads, sorted(segment_bnds), 
            '5p' if gene.strand == '+' else '3p')

        signal_cov = gene.find_coverage( tss_reads )
        for pk_start, pk_stop, peak_cov in call_peaks( 
                signal_cov, control_cov,
                '5p' if gene.strand == '+' else '3p',
                gene,
                **config.TSS_call_peaks_tuning_params):
            tss_regions.append(TranscriptBoundaryBin(
                gene.start+pk_start, gene.start+pk_stop-1, 
                "CAGE_PEAK_START", "CAGE_PEAK_STOP", "CAGE_PEAK",
                peak_cov
                ).set_tpm(tss_reads.num_reads))
    
    if tes_reads != None:
        control_cov = build_control_in_gene(
            gene, paired_rnaseq_reads, sorted(segment_bnds), 
            '3p' if gene.strand == '+' else '5p')
        signal_cov = gene.find_coverage( tes_reads )
        for pk_start, pk_stop, peak_cov in call_peaks( 
                signal_cov, control_cov,
                '3p' if gene.strand == '+' else '5p',
                gene,
                **config.TES_call_peaks_tuning_params):
            tes_regions.append(TranscriptBoundaryBin(
                gene.start+pk_start, gene.start+pk_stop-1, 
                "POLYA_PEAK_START", "POLYA_PEAK_STOP", "POLYA",
                peak_cov,
                ).set_tpm(tes_reads.num_reads))
    
    splice_graph = SpliceGraph()
    tss_segment_map = {}
    for tss_bin in tss_regions:
        tss_start = tss_bin.start if gene.strand == '+' else tss_bin.stop + 1
        segment_bnds.add(tss_start)
        segment_bnd_labels[tss_start].add('TSS')
        tss_segment_map[tss_bin] = [tss_start,]

    tes_segment_map = defaultdict(list)    
    for tes_bin in tes_regions:
        tes_start = tes_bin.stop+1 if gene.strand == '+' else tes_bin.start
        segment_bnds.add(tes_start)
        segment_bnd_labels[tes_start].add('TES')
        tes_segment_map[tes_bin] = [tes_start,]
    
    # remove boundaries that fall inside of empty regions
    empty_segments = set()
    in_empty_region = False
    for index, bnd in enumerate(sorted(segment_bnds)):
        if in_empty_region:
            assert 'EMPTY_STOP' in segment_bnd_labels[bnd], \
                str((gene, bnd, sorted(segment_bnd_labels.iteritems()), \
                     sorted(jns.iteritems())))
        empty_start = bool('EMPTY_START' in segment_bnd_labels[bnd])
        empty_stop = bool('EMPTY_STOP' in segment_bnd_labels[bnd])
        assert not (empty_start and empty_stop)
        if empty_start: in_empty_region = True
        if empty_stop: in_empty_region = False
        if in_empty_region:
            empty_segments.add(index)

    # build the exon segment connectivity graph
    segment_bnds = numpy.array(sorted(segment_bnds))
    
    for element_i in (x for x in xrange(0, len(segment_bnds)-1) 
                      if x not in empty_segments):
        start, stop = segment_bnds[element_i], segment_bnds[element_i+1]-1
        left_labels = segment_bnd_labels[start]
        right_labels = segment_bnd_labels[stop+1]
        bin = SegmentBin( start, stop, left_labels, right_labels, type=None)
        splice_graph.add_node(element_i, type='segment', bin=bin)
    
    for i in xrange(len(segment_bnds)-1-1):
        if i not in empty_segments and i+1 not in empty_segments:
            if gene.strand == '+':
                splice_graph.add_edge(i, i+1, type='adjacent', bin=None)
            else:
                splice_graph.add_edge(i+1, i, type='adjacent', bin=None)

    jn_edges = []
    for (start, stop), cnt in jns.iteritems():
        start_i = segment_bnds.searchsorted(start)-1
        assert segment_bnds[start_i+1] == start
        stop_i = segment_bnds.searchsorted(stop+1)-1+1
        assert segment_bnds[stop_i] == stop+1
        assert start_i in splice_graph
        # skip junctions that splice to the last base in the gene XXX
        if stop_i == splice_graph.number_of_nodes(): continue
        assert stop_i in splice_graph, str((stop_i, splice_graph.nodes(data=True)))
        if gene.strand == '+':
            bin = SegmentBin( start, stop, 'D_JN', 'R_JN', type='INTRON', cnt=cnt)
            splice_graph.add_edge(start_i, stop_i, type='splice', bin=bin)
        else:
            bin = SegmentBin( start, stop, 'R_JN', 'D_JN', type='INTRON', cnt=cnt)
            splice_graph.add_edge(stop_i, start_i, type='splice', bin=bin)

    for i, (tss_bin, tss_segments) in enumerate(tss_segment_map.iteritems()):
        node_id = "TSS_%i" % i
        splice_graph.add_node( node_id, type='TSS', bin=tss_bin)
        for tss_start in tss_segments:
            bin_i = segment_bnds.searchsorted(tss_start)
            assert segment_bnds[bin_i] == tss_start
            if gene.strand == '-': bin_i -= 1
            splice_graph.add_edge(node_id, bin_i, type='tss', bin=None)
    
    for i, (tes_bin, tes_segments) in enumerate(tes_segment_map.iteritems()):
        node_id = "TES_%i"% i
        splice_graph.add_node( node_id, type='TES', bin=tes_bin)
        for tes_start in tes_segments:
            bin_i = segment_bnds.searchsorted(tes_start)
            assert segment_bnds[bin_i] == tes_start
            if gene.strand == '+': bin_i -= 1
            splice_graph.add_edge(bin_i, node_id, type='tes', bin=None)
    
    return splice_graph, None
    
    config.log_statement( 
        "Binning reads in Chrm %s Strand %s Pos %i-%i" %
        (gene.chrm, gene.strand, gene.start, gene.stop) )
    
    # pre-calculate the expected and observed counts
    binned_reads = f_matrix.bin_rnaseq_reads( 
        rnaseq_reads, gene.chrm, gene.strand, segment_bnds, 
        include_read_type=False)
    
    return splice_graph, binned_reads

def fast_quantify_segment_expression(gene, splice_graph, 
                                     rnaseq_reads, cage_reads, polya_reads):
    config.log_statement( 
        "FAST Quantifying segment expression in Chrm %s Strand %s Pos %i-%i" %
        (gene.chrm, gene.strand, gene.start, gene.stop) )

    quantiles = [0.01, 0.5, 1-.01]
    avg_read_len = 0
    for ((rg, (r1_len,r2_len)), (fl_dist, marginal_frac) 
             ) in rnaseq_reads.fl_dists.iteritems():
        avg_read_len += marginal_frac*(r1_len + r2_len)/2
    
    rnaseq_cov = gene.find_coverage(rnaseq_reads)
    for element_i, data in splice_graph.nodes(data=True):
        #element = splice_graph.nodes(element_i)
        element = splice_graph.node[element_i]
        if element['type'] != 'segment': continue
        coverage = rnaseq_cov[element['bin'].start-gene.start
                              :element['bin'].stop-gene.start+1]
        n_reads_in_segment = 0 if len(coverage) == 0 else numpy.median(coverage)
        fpkms = 1e6*(1000./element['bin'].length())*beta.ppf(
            quantiles, 
            n_reads_in_segment+1e-6, 
            rnaseq_reads.num_reads-n_reads_in_segment+1e-6)
        element['bin'].set_expression(*fpkms)

    for start, stop, data in splice_graph.edges(data=True):
        if data['type'] != 'splice': continue
        n_reads_in_segment = float(data['bin'].cnt)
        effective_len = float(max(
            1, avg_read_len - 2*config.MIN_INTRON_FLANKING_SIZE))
        fpkms = 1e6*(1000./effective_len)*beta.ppf(
            quantiles, 
            n_reads_in_segment+1e-6, 
            rnaseq_reads.num_reads-n_reads_in_segment+1e-6)
        data['bin'].set_expression(*fpkms)

    return splice_graph

def quantify_segment_expression(gene, splice_graph, binned_reads ):
    config.log_statement( 
        "Quantifying segment expression in Chrm %s Strand %s Pos %i-%i" %
        (gene.chrm, gene.strand, gene.start, gene.stop) )

    # find all possible transcripts, and their association with each element
    max_fl = max(fl_dist.fl_max for (fl_dist, prb) in fl_dists.values())
    min_fl = min(fl_dist.fl_min for (fl_dist, prb) in fl_dists.values())
    transcripts = set()
    segment_transcripts_map = {}
    for segment_id, data in splice_graph.nodes(data=True):
        if 'type' not in data: 
            assert False, str((gene, splice_graph.nodes(data=True), 
                               splice_graph.node[segment_id], segment_id, data))
        if data['type'] != 'segment': continue        
        segment_transcripts = find_transcribed_fragments_covering_region(
            splice_graph, segment_id, max_fl)
        segment_transcripts = [tuple(t) for t in segment_transcripts]
        segment_transcripts_map[(segment_id,)] = segment_transcripts
        transcripts.update(segment_transcripts)
    all_transcripts = sorted(transcripts)
    
    def bin_contains_element(bin, element):
        element = tuple(sorted(element))
        try: 
            si = bin.index(element[0])
            if tuple(bin[si:si+len(element)]) == element: 
                return True
        except ValueError: 
            return False
        return False
    
    # add in the splice elements
    for start_i, stop_i, data in splice_graph.edges(data=True):
        if data['type'] != 'splice': continue
        # find transcripts that contain this splice
        segment_transcripts_map[(start_i, stop_i)] = [
            t for t in segment_transcripts_map[(start_i,)]
            if bin_contains_element(t, tuple(sorted((start_i, stop_i))))]
    
    # build the old segment bnd, label style list. We need this for read binning
    # and expected fragment count calculations - this should probably be done
    # in the splice graph class, TODO
    segment_nodes = {}
    segment_bnds = set()
    segment_bnd_labels = defaultdict(set)
    for segment_i, data in sorted(splice_graph.nodes(data=True)):
        if data['type'] != 'segment': continue
        segment_bin = data['bin']
        segment_nodes[segment_i] = segment_bin
        segment_bnds.add(segment_bin.start)
        segment_bnd_labels[segment_bin.start].update(segment_bin.left_labels)
        segment_bnds.add(segment_bin.stop+1)
        segment_bnd_labels[segment_bin.stop+1].update(segment_bin.right_labels)
    segment_bnds = numpy.array(sorted(segment_bnds))
    segment_bnd_labels = dict(segment_bnd_labels)
    exon_lens = segment_bnds[1:] - segment_bnds[:-1]
    
    #print gene
    #print segment_bnds
    #print segment_bnd_labels

    weighted_expected_cnts = defaultdict(lambda: defaultdict(float))
    for (rg, (r1_len,r2_len)), (fl_dist, marginal_frac) in fl_dists.iteritems():
        for tr, bin_cnts in f_matrix.calc_expected_cnts( 
                segment_bnds, transcripts, 
                fl_dist, r1_len, r2_len).iteritems():
            for bin, cnt in bin_cnts.iteritems():
                weighted_expected_cnts[tr][bin] += cnt*marginal_frac
        
    segment_bins, jn_bins = [], []
    for j, (element, transcripts) in enumerate(
            segment_transcripts_map.iteritems()):
        print element, transcripts
        continue
    assert False
    if True:
        if config.VERBOSE:
            config.log_statement( 
                "Estimating element expression for segment "
                + "%i/%i in Chrm %s Strand %s Pos %i-%i" %
                (j, len(segment_transcripts_map), 
                 gene.chrm, gene.strand, gene.start, gene.stop) )

        # find the normalized, expected fragment counts for this element
        exp_bin_cnts_in_segment = defaultdict(
            lambda: numpy.zeros(len(transcripts), dtype=float))
        obs_bin_cnts_in_segment = defaultdict(int)
        effective_t_lens = numpy.zeros(len(transcripts), dtype=float)
        for i, transcript in enumerate(transcripts):
            for pe_bin, weighted_n_distinct_frags in weighted_expected_cnts[
                    transcript].iteritems():
                # skip bins that don't overlap the desired element
                if not any(bin_contains_element(bin,element) for bin in pe_bin):
                    continue
                effective_t_lens[i] += weighted_n_distinct_frags
                exp_bin_cnts_in_segment[pe_bin][i] += weighted_n_distinct_frags
                if (pe_bin not in obs_bin_cnts_in_segment 
                         and pe_bin in binned_reads):
                    obs_bin_cnts_in_segment[pe_bin] = binned_reads[pe_bin]
        #if len(element) == 1:
        #    for t, e_len in zip(transcripts, effective_t_lens):
        #        assert ( e_len <= exon_lens[element[0]] + max_fl 
        #                 - f_matrix.MIN_NUM_MAPPABLE_BASES ), str((
        #                     e_len, exon_lens[element[0]], max_fl,
        #                     f_matrix.MIN_NUM_MAPPABLE_BASES ))
        """ Code for degugging the element counts
        print element, exon_lens[element]
        all_bins = sorted(set(chain(*transcripts)))
        print [(x, exon_lens[x]) for x in all_bins ]
        print [(x, segment_bnd_labels[segment_bnds[x]], 
                   segment_bnd_labels[segment_bnds[x+1]])
               for x in all_bins]
        for t, e_len in zip(transcripts, effective_t_lens):
            print t, e_len
            #print ( t, e_len, exon_lens[element] + max_fl - 1, 
            #        exon_lens[element] - max_fl + 1 + 2*max_fl - 2 )
            #print sorted(weighted_expected_cnts[t].iteritems())
            #print
        raw_input()
        """
        
        # if all of the transcript fragments have the same length, we are done
        # all equally likely so we can estimate the element frequencies directly
        # from the bin counts
        if all(x == effective_t_lens[0] for x in effective_t_lens):
            mean_t_len = effective_t_lens[0]
        # otherwise we need to estimate the transcript frequencies
        elif True:
            mean_t_len = max(effective_t_lens)
        else:
            exp_a, obs_a, zero = f_matrix.build_expected_and_observed_arrays(
                exp_bin_cnts_in_segment, obs_bin_cnts_in_segment)
            exp_a, obs_a = f_matrix.cluster_rows(exp_a, obs_a)
            try: 
                t_freqs = frequency_estimation.estimate_transcript_frequencies(
                    obs_a, exp_a)
                mean_t_len = (t_freqs*effective_t_lens).sum()
            except frequency_estimation.TooFewReadsError: 
                mean_t_len = effective_t_lens.mean()

        assert mean_t_len > 0, str((element, transcripts, exon_lens, all_transcripts, sorted(segment_bnd_labels.iteritems()), splice_graph.edges(data=True)))
        n_reads_in_segment = float(sum(obs_bin_cnts_in_segment.values()))
        quantiles = [0.01, 0.5, 1-.01]
        fpkms = 1e6*(1000./mean_t_len)*beta.ppf(
            quantiles, 
            n_reads_in_segment+1e-6, 
            rnaseq_reads.num_reads-n_reads_in_segment+1e-6)
        if len(element) == 1:
            segment_nodes[element[0]].set_expression(*fpkms)
        else:
            splice_graph[element[0]][element[1]]['bin'].set_expression(*fpkms)
    
    return splice_graph

def determine_exon_type(left_label, right_label):
    if left_label == 'TSS':
        if right_label == 'TES': 
            return 'SE_GENE'
        else:
            assert right_label == 'D_JN'
            return 'TSS_EXON'

    if right_label == 'TES':
        # we should have alrady caught this case
        assert left_label != 'TES'
        assert left_label == 'R_JN'
        return 'TES_EXON'

    assert left_label == 'R_JN' and right_label == 'D_JN'
    return 'EXON'

def build_exons_from_exon_segments(gene, splice_graph, max_min_expression):
    config.log_statement( 
        "Building Exons from Segments in Chrm %s Strand %s Pos %i-%i" %
        (gene.chrm, gene.strand, gene.start, gene.stop) )

    exon_segments = [ data['bin']
                      for node_id, data in splice_graph.nodes(data=True)
                      if data['type'] == 'segment' ]
    if gene.strand == '-':
        exon_segments = reverse_strand(exon_segments, gene.stop)
    exon_segments = sorted(exon_segments, key=lambda x:x.start)
    
    EXON_START_LABELS = ('TSS', 'R_JN')
    EXON_STOP_LABELS = ('TES', 'D_JN')
    
    # find all of the exon start bins
    exons = []
    for i in xrange(len(exon_segments)):
        # if this is not allowed to be the start of an exon
        start_segment = exon_segments[i]
        for start_label in start_segment.left_labels:
            if start_label not in EXON_START_LABELS:
                continue

            for j in xrange(i, len(exon_segments)):
                stop_segment = exon_segments[j]

                #if ( start_label == 'D_JN'
                #     and 'R_JN' in stop_segment.right_labels ):
                #    break
                
                local_max_min_expression = max_min_expression
                if j+i < len(exon_segments):
                    local_max_min_expression = max(
                        local_max_min_expression, 
                        exon_segments[j+1].fpkm_lb/config.MAX_EXPRESSION_RATIO)
                if j-i >= 0:
                    local_max_min_expression = max(
                        local_max_min_expression, 
                        exon_segments[j-1].fpkm_lb/config.MAX_EXPRESSION_RATIO)
                if stop_segment.fpkm < local_max_min_expression: 
                    break
                
                for stop_label in stop_segment.right_labels:
                    if stop_label not in EXON_STOP_LABELS:
                        continue
                    fpkm = min(segment.fpkm for segment in exon_segments[i:j+1])
                    exon_bin = TranscriptElement(
                        start_segment.start, stop_segment.stop, 
                        determine_exon_type(start_label, stop_label),
                        fpkm)
                    exons.append(exon_bin)
            
    if gene.strand == '-':
        exons = reverse_strand(exons, gene.stop)
    
    return exons

def find_exons_in_gene( gene, contig_lens, ofp,
                        ref_elements, ref_elements_to_include,
                        rnaseq_reads, cage_reads, polya_reads ):
    # extract the reference elements that we want to add in
    gene_ref_elements = defaultdict(list)
    for key, vals in ref_elements[(gene.chrm, gene.strand)].iteritems():
        if len( vals ) == 0: continue
        for start, stop in sorted(vals):
            if stop < gene.start: continue
            if start > gene.stop: break
            if gene.base_is_in_gene(start) and gene.base_is_in_gene(stop):
                gene_ref_elements[key].append((start, stop))

    config.log_statement( "Finding Exons in Chrm %s Strand %s Pos %i-%i" %
                   (gene.chrm, gene.strand, gene.start, gene.stop) )
    
    # build the transcribed segment splice graph, and bin observe rnasseq reads 
    # based upon this splice graph in this gene.
    splice_graph, binned_reads = build_splice_graph_and_binned_reads_in_gene(
        gene, rnaseq_reads, cage_reads, polya_reads, ref_elements )
    splice_graph = fast_quantify_segment_expression(
        gene, splice_graph, rnaseq_reads, cage_reads, polya_reads )
    #splice_graph = quantify_segment_expression(gene, splice_graph, binned_reads)
    # build exons, and add them to the gene
    min_max_exp_t, max_element_exp = find_widest_path(splice_graph)
    min_max_exp = max(
        max_element_exp/config.MAX_EXPRESSION_RATIO, config.MIN_EXON_FPKM)
    exons = build_exons_from_exon_segments(gene, splice_graph, min_max_exp)
    gene.elements.extend(exons)
    
    # introns are both elements and element segments
    gene.elements.extend(
        data['bin'] for n1, n2, data in splice_graph.edges(data=True)
        if data['type'] == 'splice'
        and data['bin'].fpkm > min_max_exp)

    if config.DEBUG_VERBOSE:
        gene.elements.extend(
            data['bin'] for node_id, data in splice_graph.nodes(data=True)
            if data['bin'].fpkm > min_max_exp )
    
    # add T*S's
    gene.elements.extend(
        data['bin'] for node_id, data in splice_graph.nodes(data=True)
        if data['type'] in ('TSS', 'TES')
        and data['bin'].fpkm > min_max_exp)

    # merge in the reference exons
    for tss_exon in gene_ref_elements['tss_exon']:
        gene.elements.append( Bin(tss_exon[0], tss_exon[1], 
                                  "REF_TSS_EXON_START", "REF_TSS_EXON_STOP",
                                  "TSS_EXON") )
    for tes_exon in gene_ref_elements['tes_exon']:
        gene.elements.append( Bin(tes_exon[0], tes_exon[1], 
                                  "REF_TES_EXON_START", "REF_TES_EXON_STOP",
                                  "TES_EXON") )
    
    # add the gene bin
    gene.write_elements_bed(ofp)
    return
    config.log_statement( "FINISHED Finding Exons in Chrm %s Strand %s Pos %i-%i" %
                   (gene.chrm, gene.strand, gene.start, gene.stop) )
    return None

def find_exons_worker( (genes_queue, genes_queue_lock, n_threads_running), 
                       ofp, contig_lens, ref_elements, ref_elements_to_include,
                       rnaseq_reads, cage_reads, polya_reads ):
    rnaseq_reads = rnaseq_reads.reload()
    cage_reads = cage_reads.reload() if cage_reads != None else None
    polya_reads = polya_reads.reload() if polya_reads != None else None
    
    while True:
        # try to get a gene
        with genes_queue_lock:
            try: gene = genes_queue.pop()
            except IndexError: gene = None
        
        # if there is no gene it process, but threads are still running, then 
        # wait for the queue to fill or the process to finish
        if gene == None and n_threads_running.value > 0:
            config.log_statement( 
                "Waiting for gene to process (%i)" % n_threads_running.value)
            time.sleep(0.1)
            continue
        # otherwise, take a lock to make sure that the queue is empty and no 
        # threads are running. If so, then return
        elif gene == None:
            with genes_queue_lock:
                if len(genes_queue) == 0 and n_threads_running.value == 0:
                    config.log_statement( "" )
                    return
                else: continue

        # otherwise, we have a gene to process, so process it
        assert gene != None
        with genes_queue_lock: n_threads_running.value += 1

        # find the exons and genes
        try:
            rv = find_exons_in_gene(gene, contig_lens, ofp,
                                    ref_elements, ref_elements_to_include,
                                    rnaseq_reads, cage_reads, polya_reads )
        except Exception, inst:
            config.log_statement( 
                "Uncaught exception in find_exons_in_gene", log=True )
            config.log_statement( traceback.format_exc(), log=True, display=False )
            rv = None
        
        # if the return value is new genes, then add these to the queue
        if rv != None:
            with genes_queue_lock:
                for gene in rv:
                    genes_queue.append( gene )
                n_threads_running.value -= 1
        # otherwise, we found the exons and wrote out the element data,
        # so just decrease the number of running threads
        else:
            with genes_queue_lock:
                n_threads_running.value -= 1
    
    assert False
    return

def extract_reference_elements(genes, ref_elements_to_include):
    ref_elements = defaultdict( lambda: defaultdict(set) )
    if not any(ref_elements_to_include):
        return ref_elements
    
    for gene in genes:
        elements = gene.extract_elements()
        def add_elements(key):
            for start, stop in elements[key]:
                ref_elements[(gene.chrm, gene.strand)][key].add((start, stop))

        if ref_elements_to_include.junctions:
            add_elements('intron')
        if ref_elements_to_include.promoters:
            add_elements('promoter')
        if ref_elements_to_include.polya_sites:
            add_elements('polya')
        if ref_elements_to_include.TSS:
            add_elements('tss_exon')
        if ref_elements_to_include.TES:
            add_elements('tes_exon')
        if ref_elements_to_include.TES:
            add_elements('exon')
    
    for contig_strand, elements in ref_elements.iteritems():
        for element_type, val in elements.iteritems():
            ref_elements[contig_strand][element_type] = sorted( val )
    
    return ref_elements


def find_exons( contig_lens, gene_bndry_bins, ofp,
                rnaseq_reads, cage_reads, polya_reads,
                ref_genes, ref_elements_to_include,
                junctions=None, nthreads=None):
    assert not any(ref_elements_to_include) or ref_genes != None
    if nthreads == None: nthreads = config.NTHREADS
    assert junctions == None
    
    ref_elements = extract_reference_elements( 
        ref_genes, ref_elements_to_include )
    genes_queue_lock = multiprocessing.Lock()
    threads_are_running = multiprocessing.Value('i', 0)
    if config.NTHREADS > 1:
        manager = multiprocessing.Manager()
        genes_queue = manager.list()
    else:
        genes_queue = []

    genes_queue.extend( gene_bndry_bins )
    """
    for ref_gene in ref_genes:
        gene = GeneElements(ref_gene.chrm, ref_gene.strand)
        gene.regions.append(
            SegmentBin(ref_gene.start, ref_gene.stop, 
                       ["ESTART",], ["ESTOP",], "GENE"))
        genes_queue.append(gene)
    """
    args = [ (genes_queue, genes_queue_lock, threads_are_running), 
             ofp, contig_lens, ref_elements, ref_elements_to_include,
             rnaseq_reads, cage_reads, polya_reads ]
    
    n_genes = len(genes_queue)
    if nthreads == 1:
        find_exons_worker(*args)
    else:
        config.log_statement("Waiting on exon finding children (%i/%i remain)"%(
                len(genes_queue), n_genes))
        ps = []
        for i in xrange( nthreads ):
            p = multiprocessing.Process(target=find_exons_worker, args=args)
            p.start()
            ps.append( p )
        
        while True:
            config.log_statement(
                "Waiting on exon finding children (%i/%i remain)"%(
                    len(genes_queue), n_genes))
            if all( not p.is_alive() for p in ps ):
                break
            time.sleep( 1.0 )

    config.log_statement( "" )    
    return

def find_elements( promoter_reads, rnaseq_reads, polya_reads,
                   ofname, ref_genes, ref_elements_to_include, 
                   gene_segments, region_to_use=None):
    # wrap everything in a try block so that we can with elegantly handle
    # uncaught exceptions
    try:
        ofp = ThreadSafeFile( ofname + "unfinished", "w" )
        ofp.write(
            'track name="%s" visibility=2 itemRgb="On" useScore=1\n' % ofname)
        
        contig_lens = dict(zip(*get_contigs_and_lens( 
            [ reads for reads in [rnaseq_reads, promoter_reads, polya_reads]
              if reads != None ] )))
        
        """
        # extract reference elements. This doesnt work for a few reasons,
        # so we merge in reference genes by setting exon and jn refernce
        # elements to True
        if ref_elements_to_include.genes == True:
            gene_segments = []
            for contig, contig_len in contig_lens.iteritems():
                for strand in '+-':
                    contig_gene_bndry_bins = load_gene_bndry_bins(
                        ref_genes, contig, strand, contig_len)
                    gene_segments.extend( contig_gene_bndry_bins )
            fl_dists = build_fl_dists_from_fls_dict(
                find_fls_from_annotation(ref_genes, rnaseq_reads))
        """
        
        # sort genes from longest to shortest. This should help improve the 
        # multicore performance
        gene_segments.sort( key=lambda x: x.stop-x.start, reverse=True )
        find_exons( contig_lens, gene_segments, ofp,
                    rnaseq_reads, promoter_reads, polya_reads,
                    ref_genes, ref_elements_to_include, 
                    junctions=None, nthreads=config.NTHREADS )            
    except Exception, inst:
        config.log_statement( "FATAL ERROR", log=True )
        config.log_statement( traceback.format_exc(), log=True, display=False )
        ofp.close()
        raise
    else:
        ofp.close()
    
    shutil.move(ofname + "unfinished", ofname)
    return ofname
