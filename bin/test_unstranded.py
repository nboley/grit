import os, sys

import numpy

from collections import defaultdict
from itertools import chain

import multiprocessing
import queue

from grit.genes import merge_adjacent_intervals
from grit.files.gtf import load_gtf
from grit.files.reads import RNAseqReads, MergedReads
from grit.frag_len import build_fl_dists_from_annotation
from grit.f_matrix import build_expected_and_observed_rnaseq_counts, cluster_bins, DesignMatrix, NoObservableTranscriptsError
from grit.frequency_estimation import estimate_transcript_frequencies, TooFewReadsError
from grit.transcript import Gene, Transcript
from grit.lib.multiprocessing_utils import ThreadSafeFile, fork_and_wait

def bin_and_cluster_reads(gene, reads, fl_dists):
    ( expected_rnaseq_cnts, observed_rnaseq_cnts 
     ) = build_expected_and_observed_rnaseq_counts( 
         gene, reads, fl_dists )
    exon_boundaries = numpy.array(gene.find_nonoverlapping_boundaries())
    
    clustered_bins = cluster_bins(expected_rnaseq_cnts)
    expected_cluster_cnts = {}
    observed_cluster_cnts = {}
    cluster_components = {}
    for i, cluster in enumerate(sorted(clustered_bins)):
        expected_cluster_cnts[i] = 0.0
        observed_cluster_cnts[i] = 0.0
        cluster_components[i] = []
        for bin in cluster:
            expected_cluster_cnts[i] += list(expected_rnaseq_cnts[bin].values())[0]
            try: observed_cluster_cnts[i] += observed_rnaseq_cnts[bin]
            except: pass
            cluster_components[i].append(bin)

    rv = []
    for cluster, exp_cnt in list(expected_cluster_cnts.items()):
        segments = set()
        for rd_len, read_grp, (r1_bin, r2_bin) in cluster_components[cluster]:
            segments.update(r1_bin)
            segments.update(r2_bin)
        segments = sorted(segments)
        regions = [(exon_boundaries[i], exon_boundaries[i+1]) for i in segments]
        regions = merge_adjacent_intervals(regions, 0)
        transcript = Transcript(
            "%s_%i" % (gene.id, cluster), 
            'chr' + gene.chrm, 
            gene.strand,
            regions, 
            cds_region=None,
            gene_id = gene.id )
        rv.append([transcript, exp_cnt, observed_cluster_cnts[cluster]])
    
    if len(rv) == 0:
        return [], numpy.array([]), numpy.array([])
    
    segments, exp, obs = list(zip(*rv))
    
    return segments, numpy.array(exp), numpy.array(obs)

def build_transcripts_lines_worker(genes_queue, reads, fl_dists, ofp):
    reads.reload()
    n_genes = genes_queue.qsize()
    while not genes_queue.empty():
        try: gene = genes_queue.get(0.1)
        except queue.Empty: break
        print("%i/%i\t%s\t%s:%s:%i-%i" %  (
            n_genes-genes_queue.qsize(), n_genes,
            gene.name.ljust(20), gene.chrm, gene.strand, gene.start, gene.stop), file=sys.stderr)
    
        transcrtipt_segments, exp_cnts, obs_cnts = bin_and_cluster_reads(
            gene, reads, fl_dists)
        for transcript, exp_cnt, obs_cnt in zip(
                transcrtipt_segments, exp_cnts, obs_cnts):
            # filter out all bins with expected fragment counts less than 1
            # (these are very short bins that can only rarely be observed 
            # given our fragment length distribution) 
            if exp_cnt < 1: continue
            lines = transcript.build_gtf_lines(
                    {'effective_length': round(exp_cnt, 2), 
                     'observed_cnt': round(obs_cnt, 1)})
            ofp.write(lines + "\n")
    return

def main():
    genes = load_gtf(sys.argv[1])

    reads = MergedReads([
        RNAseqReads(fname).init(reads_are_stranded=False)
        for fname in sys.argv[2:]] )

    fl_dists = build_fl_dists_from_annotation( genes, reads )

    genes_queue = multiprocessing.Queue()
    for gene in genes: genes_queue.put(gene)    
    
    with ThreadSafeFile("bins.gtf", "w") as ofp:
        ofp.write("track type=gtf name=bin_expresion_test\n")
        args = [genes_queue, reads, fl_dists, ofp]
        fork_and_wait(24, build_transcripts_lines_worker, args)
    
    return

main()
