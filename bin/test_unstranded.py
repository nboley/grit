import os, sys

import numpy

from collections import defaultdict
from itertools import chain

from grit.genes import merge_adjacent_intervals
from grit.files.gtf import load_gtf
from grit.files.reads import RNAseqReads, MergedReads
from grit.frag_len import build_fl_dists_from_annotation
from grit.f_matrix import build_expected_and_observed_rnaseq_counts, cluster_bins, DesignMatrix, NoObservableTranscriptsError
from grit.frequency_estimation import estimate_transcript_frequencies, TooFewReadsError
from grit.transcript import Gene, Transcript


chrm = '4'
start = 0
stop = 100000000

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
            expected_cluster_cnts[i] += expected_rnaseq_cnts[bin].values()[0]
            try: observed_cluster_cnts[i] += observed_rnaseq_cnts[bin]
            except: pass
            cluster_components[i].append(bin)

    rv = []
    for cluster, exp_cnt in expected_cluster_cnts.items():
        segments = set()
        for rd_len, read_grp, (r1_bin, r2_bin) in cluster_components[cluster]:
            segments.update(r1_bin)
            segments.update(r2_bin)
        segments = sorted(segments)
        regions = [(exon_boundaries[i], exon_boundaries[i+1]) for i in segments]
        regions = merge_adjacent_intervals(regions, 0)
        transcript = Transcript(
            "%s_%i" % (gene.id, cluster), 
            gene.chrm, 
            gene.strand,
            regions, 
            cds_region=None,
            gene_id = gene.id )
        rv.append([transcript, exp_cnt, observed_cluster_cnts[cluster]])
    
    if len(rv) == 0:
        return [], numpy.array([]), numpy.array([])
    
    segments, exp, obs = zip(*rv)
    
    return segments, numpy.array(exp), numpy.array(obs)

def main():
    genes = load_gtf(sys.argv[1])

    reads = MergedReads([
        RNAseqReads(fname).init(reads_are_stranded=False)
        for fname in sys.argv[2:]] )

    fl_dists = build_fl_dists_from_annotation( genes, reads )
    
    for gene in genes.iter_overlapping_genes(chrm, '.', start, stop):
        print >> sys.stderr, gene.name, gene.chrm, gene.strand, gene.start, gene.stop
        transcrtipt_segments, exp_cnts, obs_cnts = bin_and_cluster_reads(
            gene, reads, fl_dists)
        for transcript, exp_cnt, obs_cnt in zip(
                transcrtipt_segments, exp_cnts, obs_cnts):
            lines = transcript.build_gtf_lines(
                    {'expected_cnt': exp_cnt, 'observed_cnt': obs_cnt})
            print lines
    
    return

main()
