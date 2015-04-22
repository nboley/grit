import os, sys

from collections import defaultdict
from itertools import chain

from grit.files.gtf import load_gtf
from grit.files.reads import RNAseqReads, MergedReads
from grit.frag_len import build_fl_dists_from_annotation
from grit.f_matrix import build_expected_and_observed_rnaseq_counts, cluster_bins, DesignMatrix, NoObservableTranscriptsError
from grit.frequency_estimation import estimate_transcript_frequencies, TooFewReadsError

chrm = '20'
start = 0
stop = 100000000

def bin_and_cluster_reads(gene, reads, fl_dists):
    ( expected_rnaseq_cnts, observed_rnaseq_cnts 
     ) = build_expected_and_observed_rnaseq_counts( 
         gene, reads, fl_dists )

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
        rv.append([sorted(segments), exp_cnt, observed_cluster_cnts[cluster]])
    if len(rv) == 0:
        return [], [], []
    return zip(*rv)

def main():
    genes = load_gtf(sys.argv[1])

    reads = MergedReads([
        RNAseqReads(fname).init(reads_are_stranded=False)
        for fname in sys.argv[2:]] )

    fl_dists = build_fl_dists_from_annotation( genes, reads )
    
    for gene in genes.iter_overlapping_genes(chrm, '.', start, stop):
        print gene.name, gene.chrm, gene.strand, gene.start, gene.stop
        cluster_segments, exp_cnts, obs_cnts = bin_and_cluster_reads(
            gene, reads, fl_dists)
    
    return

main()
