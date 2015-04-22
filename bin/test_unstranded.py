import os, sys
import grit

from grit.files.gtf import load_gtf
from grit.files.reads import RNAseqReads, MergedReads
from grit.frag_len import build_normal_density, build_fl_dists_from_annotation
from grit.f_matrix import build_expected_and_observed_rnaseq_counts, build_expected_and_observed_arrays, DesignMatrix, NoObservableTranscriptsError
from grit.frequency_estimation import estimate_transcript_frequencies

chrm = '4'
start = 0
stop = 100000

def main():
    genes = load_gtf(sys.argv[1])

    reads = MergedReads([
        RNAseqReads(fname).init(reads_are_stranded=False)
        for fname in sys.argv[2:]] )

    fl_dists = build_fl_dists_from_annotation( genes, reads )
    
    genes = load_gtf(sys.argv[1])
    for gene in genes.iter_overlapping_genes(chrm, '.', start, stop):
        try: design_matrix = DesignMatrix(gene, fl_dists, reads, None, None)
        except NoObservableTranscriptsError: continue
        exp, obs =  design_matrix.expected_and_observed()
        print "Bin Cnts:    ", obs
        print "Exp Bin Cnts:", exp.T
        print "MLE freqs:   ", estimate_transcript_frequencies(obs, exp)
        print

    return

main()
