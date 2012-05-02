import pickle
import reads
import gene_models
from frag_len import build_robust_fl_dist_with_stats

def get_frag_lengths( reads, genes ):
    fragment_lengths = []
    for gene in genes.values():
        if len( gene.exon_bndrys ) != 1:
            continue
        if ( gene.exon_bndrys[0][1] - gene.exon_bndrys[0][0]) < 700:
            continue

        # iterate through read pairs in exon to get fragment lengths
        for read1, read2 in reads.iter_paired_reads( gene.boundaries ):
            # skip all secondary alignments
            if read1.is_secondary or read2.is_secondary:
                continue

            strand = '-' if read1.is_reverse else '+'
            if not read1.is_reverse:
                frag_len = read2.aend - read1.pos
            else:
                frag_len = read1.aend - read2.pos

            fragment_lengths.append( frag_len )

    return fragment_lengths

def estimate_fl_dists( reads, genes ):
    fragment_lengths = get_frag_lengths( reads, genes )
    fl_dist = build_robust_fl_dist_with_stats( fragment_lengths )
    print "Fragment Length Distribution Statistics:\n\tMean:\t\t\t{0[0][0]:.4f} +/- {0[1][0]:.4f}\n\t\
Standard Deviation:\t{0[0][1]:.4f} +/- {0[1][1]:.4f}\n\tSkew:\t\t\t{0[0][2]:.4f} +/- {0[1][2]:.4f}\n".format( fl_dist.stats )

    return fl_dist

def build_objs( gtf_fp, bam_fn ):
    genes = gene_models.GeneBoundaries( gtf_fp )
    gtf_fp.close()

    reads_obj = reads.Reads( bam_fn )

    return genes, reads_obj

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Produce a cached fragment length distribution ' + \
                                         'object for use with reads_simulator.py.' )
    parser.add_argument( 'gtf', type=file, \
                             help='GTF file of high read count single exon genes ' + \
                             'over which to estimate fragment length ditribution.')
    parser.add_argument( 'bam_fn',\
                             help='BAM file for which to estimate fragment length distribution.')
    args = parser.parse_args()

    return args.gtf, args.bam_fn

if __name__ == "__main__":
    seg_fp, reads_fn = parse_arguments()
    single_exon_genes, reads = build_objs( seg_fp, reads_fn )

    fl_dist = estimate_fl_dists( reads, single_exon_genes )
    with open( 'fl_dist.obj', "w" ) as fp:
        pickle.dump( fl_dist , fp )
