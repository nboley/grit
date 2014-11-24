import os, sys
import cPickle as pickle
import multiprocessing

sys.path.insert(0, "/home/nboley/grit/grit/")
from grit.merge import (
    group_overlapping_genes, 
    reduce_gene_clustered_transcripts, 
    fix_chrm_name_for_ucsc)

from grit.lib.multiprocessing_utils import ThreadSafeFile
from grit.files.gtf import load_multiple_gtfs_into_pickled_files
import grit.config as config

def parse_arguments():
    import argparse
    desc = 'Merge transcripts.'
    parser = argparse.ArgumentParser(description=desc)    
    parser.add_argument(
        "gtfs", type=file, nargs="+",
        help='GTF file contraining transcripts. (i.e. from build_transcripts)' )
    
    parser.add_argument(
        '--out-fname', '-o', 
        type=lambda fname: ThreadSafeFile(fname, 'w'), 
        default=sys.stdout,
        help='Output file. default: stdout')
    
    parser.add_argument(
        '--out-sources-fname', 
        type=lambda fname: ThreadSafeFile(fname, 'w'),
        help='File name to write the sources for each transcript. '
        + 'Default: Do not write out sources map.')
    
    parser.add_argument( '--ucsc', default=False, action='store_true',
        help='Try to format contig names in the ucsc format (typically by prepending a chr).')    
    
    parser.add_argument(
        '--min-fpkm-ub', type=float, 
        help='Filter transcripts with fpkms upper bounds below this value.')

    parser.add_argument(
        '--intrasample-max-fpkm-ratio', type=float, 
        help='For each gene cluster and sample, filter transcripts whose fpkm upper bound is this many times lower than the maximum fpkm lower bound.')

    parser.add_argument(
        '--intersample-max-fpkm-ratio', type=float, 
        help='For each gene cluster and between all samples, filter transcripts whose fpkm upper bound is this many times lower than the maximum fpkm lower bound.')
        
    parser.add_argument(
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
 
    parser.add_argument(
        '--threads', '-t', type=int, default=1,
        help='The number of threads to use.')
   
    args = parser.parse_args()
    
    config.VERBOSE = args.verbose
    
    config.NTHREADS = args.threads

    global min_upper_fpkm
    min_upper_fpkm = args.min_fpkm_ub
    global max_intrasample_fpkm_ratio
    max_intrasample_fpkm_ratio = args.intrasample_max_fpkm_ratio
    global max_intersample_fpkm_ratio
    max_intersample_fpkm_ratio = args.intersample_max_fpkm_ratio
    
    config.FIX_CHRM_NAMES_FOR_UCSC = args.ucsc
    
    return args.gtfs, args.out_fname, args.out_sources_fname

def write_reduced_gene_to_file(new_gene, merged_transcript_sources, 
                               ofp, tracking_ofp):
    gtf_lines = []
    tracking_lines = []
    for merged_transcript, sources in zip(
            new_gene.transcripts, merged_transcript_sources):
        if config.FIX_CHRM_NAMES_FOR_UCSC:
            merged_transcript.chrm = fix_chrm_name_for_ucsc(
                merged_transcript.chrm)
        gtf_lines.append(
            merged_transcript.build_gtf_lines(meta_data={}, source='GRIT'))
        for old_t, source_key in sources:
            line = "\t".join(
                (old_t.id, source_key, merged_transcript.id))
            tracking_lines.append(line)
    
    ofp.write("\n".join(gtf_lines)+"\n")
    if tracking_ofp != None:
        tracking_ofp.write("\n".join(tracking_lines)+"\n")
    return

def merge_clustered_genes_worker(clustered_genes, clustered_genes_lock, 
                                 ofp, tracking_ofp,
                                 gene_id_cntr):
    while True:
        with clustered_genes_lock:
            if len(clustered_genes) == 0: 
                return
            genes = clustered_genes.pop()
        with gene_id_cntr.get_lock():
            gene_id_cntr.value += 1
            new_gene_id = "GENE_%i" % gene_id_cntr.value
        new_gene, merged_transcript_sources = reduce_gene_clustered_transcripts(
            genes, new_gene_id, 
            min_upper_fpkm=min_upper_fpkm,
            max_intrasample_fpkm_ratio=max_intrasample_fpkm_ratio, 
            max_intersample_fpkm_ratio=max_intersample_fpkm_ratio,
            max_cluster_gap=500)
        write_reduced_gene_to_file(
            new_gene, merged_transcript_sources, ofp, tracking_ofp)
    
    return

def merge_genes(all_sources_and_genes, ofp, sources_ofp):
    # write the gtf header
    ofp.write("track name=%s\n" % ofp.name)
    
    # group overlapping genes
    config.log_statement("Grouping genes", log=True)
    manager = multiprocessing.Manager()
    grpd_genes = manager.list()
    grpd_genes_lock = multiprocessing.Lock()

    for genes in group_overlapping_genes(all_sources_and_genes):
        grpd_genes.append(genes)
    
    # merge gene clustered transcripts
    config.log_statement("Merging transcripts", log=True)
    
    gene_id_cntr = multiprocessing.Value('i', 0)

    if config.NTHREADS == 1:
        merge_clustered_genes_worker(
            grpd_genes, grpd_genes_lock, 
            ofp, sources_ofp, gene_id_cntr)
    else:
        pids = []
        for i in xrange(config.NTHREADS):
            pid = os.fork()
            if pid == 0:
                merge_clustered_genes_worker(
                    grpd_genes, grpd_genes_lock, 
                    ofp, sources_ofp, gene_id_cntr)
                os._exit(0)
            else:
                pids.append(pid)
        for pid in pids:
            os.waitpid(pid, 0)

    return

def main():
    gtf_fps, ofp, sources_ofp = parse_arguments()
    gtf_fnames = [os.path.abspath(fp.name) for fp in gtf_fps]
        
    config.log_statement("Loading gtfs")
    all_genes_and_fnames = load_multiple_gtfs_into_pickled_files(gtf_fnames)
    
    merge_genes(all_genes_and_fnames, ofp, sources_ofp)
    
    ofp.close()
    if sources_ofp != None:
        sources_ofp.close()

main()
