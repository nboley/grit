import os, sys
import cPickle as pickle
import multiprocessing
import copy
import random

from itertools import izip, chain
from collections import defaultdict

import numpy
from scipy.cluster.hierarchy import fclusterdata

sys.path.insert(0, "/home/nboley/grit/grit/")
from grit.files.gtf import load_gtf_into_pickled_files
from grit.files.reads import fix_chrm_name_for_ucsc
from grit.transcript import Transcript
from grit.lib.multiprocessing_utils import ThreadSafeFile
import grit.config as config

NTHREADS = 32
VERBOSE = True
MERGE_DISTAL_ENDS = False
SINGLE_LINKAGE_CLUSTER = True
LINKAGE_CLUSTER_GAP = 500
LINKAGE_CLUSTER_MAX_DIFF = 2000
FIX_CHRM_NAMES_FOR_UCSC = True

min_upper_fpkm = 0
max_intrasample_fpkm_ratio = 1e6
max_intersample_fpkm_ratio = 1e6

def build_merged_transcript(gene_id, clustered_transcripts):
    # find hte transcript bounds
    start, stop = 1e20, 0
    for transcript in clustered_transcripts:
        start = min( start, transcript.exons[0][0] )
        stop = max( stop, transcript.exons[-1][-1] )
    
    # merge the promoters
    try:
        new_promoter = ( min(t.promoter[0] for t in clustered_transcripts 
                             if t.promoter != None),
                         max(t.promoter[1] for t in clustered_transcripts
                             if t.promoter != None) )
    except ValueError:
        new_promoter = None
    
    # merge the polyas
    try:
        new_polya = ( min(t.polya_region[0] for t in clustered_transcripts 
                          if t.polya_region != None),
                      max(t.polya_region[1] for t in clustered_transcripts
                          if t.polya_region != None) )
    except ValueError:
        new_polya = None
    
    # choose a tempalte transcript, and make sure that all of the 
    # clustered transcripts have the same internal structure (
    # this should be guaranteed by the calling function )
    bt = clustered_transcripts[0]
    assert all( t.IB_key() == bt.IB_key() for t in clustered_transcripts )
    new_exons = list(bt.exons)
    new_exons[0] = ( start, new_exons[0][1] )
    new_exons[-1] = ( new_exons[-1][0], stop )
    # choose a random id - this should be renamed in the next step
    new_trans_id = gene_id + "_RNDM_%i" % random.randrange(1e9)
    new_transcript = Transcript( new_trans_id,
                                 bt.chrm, bt.strand, 
                                 new_exons, bt.cds_region, gene_id,
                                 name=bt.name,
                                 gene_name=bt.gene_name,
                                 promoter=new_promoter,
                                 polya_region=new_polya)
    
    return new_transcript

def reduce_internal_clustered_transcripts( internal_grpd_transcripts, gene_id ):
    """Take a set of clustered transcripts and reduce them into 
    a set of canonical transcripts, and associated sources.
    
    """
    # if there is only a single trnascript, clustering doesnt make sense
    if len( internal_grpd_transcripts ) == 1: 
        new_t = copy.copy(internal_grpd_transcripts[0][0])
        new_t.gene_id = gene_id
        new_t.id = new_t.gene_id + "_1"
        yield ( new_t, 
                [internal_grpd_transcripts[0][0],], 
                [internal_grpd_transcripts[0][1],] )
        return

    # 2 transcripts are in the same cluster if both their 5' and 3' ends
    # are within 50 bp's of each other. Use the scipy cluster machinery 
    # to do this for us
    transcript_ends = numpy.array( [(t.exons[0][0], t.exons[-1][1])
                                    for t, s in internal_grpd_transcripts])    
    cluster_indices = fclusterdata( transcript_ends, t=LINKAGE_CLUSTER_GAP,
                                    criterion='distance', metric='chebyshev' )
    
    # convert the incdices returned by flclusterdata into lists of transcript
    # source pairs
    clustered_transcript_grps = defaultdict( list )
    clustered_transcript_grp_sources = defaultdict( list )
    for cluster_index, ( trans, src ) in \
            izip(cluster_indices, internal_grpd_transcripts):
        clustered_transcript_grps[cluster_index].append( trans )
        clustered_transcript_grp_sources[cluster_index].append( src )
    
    # finally, decide upon the 'canonical' transcript for each cluster, and 
    # add it and it's sources
    for cluster_index in clustered_transcript_grps.keys():
        clustered_transcripts = clustered_transcript_grps[cluster_index]
        clustered_transcripts_sources = clustered_transcript_grp_sources[
            cluster_index]
        merged_transcript = build_merged_transcript(
                gene_id, clustered_transcripts)
        yield ( merged_transcript, 
                clustered_transcripts, 
                clustered_transcripts_sources)
    
    return

def reduce_gene_clustered_transcripts(genes, gene_id_cntr,
                                      all_genes_and_fnames):
    with gene_id_cntr.get_lock():
        gene_id_cntr.value += 1
        new_gene_id = "GENE_%i" % gene_id_cntr.value

    # unpickle the genes, and find the expression thresholds    
    unpickled_genes = []
    max_fpkm_lb_across_samples = 0
    max_fpkm_lb_in_sample = [0.0]*len(all_genes_and_fnames)
    
    for sample_i, gene_id in genes:
        region, fname = all_genes_and_fnames[sample_i][gene_id]
        with open(fname) as fp:
            gene = pickle.load(fp)
        unpickled_genes.append((sample_i, gene))
        try: 
            max_fpkm_lb_in_gene = max( 
                transcript.conf_lo for transcript in gene.transcripts
                if transcript.conf_lo != None )
        except ValueError:
                continue
        max_fpkm_lb_across_samples = max( 
                max_fpkm_lb_across_samples, max_fpkm_lb_in_gene )
        max_fpkm_lb_in_sample[sample_i] = max(
                max_fpkm_lb_in_gene, max_fpkm_lb_in_sample[sample_i])
    genes = unpickled_genes
    del unpickled_genes
    
    # group transcript by their internal structure
    internal_clustered_transcript_groups = defaultdict( list )
    for sample_i, gene in genes:
        min_max_fpkm = max(
                min_upper_fpkm, 
                max_fpkm_lb_in_sample[sample_i]/max_intrasample_fpkm_ratio,
                max_fpkm_lb_across_samples/max_intersample_fpkm_ratio)
        
        region, fname = all_genes_and_fnames[sample_i][gene.id]
        for transcript in gene.transcripts:
            if transcript.conf_hi < min_max_fpkm: continue
            IB_key = transcript.IB_key()
            internal_clustered_transcript_groups[IB_key].append(
                (transcript, sample_i))
    
    # reduce the non single exon genes
    reduced_transcripts = []
    transcript_id = 1
    for internal_clustered_transcripts in \
            internal_clustered_transcript_groups.values():
        for (merged_transcript, old_transcripts, sources
             ) in reduce_internal_clustered_transcripts( 
            internal_clustered_transcripts, new_gene_id ):
             merged_transcript.id = new_gene_id + "_%i" % transcript_id
             transcript_id += 1
             reduced_transcripts.append(
                (merged_transcript, new_gene_id, old_transcripts, sources))
    
    return reduced_transcripts

def write_reduced_gene_to_file(reduced_transcripts, fnames, ofp, tracking_ofp):
    gtf_lines = []
    tracking_lines = []
    for merged_transcript, gene_id, old_ts, sources in reduced_transcripts:
        if FIX_CHRM_NAMES_FOR_UCSC:
            merged_transcript.chrm = fix_chrm_name_for_ucsc(
                merged_transcript.chrm)
        gtf_lines.append(
            merged_transcript.build_gtf_lines(meta_data={}, source='GRIT'))
        for old_t, source_i in zip(old_ts, sources):
            line = "\t".join(
                (old_t.id, fnames[source_i], merged_transcript.id))
            tracking_lines.append(line)
    ofp.write("\n".join(gtf_lines)+"\n")
    tracking_ofp.write("\n".join(tracking_lines)+"\n")
    return

def merge_clustered_genes_worker(clustered_genes, clustered_genes_lock, 
                                 all_genes_and_fnames, fnames,
                                 ofp, tracking_ofp,
                                 gene_id_cntr):
    while True:
        with clustered_genes_lock:
            if len(clustered_genes) == 0: 
                return
            genes = clustered_genes.pop()
        reduced_transcripts = reduce_gene_clustered_transcripts(
            genes, gene_id_cntr, all_genes_and_fnames )
        write_reduced_gene_to_file(
            reduced_transcripts, fnames, ofp, tracking_ofp)
   
    return


def group_overlapping_genes(all_genes_and_fnames):
    chrm_grpd_genes = defaultdict(list)
    for sample_i, genes_and_fnames in enumerate(all_genes_and_fnames):
        for gene_id, (region, fname) in genes_and_fnames.iteritems():
            chrm, strand, start, stop = region
            chrm_grpd_genes[(chrm, strand)].append(
                (start, stop, sample_i, gene_id))
    
    grpd_genes = []
    for (chrm, strand), gene_regions in chrm_grpd_genes.iteritems():
        gene_regions.sort()
        curr_stop = -1
        for start, stop, sample_i, gene_id in gene_regions:
            if start > curr_stop: 
                grpd_genes.append([])
            grpd_genes[-1].append((sample_i, gene_id))
            curr_stop = stop
    
    return grpd_genes

def load_all_gtfs(fnames):
    # load all the gtfs
    config.log_statement("Loading gtfs")
    manager = multiprocessing.Manager()
    all_genes_and_fnames = manager.list()
    all_genes_and_fnames_lock = multiprocessing.Lock()
    pids = []
    for fname in fnames:
        pid = os.fork()
        if pid == 0:
            config.log_statement("Loading %s" % fname)
            sample_type = os.path.basename(fname).split('.')[0]
            expression_fnames = [ 
                    os.path.join(os.path.dirname(fname), f) 
                    for f in os.listdir(os.path.dirname(fname)) 
                    if os.path.basename(f).startswith(sample_type)
                    and os.path.basename(f).endswith("expression_tracking") ]
            genes_and_fnames = load_gtf_into_pickled_files(
                    fname, expression_fnames=expression_fnames)
            with all_genes_and_fnames_lock:
                all_genes_and_fnames.append(genes_and_fnames)
            config.log_statement("FINISHED Loading %s" % fname)
            os._exit(0)
        else:
            pids.append(pid)
    
    for pid in pids:
        os.waitpid(pid, 0)
    del all_genes_and_fnames_lock
    all_genes_and_fnames = list(all_genes_and_fnames)
    manager.shutdown()
    
    return all_genes_and_fnames

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
    
    parser.add_argument(
        '--min-fpkm-ub', default=0.0, type=float,
        help='Filter transcripts with fpkms upper bounds below this value.')

    parser.add_argument(
        '--intrasample-max-fpkm-ratio', type=int,
        help='For each gene cluster and sample, filter transcripts whose fpkm upper bound is this many times lower than the maximum fpkm lower bound.')

    parser.add_argument(
        '--intersample-max-fpkm-ratio', type=int,
        help='For each gene cluster and between all samples, filter transcripts whose fpkm upper bound is this many times lower than the maximum fpkm lower bound.')
        
    parser.add_argument(
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
 
    parser.add_argument(
        '--threads', '-t', type=int, default=1,
        help='The number of threads to use.')
   
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    global NTHREADS
    NTHREADS = args.threads

    global min_upper_fpkm
    min_upper_fpkm = args.min_fpkm_ub
    global max_intrasample_fpkm_ratio
    max_intrasample_fpkm_ratio = args.intrasample_max_fpkm_ratio
    global max_intersample_fpkm_ratio
    max_intersample_fpkm_ratio = args.intersample_max_fpkm_ratio
    
    return args.gtfs, args.out_fname, args.out_sources_fname

def main():
    gtf_fps, ofp, sources_ofp = parse_arguments()
    gtf_fnames = [os.path.abspath(fp.name) for fp in gtf_fps]
        
    all_genes_and_fnames = load_all_gtfs(gtf_fnames)
    
    # write hte gtf header
    ofp.write("track name=%s\n" % ofp.name)
    
    # group overlapping genes
    config.log_statement("Grouping genes", log=True)
    manager = multiprocessing.Manager()
    grpd_genes = manager.list()
    grpd_genes_lock = multiprocessing.Lock()
    for genes in group_overlapping_genes(all_genes_and_fnames):
        grpd_genes.append(genes)
    
    # merge gene clustered transcripts
    config.log_statement("Merging transcripts", log=True)
    
    gene_id_cntr = multiprocessing.Value('i', 0)

    if False and NTHREADS == 1:
        merge_clustered_genes_worker(
            grpd_genes, grpd_genes_lock, 
            all_genes_and_fnames, gtf_fnames,
            ofp, sources_ofp, gene_id_cntr)
    else:
        pids = []
        for i in xrange(NTHREADS):
            pid = os.fork()
            if pid == 0:
                merge_clustered_genes_worker(
                    grpd_genes, grpd_genes_lock, 
                    all_genes_and_fnames, gtf_fnames,
                    ofp, sources_ofp, gene_id_cntr)
                os._exit(0)
            else:
                pids.append(pid)
        for pid in pids:
            os.waitpid(pid, 0)
        
    ofp.close()
    sources_ofp.close()

main()
