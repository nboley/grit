import os, sys
import cPickle as pickle

from collections import defaultdict

import numpy
from scipy.cluster.hierarchy import fclusterdata

from itertools import izip, chain

import copy

sys.path.insert(0, "/home/nboley/grit/grit/")
from grit.files.gtf import load_gtf_into_pickled_files
from grit.files.reads import fix_chrm_name_for_ucsc
from grit.transcript import Transcript
from grit.lib.multiprocessing_utils import ThreadSafeFile

VERBOSE = True
MERGE_DISTAL_ENDS = False
SINGLE_LINKAGE_CLUSTER = True
LINKAGE_CLUSTER_GAP = 500
LINKAGE_CLUSTER_MAX_DIFF = 2000
FIX_CHRM_NAMES_FOR_UCSC = True

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
    for new_t_id, ( cluster_index, clustered_transcripts) in \
            enumerate(clustered_transcript_grps.iteritems()):
            
        start, stop = 1e20, 0
        for transcript in clustered_transcripts:
            start = min( start, transcript.exons[0][0] )
            stop = max( stop, transcript.exons[-1][-1] )
        try:
            new_promoter = ( min(t.promoter[0] for t in clustered_transcripts 
                                 if t.promoter != None),
                             max(t.promoter[1] for t in clustered_transcripts
                                 if t.promoter != None) )
        except ValueError:
            new_promoter = None
        
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
        new_transcript = Transcript( gene_id + "_%i" % new_t_id, 
                                     bt.chrm, bt.strand, 
                                     new_exons, bt.cds_region, gene_id,
                                     name=bt.name,
                                     gene_name=bt.gene_name,
                                     promoter=new_promoter,
                                     polya_region=new_polya)
        yield ( new_transcript, 
                clustered_transcript_grps[cluster_index],
                clustered_transcript_grp_sources[cluster_index] )
                
    
    return

def reduce_gene_clustered_transcripts(clustered_genes, clustered_genes_lock, 
                                      all_genes_and_fnames, fnames,
                                      ofp, tracking_ofp,
                                      gene_id_cntr):
    while True:
        with clustered_genes_lock:
            if len(clustered_genes) == 0: 
                return
            genes = clustered_genes.pop()
        
        with gene_id_cntr.get_lock():
            gene_id_cntr.value += 1
        
        # group transcript by their internal structure
        internal_clustered_transcript_groups = defaultdict( list )
        for sample_i, gene_id in genes:
            region, fname = all_genes_and_fnames[sample_i][gene_id]
            with open(fname) as fp:
                gene = pickle.load(fp)
            for transcript in gene.transcripts:
                IB_key = transcript.IB_key()
                internal_clustered_transcript_groups[IB_key].append(
                    (transcript, sample_i))
        
        # reduce the non single exon genes
        reduced_transcripts = []
        new_gene_id = "GENE_%i" % gene_id_cntr.value
        for internal_clustered_transcripts in \
                internal_clustered_transcript_groups.values():
            for transcript, old_transcripts, sources in \
                    reduce_internal_clustered_transcripts( 
                        internal_clustered_transcripts, new_gene_id ):
                reduced_transcripts.append(
                    (transcript, new_gene_id, old_transcripts, sources))

        gtf_lines = []
        tracking_lines = []
        for transcript, gene_id, old_ts, sources in reduced_transcripts:
            if FIX_CHRM_NAMES_FOR_UCSC:
                transcript.chrm = fix_chrm_name_for_ucsc(transcript.chrm)
            gtf_lines.append(
                transcript.build_gtf_lines(meta_data={}, source='GRIT'))
            for old_t, source_i in zip(old_ts, sources):
                tracking_lines.append(
                    "\t".join((old_t.id, fnames[source_i], transcript.id)))
        ofp.write("\n".join(gtf_lines)+"\n")
        tracking_ofp.write("\n".join(tracking_lines)+"\n")
    
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

import multiprocessing
NTHREADS = 32

def main():
    manager = multiprocessing.Manager()
    all_genes_and_fnames = manager.list()
    
    # load all the gtfs
    print >> sys.stderr, "Loading gtfs"
    manager = multiprocessing.Manager()
    all_genes_and_fnames = manager.list()
    all_genes_and_fnames_lock = multiprocessing.Lock()
    fnames = sys.argv[1:]
    pids = []
    for fname in fnames:
        pid = os.fork()
        if pid == 0:
            genes_and_fnames = load_gtf_into_pickled_files(fname)
            with all_genes_and_fnames_lock:
                all_genes_and_fnames.append(genes_and_fnames)
            os._exit(0)
        else:
            pids.append(pid)
    for pid in pids:
        os.waitpid(pid, 0)
    all_genes_and_fnames = list(all_genes_and_fnames)
    del all_genes_and_fnames_lock
    
    # group overlapping genes
    print >> sys.stderr, "Grouping genes"
    grpd_genes = manager.list()
    grpd_genes_lock = multiprocessing.Lock()
    for genes in group_overlapping_genes(all_genes_and_fnames):
        grpd_genes.append(genes)
    
    # merge gene clustered transcripts
    print >> sys.stderr, "Merging transcripts"
    ofp = ThreadSafeFile("merged.gtf", "w")
    tracking_ofp = ThreadSafeFile("tracking.txt", "w")
    gene_id_cntr = multiprocessing.Value('i', 0)

    """
    reduce_gene_clustered_transcripts(
        grpd_genes, grpd_genes_lock, 
        all_genes_and_fnames, fnames,
        ofp, tracking_ofp, gene_id_cntr)
    """
    
    pids = []
    for i in xrange(NTHREADS):
        pid = os.fork()
        if pid == 0:
            reduce_gene_clustered_transcripts(
                grpd_genes, grpd_genes_lock, 
                all_genes_and_fnames, fnames,
                ofp, tracking_ofp, gene_id_cntr)
            os._exit(0)
        else:
            pids.append(pid)
    for pid in pids:
        os.waitpid(pid, 0)
    
    
    ofp.close()
    tracking_ofp.close()

main()
