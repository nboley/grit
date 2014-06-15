import os, sys
import cPickle as pickle
import multiprocessing
import copy
import random

from itertools import izip, chain
from collections import defaultdict

import numpy
from scipy.cluster.hierarchy import fclusterdata

from files.gtf import load_gtf_into_pickled_files
from files.reads import fix_chrm_name_for_ucsc
from transcript import Transcript, Gene
from lib.multiprocessing_utils import ThreadSafeFile
import config

LINKAGE_CLUSTER_GAP = 500

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

def reduce_internal_clustered_transcripts( 
        internal_grpd_transcripts, gene_id, max_cluster_gap ):
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
    cluster_indices = fclusterdata( transcript_ends, t=max_cluster_gap,
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

def reduce_gene_clustered_transcripts(
        genes, new_gene_id,
        min_upper_fpkm, 
        max_intrasample_fpkm_ratio, 
        max_intersample_fpkm_ratio,
        max_cluster_gap):
    # unpickle the genes, and find the expression thresholds    
    unpickled_genes = []
    max_fpkm_lb_across_samples = 0
    max_fpkm_lb_in_sample = defaultdict(float)
    
    for gtf_fname, pickled_fname in genes:
        with open(pickled_fname) as fp:
            gene = pickle.load(fp)
        unpickled_genes.append((gtf_fname, gene))
        try: 
            max_fpkm_lb_in_gene = max( 
                transcript.conf_lo for transcript in gene.transcripts
                if transcript.conf_lo != None )
        except ValueError:
                continue
        max_fpkm_lb_across_samples = max( 
                max_fpkm_lb_across_samples, max_fpkm_lb_in_gene )
        max_fpkm_lb_in_sample[gtf_fname] = max(
                max_fpkm_lb_in_gene, max_fpkm_lb_in_sample[gtf_fname])
    genes = unpickled_genes
    del unpickled_genes
    
    # group transcript by their internal structure
    internal_clustered_transcript_groups = defaultdict( list )
    for gtf_fname, gene in genes:
        min_max_fpkm = max(
                min_upper_fpkm, 
                max_fpkm_lb_in_sample[gtf_fname]/max_intrasample_fpkm_ratio,
                max_fpkm_lb_across_samples/max_intersample_fpkm_ratio)
        
        for transcript in gene.transcripts:
            if transcript.conf_hi < min_max_fpkm: continue
            IB_key = transcript.IB_key()
            internal_clustered_transcript_groups[IB_key].append(
                (transcript, gtf_fname))
    
    # reduce the non single exon genes
    merged_transcripts = []
    merged_transcript_sources = []
    transcript_id = 1
    for internal_clustered_transcripts in \
            internal_clustered_transcript_groups.values():
        for (merged_transcript, old_transcripts, sources
                ) in reduce_internal_clustered_transcripts( 
                internal_clustered_transcripts, new_gene_id, max_cluster_gap ):
            merged_transcript.id = new_gene_id + "_%i" % transcript_id
            transcript_id += 1
            merged_transcripts.append(merged_transcript)
            merged_transcript_sources.append(
                zip(old_transcripts, sources))

    new_gene = Gene(new_gene_id, None,
                    chrm=merged_transcripts[0].chrm,
                    strand=merged_transcripts[0].strand,
                    start=merged_transcripts[0].start,
                    stop=merged_transcripts[0].stop,
                    transcripts=merged_transcripts )
    
    return new_gene, merged_transcript_sources

def group_overlapping_genes(all_sources_and_pickled_gene_fnames):
    chrm_grpd_genes = defaultdict(list)
    for gtf_fname, pickled_gene_fnames in all_sources_and_pickled_gene_fnames:
        for pickled_gene_fname in pickled_gene_fnames:
            with open(pickled_gene_fname) as fp:
                gene = pickle.load(fp)
            chrm_grpd_genes[(gene.chrm, gene.strand)].append(
                (gene.start, gene.stop, pickled_gene_fname))
    
    grpd_genes = []
    for (chrm, strand), gene_regions in chrm_grpd_genes.iteritems():
        gene_regions.sort()
        curr_stop = -1
        for start, stop, pickled_gene_fname in gene_regions:
            if start > curr_stop: 
                grpd_genes.append([])
            grpd_genes[-1].append((gtf_fname, pickled_gene_fname))
            curr_stop = stop
    
    return grpd_genes
