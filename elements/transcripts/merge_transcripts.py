# Copyright (c) 2011-2012 Nathan Boley

import sys, os

import numpy
from collections import defaultdict, OrderedDict
from itertools import izip, repeat
import networkx as nx

sys.path.append( os.path.join( os.path.dirname( __file__ ),
                               "..", "..", "file_types", "fast_gtf_parser" ) )
from gtf import load_gtfs, Transcript

sys.path.append( os.path.join( os.path.dirname( __file__ ), ".."  ) )
from cluster_exons import find_overlapping_exons

from scipy.cluster.hierarchy import fclusterdata

VERBOSE = True
MERGE_DISTAL_ENDS = False
SINGLE_LINKAGE_CLUSTER = True
LINKAGE_CLUSTER_GAP = 50
LINKAGE_CLUSTER_MAX_DIFF = 200

def build_exon_to_transcripts_map( transcripts, source_fnames ):
    exon_to_transcripts_map = defaultdict( lambda: defaultdict(list) )
    for source, source_transcripts in izip( source_fnames, transcripts ):
        for cluster_id, chrm, strand, start, stop, transcripts \
                in source_transcripts:
            for transcript in transcripts:
                for exon in transcript.exons:
                    exon_to_transcripts_map[(chrm, strand)][exon].append( 
                        (transcript, source) )
    
    return exon_to_transcripts_map

def cluster_transcripts( genes, sources ):
    """Two transcripts overlap if they share overlapping exons So, we  
    use the following algorithm:
    
    1) build a mapping from exons to transcripts
    2) build a mapping from transcripts to exons
    3) cluster exons
    4) build edges through exons

    """
    exons_to_transcripts = defaultdict( lambda: defaultdict(list) )
    transcripts_to_exons = defaultdict( lambda: defaultdict(list) )
    for source_genes, source in izip( genes, sources ):
        for gene in source_genes:
            for transcript in gene[-1]:
                contig = ( transcript.chrm, transcript.strand )
                for exon in transcript.exons:
                    exons_to_transcripts[contig][exon].append((transcript, source))
                    transcripts_to_exons[contig][(transcript, source)].append(exon)
    
    clustered_transcripts = {}
    for contig, exons in exons_to_transcripts.iteritems():
        transcripts_graph = nx.Graph()
        # add all of the transcripts
        transcripts_graph.add_nodes_from( 
            transcripts_to_exons[contig].iterkeys() )
        
        # cluster overlapping exons
        sorted_exons = sorted(exons.keys())
        clustered_exons = find_overlapping_exons( sorted_exons )
        for exon_cluster in clustered_exons:
            connected_transcripts = []
            for exon_i in exon_cluster:
                connected_transcripts.extend( 
                    exons_to_transcripts[contig][sorted_exons[exon_i]] )
            
            # add edges between transcripts that contain overlapping exons
            transcripts_graph.add_edges_from( 
                izip(connected_transcripts[:-1], connected_transcripts[1:] ) )

        # we also need to add exons that aren't clustered, because we could
        # have two identical transcripts that share the same exons exactly
        for exon in sorted_exons:
            connected_transcripts.extend( 
                exons_to_transcripts[contig][exon] )
        
        transcripts_graph.add_edges_from( 
            izip(connected_transcripts[:-1], connected_transcripts[1:] ) )
            
        
        # finally, find all of the overlapping components
        clustered_transcripts[contig] = nx.connected_components(
            transcripts_graph )
    
    return clustered_transcripts

def reduce_internal_clustered_transcripts( internal_grpd_transcripts ):
    """Take a set of clustered transcripts and reduce them into 
    a set of canonical transcripts, and associated sources.
    
    """
    # if there is only a single trnascript, clustering doesnt make sense
    if len( internal_grpd_transcripts ) == 1: 
        yield internal_grpd_transcripts[0][0],[internal_grpd_transcripts[0][1],]
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
    for cluster_index, clustered_transcripts in \
            clustered_transcript_grps.iteritems():
        start, stop = 1e20, 0
        for transcript in clustered_transcripts:
            start = min( start, transcript.exons[0][0] )
            stop = max( stop, transcript.exons[-1][-1] )
        
        # choose a tempalte transcript, and make sure that all of the 
        # clustered transcripts have the same internal structure (
        # this should be guaranteed by the calling function )
        bt = clustered_transcripts[0]
        assert all( t.IB_key() == bt.IB_key() for t in clustered_transcripts )
        new_trans_id = "_".join( t.id for t in clustered_transcripts )
        new_exon_bnds = bt.exon_bnds
        new_exon_bnds[0] = start
        new_exon_bnds[-1] = stop
        yield ( Transcript( new_trans_id, bt.chrm, bt.strand, 
                            new_exon_bnds, bt.cds_region ), 
                clustered_transcript_grp_sources[cluster_index] )
    
    return

def reduce_clustered_transcripts(clustered_transcripts):
    reduced_transcripts = []
    gene_id_cntr = 1
    for (chr, strand), contig_transcripts in clustered_transcripts.iteritems():
        for gene_transcripts in contig_transcripts:
            gene_id = "GENE_%i" % gene_id_cntr
            
            # group transcript by their internal structure
            internal_clustered_transcript_groups = defaultdict( list )
            for transcript, source in gene_transcripts:
                IB_key = transcript.IB_key()
                internal_clustered_transcript_groups[IB_key].append(
                    (transcript, source))

            # reduce the non single exon genes
            for internal_clustered_transcripts in \
                    internal_clustered_transcript_groups.values():
                for transcript, sources in reduce_internal_clustered_transcripts( 
                        internal_clustered_transcripts ):
                    reduced_transcripts.append( (transcript, gene_id, sources) )
            
            gene_id_cntr += 1
    
    return reduced_transcripts

def parse_arguments():
    import argparse
    desc = 'Merge transcripts.'
    parser = argparse.ArgumentParser(description=desc)    
    parser.add_argument(
        "gtfs", type=str, nargs="+",
        help='GTF file contraining transcripts. (i.e. from build_transcripts)' )
    
    parser.add_argument(
        '--sources-fname', type=argparse.FileType( 'w' ),
        help='File name to write the sources for each transcript. '
        + 'Default: Do not write out sources map.')
    
    parser.add_argument(
        '--out-fname', '-o', type=argparse.FileType('w'), default=sys.stdout,
        help='Output file. default: stdout')
    
    parser.add_argument(
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
 
    parser.add_argument(
        '--threads', '-t', type=int, default=1,
        help='The number of threads to use.')
   
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.gtfs, args.out_fname, args.sources_fname, \
        args.threads

def main():
    gtf_fnames, ofp, sources_fp, n_threads = parse_arguments()
    
    # load all transcripts file and store corresponding data
    transcriptomes = load_gtfs( gtf_fnames, n_threads )
    
    # cluster transcripts by shared overlapping exons
    clustered_transcripts = cluster_transcripts( transcriptomes, gtf_fnames )
    
    # group transcripts by matching internal boundaries
    reduced_transcripts = reduce_clustered_transcripts( clustered_transcripts )
    
    transcript_ids = defaultdict( lambda: 1 )
    for transcript, gene_id, sources in reduced_transcripts:
        transcript.id = "%s_%i" % ( gene_id, transcript_ids[gene_id] )
        transcript_ids[gene_id] += 1        
        ofp.write(transcript.build_gtf_lines(gene_id, {}, source="grit")+"\n")
        line = "\t".join((gene_id, transcript.id, ",".join(sources)))
        sources_fp.write(line+"\n")
    
    return 

if __name__ == "__main__":
    main()
