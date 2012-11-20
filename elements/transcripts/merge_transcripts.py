# Copyright (c) 2011-2012 Nathan Boley

import sys, os

import numpy
from collections import defaultdict, OrderedDict
from itertools import izip, repeat
import networkx as nx

sys.path.append( os.path.join( os.path.dirname( __file__ ),
                               "..", "..", "file_types", "fast_gtf_parser" ) )
from gtf import load_gtfs

sys.path.append( os.path.join( os.path.dirname( __file__ ), ".."  ) )
from cluster_exons import find_overlapping_exons


VERBOSE = True
MERGE_DISTAL_ENDS = False
SINGLE_LINKAGE_CLUSTER = True
LINKAGE_CLUSTER_GAP = 50
LINKAGE_CLUSTER_MAX_DIFF = 200

def build_gtf_lines( chrm, strand, gene_name, trans_name, exon_bndries ):
    feature="exon"
    source="merge_transcripts"
    gene_name_iter = repeat( gene_name )
    trans_name_iter = repeat( trans_name )
    
    # create an iterator over every exon
    exons_iter = enumerate( izip( exon_bndries[0::2], exon_bndries[1::2] ) )
    regions_iter = ( GenomicInterval( chrm, strand, start, stop )
                     for exon_num, (start, stop) in exons_iter )
    
    return iter_gtf_lines( regions_iter, gene_name_iter, trans_name_iter,
                           feature=feature, source=source )

def write_all_transcripts( clustered_trans, ofp, sources_fp ):
    # iterate through unique transcripts, and write them to the output file
    for (chrm, strand), genes in clustered_trans.iteritems():
        for gene_id, transcripts in genes.iteritems():
            gene_name = "merged_cluster_%s" % gene_id
            for trans_id, exons in enumerate( transcripts ):
                trans_name = gene_name + "_%i" % trans_id
                
                if sources_fp != None:
                    sources = [ os.path.basename(source) for source in
                                transcript_maps[ ( chrm, strand ) ][ exons ] ]
                    source_str = ",".join( sources )
                    sources_fp.write( '\t'.join(
                            (gene_name, trans_name, source_str ) ) + '\n' )
                
                gtf_lines = build_gtf_lines(
                    chrm, strand, gene_name, trans_name, exons )
                
                ofp.write( "\n".join( gtf_lines ) + "\n" )
    
    ofp.close()
    if sources_fp != None:
        sources_fp.close()
    
    return


def cluster_ends( ends_and_sources ):
    def create_cluster_map( coords, map_to_index ):
        coords = sorted( coords )
        
        coord_map = {}
        coord_range = [coords[0], coords[0]]
        cluster_coords = set()
        
        for coord in coords:
            if coord < coord_range[0] - LINKAGE_CLUSTER_GAP or \
                    coord > coord_range[1] + LINKAGE_CLUSTER_GAP:    
                for i_coord in cluster_coords:
                    coord_map[i_coord] = coord_range[ map_to_index ]
                
                coord_range = [coord, coord]
                cluster_coords = set( (coord,) )
            else:
                coord_range[0] = min( coord_range[0], coord )
                coord_range[1] = max( coord_range[1], coord )
                cluster_coords.add( coord )
        
        for i_coord in cluster_coords:
            coord_map[i_coord] = coord_range[ map_to_index ]
        
        return coord_map
    
    starts_map = create_cluster_map( zip( *ends_and_sources )[0], 0 )
    stops_map = create_cluster_map( zip( *ends_and_sources )[1], 1 )
    
    clustered_ends_and_sources = []
    clusters_map = defaultdict( set )
    for start, stop, source in ends_and_sources:
        if start > starts_map[start] + LINKAGE_CLUSTER_MAX_DIFF or \
                stop < stops_map[stop] - LINKAGE_CLUSTER_MAX_DIFF:
            # should really add these to a new list and recluster 
            # to avoid wierd regions near long linkage chains
            clustered_ends_and_sources.append( (start, stop, source) )
        else:
            clusters_map[ (starts_map[start], stops_map[stop] ) ].add( source )
    
    for (start, stop), sources in clusters_map.iteritems():
        clustered_ends_and_sources.append( (start, stop, sources) )
    
    return clustered_ends_and_sources

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

def build_exon_to_transcripts_map( transcripts, source_fnames ):
    exon_to_transcripts_map = defaultdict( lambda: defaultdict(list) )
    for source, source_transcripts in izip( source_fnames, transcripts ):
        for cluster_id, chrm, strand, start, stop, transcripts in source_transcripts:
            for transcript in transcripts:
                for exon in transcript.exons:
                    exon_to_transcripts_map[(chrm, strand)][exon].append( 
                        (transcript, source) )
    
    return exon_to_transcripts_map

def cluster_transcripts( transcripts, sources ):
    """Two transcripts overlap if they share overlapping exons So, we  
    use the following algorithm:
    
    1) build a mapping from exons to transcripts
    2) build a mapping from transcripts to exons
    3) cluster exons
    4) build edges through exons

    """
    exons_to_transcripts = defaultdict( lambda: defaultdict(list) )
    transcripts_to_exons = defaultdict( lambda: defaultdict(list) )
    for source_genes, source in izip( transcripts, sources ):
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
        
        sorted_exons = sorted(exons.keys())
        clustered_exons = find_overlapping_exons( sorted_exons )
        for exon_cluster in clustered_exons:
            connected_transcripts = []
            for exon_i in exon_cluster:
                connected_transcripts.extend( 
                    exons_to_transcripts[contig][sorted_exons[exon_i]] )
            
            transcripts_graph.add_edges_from( 
                zip(connected_transcripts[:-1], connected_transcripts[1:] ) )
            
            
        clustered_transcripts[contig] = nx.connected_components(
            transcripts_graph )
    
    return clustered_transcripts

def reduce_internal_clustered_transcripts( internal_grpd_transcripts ):
    """Take a set of clustered transcripts and reduce them into 
    a set of canonical transcripts, and associated sources.
    
    """
    reduced_transcripts = defaultdict( list )
    for grpd_transcripts in internal_grpd_transcripts:
        for transcript, source in grpd_transcripts:
            reduced_transcripts[transcript].append( source )
    
    return sorted( reduced_transcripts.iteritems() )

def reduce_single_exon_transcripts(se_transcripts):
    reduced_transcripts = defaultdict( list )
    for transcript, source in se_transcripts:
        transcripts[transcript].append( source )
    
    return sorted( reduced_transcripts.iteritems() )

def reduce_clustered_transcripts(clustered_transcripts):
    reduced_transcripts = []
    gene_id_cntr = 1
    for (chr, strand), contig_transcripts in clustered_transcripts.iteritems():
        for gene_transcripts in contig_transcripts:
            gene_id = "GENE_%i" % gene_id_cntr
            
            # group transcript by their internal structure
            internal_clustered_transcripts = defaultdict( list )
            single_exon_transcripts = []
            for transcript, source in gene_transcripts:
                IB_key = transcript.IB_key(error_on_SE_genes=False)
                if IB_key[0] == 'SE_GENE':
                    single_exon_transcripts.append( (transcript, source) )
                else:
                    internal_clustered_transcripts[IB_key].append(
                        (transcript, source))

            # reduce the non single exon genes
            for transcript, sources in reduce_internal_clustered_transcripts( 
                    internal_clustered_transcripts.values() ):
                reduced_transcripts.append( (transcript, gene_id, sources) )

            # reduce the single exon genes
            for transcript, sources in reduce_single_exon_transcripts( 
                    single_exon_transcripts ):
                reduced_transcripts.append( (transcript, gene_id, sources) )
            
            gene_id_cntr += 1
    
    return reduced_transcripts

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
        sources_fp.write("\t".join((gene_id, transcript.id, ",".join(sources)))+"\n")
    
    return 

if __name__ == "__main__":
    main()
