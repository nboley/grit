import sys

from itertools import chain
from collections import namedtuple, defaultdict
CategorizedExons = namedtuple( "CategorizedExons", ["TSS", "TES", "internal", 
                                                    "full_transcript"] )

import igraph

VERBOSE = False

def find_paths(g, n, m, path=[]):
    "Find paths from node index n to m using adjacency list a."
    path = path + [n]
    if n == m:
         return [path]
    
    paths = []
    for child in g.successors( n ):
      if child not in path:
        child_paths = find_paths(g, child, m, path)
        for child_path in child_paths:
          paths.append(child_path)
    
    return paths

def find_overlapping_exons(exons):
    overlapping_exons_mapping = set()
    for o_i, (o_start, o_stop) in enumerate(exons):
        for i_i, (i_start, i_stop) in enumerate(exons):
            if i_i > o_i: break
            if not (i_stop < o_start or i_start > o_stop):
                overlapping_exons_mapping.add( (min(i_i, o_i), max(i_i, o_i)) )
    
    return overlapping_exons_mapping


def find_jn_connected_exons(exons, jns, strand, use_names=False):
    edges = set()
    
    # build mappings from exon starts to indices
    exon_starts_map = defaultdict(list)
    exon_stops_map = defaultdict(list)
    for i, (start, stop) in enumerate(exons):
        exon_starts_map[start].append( i )
        exon_stops_map[stop].append( i )
    
    for jn in jns:
        for start in exon_stops_map[jn[0]-1]:
            for stop in exon_starts_map[jn[1]+1]:
                if use_names:
                  if strand == '+':
                    edges.add((exons[start], exons[stop]))
                  else:
                    edges.add((exons[stop], exons[start]))
                else:
                  if strand == '+':
                    edges.add((start, stop))
                  else:
                    edges.add((stop, start))
    
    return edges

def iter_nonoverlapping_exons(exons):
    genes_graph = igraph.Graph()
    genes_graph.add_vertices(len(exons)-1)
    genes_graph.add_edges( find_overlapping_exons(exons) )
    genes = genes_graph.components()
    for gene in genes:
        clustered_exons = [ exons[exon_i] for exon_i in gene ]
        if len( clustered_exons ) == 1:
            yield clustered_exons[0]
    
    return

def cluster_exons( tss_exons, internal_exons, tes_exons, se_transcripts, 
                   jns, strand ):
    assert isinstance( tss_exons, set )
    assert isinstance( internal_exons, set )
    assert isinstance( tes_exons, set )
    assert isinstance( se_transcripts, set )
    
    all_exons = sorted( chain(tss_exons, internal_exons, 
                              tes_exons, se_transcripts) )
    
    genes_graph = igraph.Graph()
    genes_graph.add_vertices( len(all_exons)-1 )
    genes_graph.add_edges( find_overlapping_exons(all_exons) )
    genes_graph.add_edges( find_jn_connected_exons(all_exons, jns, strand ) )
    
    genes = genes_graph.components()
    for gene in genes:
        exons = [ all_exons[exon_i] for exon_i in gene ]
        yield ( tss_exons.intersection( exons ),
                tes_exons.intersection( exons ),
                internal_exons.intersection( exons ),
                se_transcripts.intersection( exons ) )
    
    return

def build_transcripts( tss_exons, internal_exons, tes_exons, se_transcripts, 
                       jns, strand, max_num_transcripts=100 ):
    import networkx as nx
    # build a directed graph, with edges leading from exon to exon via junctions
    all_exons = sorted(chain(tss_exons, internal_exons, tes_exons))
    graph = nx.DiGraph()
    graph.add_nodes_from( tss_exons )
    graph.add_nodes_from( internal_exons )
    graph.add_nodes_from( tes_exons )
    
    edges = find_jn_connected_exons(all_exons, jns, strand, use_names=True )
    graph.add_edges_from( edges  )
    
    transcripts = []
    num_transcripts = 0
    for tss in tss_exons:
        for tes in tes_exons:
            for transcript in nx.all_simple_paths(graph, tss, tes):
                num_transcripts += 1
                if num_transcripts > max_num_transcripts: 
                    if VERBOSE: print >> sys.stderr, "TOO COMPLEX"
                    return []

                transcripts.append( sorted(transcript) )
    
    return transcripts + list(se_transcripts)
