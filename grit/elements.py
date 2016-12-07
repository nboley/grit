"""
Copyright (c) 2011-2015 Nathan Boley

This file is part of GRIT.

GRIT is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GRIT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GRIT.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy
import networkx as nx

from collections import defaultdict, namedtuple
from itertools import chain

RefElementsToInclude = namedtuple(
    'RefElementsToInclude', 
    ['genes', 'junctions', 'TSS', 'TES', 'promoters', 'polya_sites', 'exons'])


def convert_elements_to_arrays(all_elements):
    # convert into array
    all_array_elements = defaultdict( 
        lambda: defaultdict(lambda: numpy.zeros(0)) )
    for key, elements in all_elements.items():
        for element_type, contig_elements in elements.items():
            all_array_elements[key][element_type] \
                = numpy.array( sorted( contig_elements ) )

    return all_array_elements

def load_elements( fp ):
    all_elements = defaultdict( lambda: defaultdict(set) )
    for line in fp:
        if line.startswith( 'track' ): continue
        chrm, start, stop, element_type, score, strand = line.split()[:6]
        # subtract 1 from stop becausee beds are closed open, and we 
        # wnat everything in 0-based closed-closed
        all_elements[(chrm, strand)][element_type].add( 
            (int(start), int(stop)-1) )
    
    return convert_elements_to_arrays(all_elements)

def find_overlapping_exons(exons):
    overlapping_exons_mapping = set()
    for o_i, (o_start, o_stop) in enumerate(exons):
        for i_i, (i_start, i_stop) in enumerate(exons):
            if i_i > o_i: break
            if not (i_stop < o_start or i_start > o_stop):
                overlapping_exons_mapping.add( 
                    ((o_start, o_stop), (i_start, i_stop)) )
                overlapping_exons_mapping.add( 
                    ((i_start, i_stop), (o_start, o_stop)) )
    
    return list(overlapping_exons_mapping)


def find_jn_connected_exons(exons, jns, strand):
    edges = set()
    
    # build mappings from exon starts to indices
    exon_starts_map = defaultdict(list)
    exon_stops_map = defaultdict(list)
    for start, stop in exons:
        exon_starts_map[start].append( (start, stop) )
        exon_stops_map[stop].append( (start, stop ) )
    
    for jn in jns:
        for start_exon in exon_stops_map[jn[0]-1]:
            for stop_exon in exon_starts_map[jn[1]+1]:
                if strand == '+':
                    edges.add((tuple(jn), start_exon, stop_exon))
                else:
                    edges.add((tuple(jn), stop_exon, start_exon))
    
    return edges

def iter_nonoverlapping_exons(exons):
    if len(exons) == 0: return
    try: exons = exons.tolist()
    except AttributeError: pass
    exons = [tuple(x) for x in exons]
    G = nx.Graph()
    G.add_nodes_from(exons)
    overlapping_exons = find_overlapping_exons(exons)
    G.add_edges_from( overlapping_exons )
    
    for clustered_exons in nx.connected_components(G):
        if len( clustered_exons ) == 1:
            yield clustered_exons[0]
    
    return

def cluster_elements( tss_exons, internal_exons, tes_exons, se_transcripts, 
                      promoters, polyas, jns, strand ):
    assert isinstance( tss_exons, set )
    assert isinstance( internal_exons, set )
    assert isinstance( tes_exons, set )
    assert isinstance( se_transcripts, set )
    assert isinstance( promoters, set )
    assert isinstance( polyas, set )
    
    all_exons = sorted( chain(tss_exons, internal_exons, 
                              tes_exons, se_transcripts,
                              promoters, polyas) )
    if len(all_exons) == 0: return
    
    G = nx.Graph()
    G.add_nodes_from(all_exons)
    overlapping_exons = find_overlapping_exons(all_exons)
    G.add_edges_from( overlapping_exons )
    jns_and_connected_exons = find_jn_connected_exons(all_exons, jns, strand)
    edges = []
    observed_jns = set()
    for jn, start, stop in jns_and_connected_exons:
        observed_jns.add(jn)
        edges.append((start, stop))
    G.add_edges_from(edges)
    
    for gene in nx.connected_component_subgraphs(G):
        exons = gene.nodes()
        yield ( tss_exons.intersection( exons ),
                tes_exons.intersection( exons ),
                internal_exons.intersection( exons ),
                se_transcripts.intersection( exons ),
                promoters.intersection( promoters ),
                polyas.intersection( polyas ),
                sorted(observed_jns) 
              )
    
    return
