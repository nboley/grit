from collections import defaultdict
from itertools import product

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
