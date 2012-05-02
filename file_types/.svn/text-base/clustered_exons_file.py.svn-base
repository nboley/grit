import sys
from collections import defaultdict
import numpy

from gtf_file import parse_gff_line

def parse_clustered_exons_file( exons_fp ):
    clusters = defaultdict( set )
    for line in exons_fp:
        gff = parse_gff_line( line )
        if gff == None: continue
        
        key = ( gff.group, gff.region.chr, gff.region.strand )
        clusters[ key ].add( (gff.region.start, gff.region.stop) )
    
    sorted_clusters = defaultdict( list )
    for key, exons in clusters.iteritems():
        gene_name, chrm, strand = key
        sorted_clusters[ (chrm, strand) ].append( \
            sorted( exons ) )
    
    return sorted_clusters 

