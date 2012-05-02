from collections import defaultdict
import numpy

from gtf_file import parse_gff_line

def _load_regions( exons_fp ):
    all_exons = defaultdict( set )
    for line in exons_fp:
        gff = parse_gff_line( line )
        if gff == None: continue
        all_exons[(gff.region.chr, gff.region.strand)].add( \
            (gff.region.start, gff.region.stop) )
    
    return all_exons

def _convert_into_sorted_numpy_array( all_exons ):
    sorted_exons = {}
    for key, exons in all_exons.iteritems():
        sorted_exons[key] = numpy.array( sorted( exons ) )
    
    return sorted_exons

def parse_exons_file( exons_fp, merge_5p_bndry=True, merge_3p_bndry=True ):
    all_exons = _load_regions( exons_fp )
    return _convert_into_sorted_numpy_array( all_exons )

def parse_exons_files( exons_fps ):
    all_exons = defaultdict( set )
    for exons_fp in exons_fps:
        inner_exons = _load_regions( exons_fp )
        for key, exons in inner_exons.iteritems():
            all_exons[ key ].update( exons )
    
    return _convert_into_sorted_numpy_array( all_exons )
