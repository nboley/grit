#!/usr/bin/python

VERBOSE = False

import sys
import os

from collections import defaultdict

# import slide modules
sys.path.append( os.path.join( os.path.dirname(__file__), "../../", "file_types" ) )
from junctions_file import parse_junctions_file, write_junctions

def get_jns_set( gff_fp ):
    all_jns = set()
    jns = parse_junctions_file( gff_fp )
    for jn, cnt, grp in jns.iter_jns_and_cnts_and_grps():
        all_jns.add( jn )
    
    return all_jns

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description='Get exon from a variety of sources.')
    parser.add_argument( 'reference_gff', type=file, \
        help='Reference junctions.')
    parser.add_argument( 'sample_gff', type=file, \
        help='GTF file which is to be intersected with the reference gff file.'  )
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', \
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.reference_gff, args.sample_gff

def main():
    ref_jns_fp, sample_jns_fp = parse_arguments()
    ref_jns_set = get_jns_set( ref_jns_fp )    

    jns = parse_junctions_file( sample_jns_fp )
    if len( jns ) == 0:
        return
    
    # group junctions by introns
    grped_jns = defaultdict( int )
    for jn, cnt, grp in jns.iter_jns_and_cnts_and_grps():
        grped_jns[ jn ] += cnt    
    
    jns, cnts = zip( *sorted(grped_jns.iteritems()) )
    
    
    write_junctions( jns, sys.stdout, cnts )
    
    return

if __name__ == "__main__":
    main()

