# Copyright (c) 2011-2012 Nathan Boley

#!/usr/bin/python

VERBOSE = False

import sys
import os

from collections import defaultdict

# import slide modules
sys.path.append( os.path.join( os.path.dirname(__file__), "../../", "file_types" ) )
from junctions_file import parse_jn_gff

def get_jns_set( fname ):
    all_jns = set( jn.region for jn in  parse_jn_gff( fname ) )

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
    
    ref_jns_set = set( jn.region for jn in parse_jn_gff( ref_jns_fp ) )
    sample_jns = parse_jn_gff( sample_jns_fp )
    jn_id = 0
    for jn in sample_jns:
        if jn.region in ref_jns_set:
            jn_id += 1
            print jn.build_gff_line( jn_id )
    
    """
    # group junctions by introns
    grped_jns = defaultdict( int )
    for jn, cnt, grp in jns.iter_jns_and_cnts_and_grps():
        grped_jns[ jn ] += cnt    
    """
    
    return

if __name__ == "__main__":
    main()

