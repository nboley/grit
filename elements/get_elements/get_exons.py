# Copyright (c) 2011-2012 Nathan Boley

#!/usr/bin/python

VERBOSE = False
MIN_VERBOSE = True

import sys
import os

# import slide modules
sys.path.append( os.path.join( os.path.dirname(__file__), ".." ) )
from gene_models import parse_gff_line

def create_exon_gff_line( exon ):
    gff_line = ['chr' + exon.chr, 'get_elements', 'exon', str(exon.start), \
                    str(exon.stop), '.', exon.strand, '.' ]
    return '\t'.join( gff_line )

def write_exons( exons, out_fp ):
    for exon in sorted( exons ):
        out_fp.write( create_exon_gff_line( exon ) + '\n' )
    out_fp.close()

def get_exons_from_gtf( gtf_fp ):
    exons = set()
    
    for line in gtf_fp:
        gene_name, trans_name, feature_type, region = parse_gff_line( line )
        if region != None:
            # store exon information in exons set
            exons.add( region )
    
    return exons

def get_exons( trans_gtfs ):
    # create empty set of exons
    exons = set()
    
    if trans_gtfs:
        # process transcript gtf files
        for trans_fp in trans_gtfs:
            exons = exons.union( get_exons_from_gtf( trans_fp ) )
            trans_fp.close()
    
    return exons

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description='Get exon from a variety of sources.')
    parser.add_argument( '--transcripts_gtfs', '-t', type=file, nargs='*', \
        help='GTF format file which contains exons associated with transcripts. (eg. flybase, ginormous)')
    
    parser.add_argument( '--out_fname', '-o', \
        help='Output file will be written to default. default: stdout')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', \
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    out_fp = open( args.out_fname, "w" ) if args.out_fname else sys.stdout
    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.transcripts_gtfs, out_fp

def main():
    trans_gtfs, out_fp = parse_arguments()
    
    if MIN_VERBOSE:
        print 'Getting exons...'
    exons = get_exons( trans_gtfs )
    
    # write exons out to an exons file
    write_exons( exons, out_fp )
    
    return

if __name__ == "__main__":
    main()

