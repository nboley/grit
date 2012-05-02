#!/usr/bin/python

VERBOSE = True
MIN_VERBOSE = True

import sys
import os

# import slide modules
sys.path.append( os.path.join( os.path.dirname(__file__), ".." ) )
import slide
slide.VERBOSE = VERBOSE
from slide import get_raw_transcripts
from gene_models import parse_gff_line, GenomicInterval

def create_junction_gff_line( junct ):
    gff_line = [ junct.chr, 'elements', 'junction', str(junct.start), \
                     str(junct.stop ), '.', junct.strand, '.' ]
    return '\t'.join( gff_line )

def write_junctions( junctions, out_fp ):
    for junct in sorted( junctions ):
        out_fp.write( create_junction_gff_line( junct ) + '\n' )
    out_fp.close()
            
def get_juncts_from_trans( trans_fp ):
    junctions = set()

    # get the gene -> trans ->  exons structure from the gtf file
    # verification of transcripts is done from within get_raw_transcripts
    raw_gene_transcripts = get_raw_transcripts( trans_fp )
    for gene_name, transcripts in raw_gene_transcripts.iteritems():
        for trans_name, exons in transcripts.iteritems():
            for i, exon in enumerate( exons[:-1] ):
                junction = GenomicInterval( exon.chr, exon.strand, exon.stop + 1, \
                                                exons[i+1].start - 1 )
                junctions.add( junction )

    return junctions

def get_juncts_from_junct( junct_fp ):
    junctions = set()
    for line in junct_fp:
        gene_name, trans_name, feature_type, region = parse_gff_line( line )
        junctions.add( region )

    return junctions

def get_junctions( trans_gtfs, junct_gffs ):
    # create empty set of exons
    junctions = set()

    if trans_gtfs:
        # process transcript gtf files
        for trans_fp in trans_gtfs:
            junctions = junctions.union( get_juncts_from_trans( trans_fp ) )
            trans_fp.close()

    if junct_gffs:
        # process junctions gff files
        for junct_fp in junct_gffs:
            junctions = junctions.union( get_juncts_from_junct( junct_fp ) )
            junct_fp.close()
    
    return junctions

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description='Get exon from a variety of sources.')
    parser.add_argument( '--transcripts_gtfs', '-t', type=file, nargs='*', \
        help='GTF file which contains exons associated with transcripts. Introns will be inferred from exon structure. (eg. flybase, ginormous)')
    parser.add_argument( '--junctions_gffs', '-j', type=file, nargs='*', \
        help='GFF file which contains intron structure. (eg. brenton\'s, all_junts)')
    parser.add_argument( '--out_fname', '-o', \
        help='Output file will be written to default. default: stdout')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', \
        help='Whether or not to print status information.')
    args = parser.parse_args()

    out_fp = open( args.out_fname, "w" ) if args.out_fname else sys.stdout

    global VERBOSE
    VERBOSE = args.verbose

    return args.transcripts_gtfs, args.junctions_gffs, out_fp

if __name__ == "__main__":
    trans_gtfs, juncts_gffs, out_fp = parse_arguments()

    if MIN_VERBOSE:
        print 'Getting junctions...'
    junctions = get_junctions( trans_gtfs, juncts_gffs )
    
    # write exons out to an exons file
    write_junctions( junctions, out_fp )
