#!/usr/bin/python

VERBOSE = False

import os
import sys

from pysam import Samfile, Fastafile
from collections import defaultdict

sys.path.append( os.path.join( os.path.dirname(__file__), \
                                   "..", "..", "file_types" ) )
from junctions_file import Junction, GenomicInterval, get_jn_type

sys.path.append( os.path.join( os.path.dirname(__file__), ".." ) )
from build_read_coverage_bedgraph import clean_chr_name

def read_spans_single_intron( read ):
    # quickly check if read could spans a single intron
    if len( read.cigar ) != 3:
        return False
    
    # check that cigar regions are alternating matching(exon) and \
    # skipping(intron) regions
    for i, cigar_region in enumerate(read.cigar):
        if (i % 2 == 0 and cigar_region[0] != 0) or \
                (i % 2 == 1 and cigar_region[0] != 3):
            return False
    
    # if the read is not primary then do not include it in the counts
    if read.is_secondary:
        return False
    
    return True

def iter_junctions( reads_fn, fasta_fn, stranded, reverse_strand ):
    """Get all of the junctions represented in a reads object
    """
    def find_junctions_in_chrm( chrm, all_junctions, unique_junctions ):
        # store how many times we have observed each read
        query_name_cnts = defaultdict( int )
        jn_reads = []
        for read in reads.fetch(chrm):
            # increment the number of times we've seen this read
            query_name_cnts[ read.qname ] += 1
            if read_spans_single_intron( read ):
                jn_reads.append( read )

        for read in jn_reads:
            # obtain chrm from read pysam object
            chrm = clean_chr_name( reads.getrname( read.tid ) )

            # add one to left_intron since bam files are 0-based
            upstrm_intron_pos = read.pos + read.cigar[0][1] + 1
            dnstrm_intron_pos = upstrm_intron_pos + read.cigar[1][1] - 1

            # if the protocol is stranded, get the strand from the read mapping 
            # directions
            if stranded:
                if (read.is_read1 and not read.is_reverse) \
                        or (read.is_read2 and read.is_reverse):
                    strand = '+'
                else:
                    strand = '-'

                # if the protocol strand is reversed, then reverse the strand
                if reverse_strand:
                    if strand == '+': 
                        strand = '-'
                    else: 
                        strand = '+'
            else:
                jn_type, strand = get_jn_type( \
                    chrm, upstrm_intron_pos, dnstrm_intron_pos, fasta_obj )
                if jn_type == 'non-canonical': 
                    continue

            # increment count of junction reads at this read position
            # for this intron or initialize it to 1
            jn_type = ( chrm, strand, upstrm_intron_pos, \
                            dnstrm_intron_pos, read.cigar[0][1] )

            all_junctions[ jn_type ] += 1

            if query_name_cnts[ read.qname ] == 2:
                unique_junctions[ jn_type ] += 1
    
    # build reads object
    reads = Samfile( reads_fn, "rb" )
    fasta_obj = None if fasta_fn == None else Fastafile( fasta_fn )
    
    # store the number of reads across the junction for each relative 
    # read position
    all_junctions = defaultdict(int)
    unique_junctions = defaultdict(int)
    
    for ref in reads.references:
        find_junctions_in_chrm( ref, all_junctions, unique_junctions )
    
    jns = []
    read_offsets = []
    cnts = []
    uniq_cnts = []
    for (chrm, strand, start, stop, read_offset ), cnt \
            in sorted(all_junctions.iteritems()):

        uniq_cnt = unique_junctions[ (chrm, strand, start, stop, read_offset ) ]
        jn = Junction( GenomicInterval(chrm, strand, start, stop), 
                       source_read_offset=read_offset, 
                       cnt=cnt, uniq_cnt=uniq_cnt )
        yield jn
        
    return


def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(description=\
       'Parse junctions from provided bam file and write to std out.')
    parser.add_argument( 'bam_fn', help="BAM file to parse jns from." )
    parser.add_argument( '--fasta', '-f', default=None,\
       help='Fasta file used to determine the consensus intron type' \
           + '(fasta file expects to be indexed: use samtools faidx).')    
    parser.add_argument( '--verbose', '-v', default=False, \
                             action='store_true', \
       help='Whether or not to print status information.' )
    parser.add_argument( '--reverse-strand', '-r', default=False, \
                             action='store_true', \
       help='Whether or not to reverse the junction strand.' )
    parser.add_argument( '--stranded', '-s', default=False, \
                             action='store_true', \
       help='Needs to be set if the jn file type is unstranded.' )
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose

    if not args.stranded:
        assert not args.reverse_strand
        assert args.fasta != None
    
    return args.bam_fn, args.fasta, args.stranded, args.reverse_strand

def main():
    reads_fn, fasta_fn, stranded, reverse_strand = parse_arguments()
    
    # get junctions
    for jn in iter_junctions(reads_fn, fasta_fn, stranded, reverse_strand):
        print jn.build_gff_line()

if __name__ == "__main__":
    main()
