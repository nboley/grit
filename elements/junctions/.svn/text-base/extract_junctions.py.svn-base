#!/usr/bin/python

VERBOSE = False

import os
import sys

from pysam import Samfile
from collections import defaultdict
from operator import itemgetter

sys.path.append( os.path.join( os.path.dirname(__file__), \
                                   "..", "..", "file_types" ) )
from junctions_file import write_junctions
from genomic_intervals import GenomicInterval

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

def get_junctions( reads_fn, fasta ):
    """Get all of the junctions represented in a reads object
    """
    # build reads object
    reads = Samfile( reads_fn, "rb" )
    
    # store the number of reads across the junction for each relative 
    # read position
    junctions = defaultdict(int)
    
    for read in reads.fetch():
        # This could be changed to incorporate more complex reads 
        # indicating an intron such as including reads with clipped, 
        # inserted, or deleted regions or spanning multiple introns
        if not read_spans_single_intron( read ):
            continue
        
        # obtain chrm from read pysam object
        chrm = clean_chr_name( reads.getrname( read.tid ) )
        
        # add one to left_intron since bam files are 0-based
        upstrm_intron_pos = read.pos + read.cigar[0][1] + 1
        dnstrm_intron_pos = upstrm_intron_pos + read.cigar[1][1] - 1
        # get strand from stranded read information
        if (read.is_read1 and not read.is_reverse) \
                or (read.is_read2 and read.is_reverse):
            strand = '+'
        else:
            strand = '-'
        
        # increment count of junction reads at this read position
        # for this intron or initialize it to 1
        jn_type = ( chrm, strand, upstrm_intron_pos, \
                        dnstrm_intron_pos, read.cigar[0][1] )
        junctions[ jn_type ] += 1

    jns = []
    grps = []
    cnts = []
    for (chrm, strand, start, stop, grp ), cnt \
            in sorted(junctions.iteritems()):
        jns.append( GenomicInterval( chrm, strand, start, stop ) )
        grps.append( grp )
        cnts.append( cnt )
    
    return jns, grps, cnts


def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(description=\
       'Parse junctions from provided bam file and write to std out.')
    parser.add_argument( 'bam_fn', help="BAM file to parse jns from." )
    parser.add_argument( '--fasta', '-s', default=None,\
       help='Fasta file used to determine the consensus intron type' \
           + '(fasta file expects to be indexed: use samtools faidx).')    
    parser.add_argument( '--verbose', '-v', default=False, \
                             action='store_true', \
       help='Whether or not to print status information.' )
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose

    return args.bam_fn, args.fasta

def main():
    reads_fn, fasta_fn = parse_arguments()
    
    # get junctions
    jns, grps, cnts = get_junctions( reads_fn, fasta_fn )
   # write the junctions out to file
    write_junctions( jns, sys.stdout, cnts, grps, fasta_fn )

if __name__ == "__main__":
    main()
