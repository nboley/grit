# Copyright (c) 2011-2012 Nathan Boley

import sys, os
import pysam

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../file_types/" ) )
from wiggle import Wiggle
from reads import iter_coverage_regions_for_read

def populate_wiggle( reads, rd1_wig, rd2_wig, reverse_read_strand=False ):
    """Get all of the junctions represented in a reads object
    """
    for read in reads.fetch():
        # skip reads that aren't mapped in pair
        if not read.is_proper_pair:
            continue

        # find which wiggle to write the read to
        if rd2_wig == None:
            wiggle = rd1_wig
        else:
            if read.is_read1 or ( read.is_read2 and reverse_read_strand ):
                wiggle = rd1_wig
            else:
                wiggle = rd2_wig
                
        for chrm, strand, start, stop in iter_coverage_regions_for_read( 
                read, reads, reverse_read_strand ):
            wiggle.add_cvg( chrm, strand, start, stop, 1.0 )
    
    reads.close()
    
    return

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(
        description='Get coverage bedgraphs from aligned reads.')
    parser.add_argument( '--mapped-reads-fname', required=True,
        help='BAM or SAM file containing the mapped reads.')
    parser.add_argument( '--chrm-sizes', '-c', required=True, type=file, \
        help='File with chromosome names and sizes.')
    
    parser.add_argument( '--reverse-read-strand', '-r', default=False, \
                             action='store_true', \
        help='Whether or not to reverse the strand of the read.')
    parser.add_argument( '--merge-read-ends', '-m', default=False, \
                             action='store_true', \
        help='Whether or not to merge pair1 and pair2 reads for paired reads.')
    
    parser.add_argument( '--out-fname-prefix', '-o', required=True,\
        help='Output files will be named outprefix.(plus,minus).bedGraph')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.mapped_reads_fname, args.chrm_sizes, args.out_fname_prefix, \
        args.reverse_read_strand, args.merge_read_ends

def main():
    reads_fname, chrm_sizes, op_prefix, reverse_read_strand, merge_read_ends \
        = parse_arguments()
    
    # open the reads object
    reads = pysam.Samfile( reads_fname, "rb" )
    # build a new wiggle object to put the reads in
    rd1_wig = Wiggle( chrm_sizes )
    rd2_wig = Wiggle( chrm_sizes ) if not merge_read_ends else None
    
    # populate the wiggle from the bam file
    populate_wiggle( reads, rd1_wig, rd2_wig, reverse_read_strand )
    
    # write the wiggle to disk
    if merge_read_ends:
        rd1_wig.write_wiggles( "{0}.{1}.bedGraph".format( op_prefix, "plus"),
                               "{0}.{1}.bedGraph".format( op_prefix, "minus"),
                              ignore_zeros=True)
    else:
        rd1_wig.write_wiggles( "{0}.{1}.{2}.bedGraph".format( op_prefix, "rd1", "plus"),
                               "{0}.{1}.{2}.bedGraph".format( op_prefix, "rd1", "minus"),
                              ignore_zeros=True)
        rd2_wig.write_wiggles( "{0}.{1}.{2}.bedGraph".format( op_prefix, "rd2", "plus"),
                               "{0}.{1}.{2}.bedGraph".format( op_prefix, "rd2", "minus"),
                              ignore_zeros=True)
    
    return

if __name__ == "__main__":
    main()

