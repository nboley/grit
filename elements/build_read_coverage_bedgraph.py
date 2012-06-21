# Copyright (c) 2011-2012 Nathan Boley

import sys, os
import pysam

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../file_types/" ) )
from wiggle import Wiggle

def clean_chr_name( chrm ):
    if chrm.startswith( "chr" ):
        chrm = chrm[3:]
    # convert the dmel chrm to M, to be consistent with ucsc
    if chrm == 'dmel_mitochondrion_genome':
        chrm = "M"
    return chrm

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
        
        # make sure that the read is on the correct strand
        if ( read.is_read1 and not read.is_reverse ) \
                or ( read.is_read2 and read.is_reverse ):
            strand = '+'
        else:
            strand = '-'
        
        if reverse_read_strand:
            if strand == '+':
                strand = '-'
            else:
                strand = '+'
            
        
        # get the chromosome, correcting for alternate chrm names
        chrm = reads.getrname( read.tid )
        chrm = clean_chr_name( chrm )
        
        # we loop through each contig in the cigar string to deal
        # with junctions reads.
        # add 1 to the start because bam files are 0 based
        start = read.pos + 1
        for contig_type, length in read.cigar:
            # if this is a match, add it 
            if contig_type == 0:
                wiggle.add_cvg( chrm, strand, start, start + length - 1, 1.0 )
                start += length
            # skip reference insertions
            elif contig_type == 1:
                pass
            # move past refernce deletions
            elif contig_type == 2:
                start += length
            # skip past skipped regions
            elif contig_type == 3:
                start += length
            # skip past soft clipped regions, because the
            # actual aligned sequence doesnt start until we've moved
            # past the clipped region
            elif contig_type == 4:
                start += length
            # hard clipped regions are not present int he aligned 
            # sequence, so do nothing
            elif contig_type == 5:
                pass
            else:
                print >> sys.stderr, "Unrecognized cigar format:", read.cigar

        
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

