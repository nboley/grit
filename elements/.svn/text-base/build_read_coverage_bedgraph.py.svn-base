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

def populate_wiggle( reads, wiggle, reverse_read_strand=False ):
    """Get all of the junctions represented in a reads object
    """
    for read in reads.fetch():
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
    parser.add_argument( '--out-fname-prefix', '-o', required=True,\
        help='Output files will be named outprefix.(plus,minus).bedGraph')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', \
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.mapped_reads_fname, args.chrm_sizes, args.out_fname_prefix, \
        args.reverse_read_strand

def main():
    reads_fname, chrm_sizes, op_prefix, reverse_read_strand = parse_arguments()
    
    # open the reads object
    reads = pysam.Samfile( reads_fname, "rb" )
    # build a new wiggle object to put the reads in
    wiggle = Wiggle( chrm_sizes )
    # populate the wiggle from the bam file
    populate_wiggle( reads, wiggle, reverse_read_strand )
    # write the wiggle to disk
    wiggle.write_wiggles( "{0}.{1}.bedGraph".format( op_prefix, "plus"), \
                          "{0}.{1}.bedGraph".format( op_prefix, "minus"),
                          ignore_zeros=True)
    return

if __name__ == "__main__":
    main()

