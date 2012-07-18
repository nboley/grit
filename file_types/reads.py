import sys, os

def clean_chr_name( chrm ):
    if chrm.startswith( "chr" ):
        chrm = chrm[3:]
    # convert the dmel chrm to M, to be consistent with ucsc
    if chrm == 'dmel_mitochondrion_genome':
        chrm = "M"
    return chrm

def get_chrm( read, bam_obj ):
    chrm = bam_obj.getrname( read.tid )
    chrm = clean_chr_name( chrm )
    return chrm

def get_strand( read, reverse_read_strand, pairs_are_opp_strand ):
    # make sure that the read is on the correct strand
    if pairs_are_opp_strand  \
            and ( ( read.is_read1 and not read.is_reverse ) \
                      or ( read.is_read2 and read.is_reverse ) ):
        strand = '+'
    elif not pairs_are_opp_strand and not read.is_reverse:
        strand = '+'
    else:
        strand = '-'

    if reverse_read_strand:
        if strand == '+':
            strand = '-'
        else:
            strand = '+'

    return strand

def iter_coverage_regions_for_read( 
    read, bam_obj, reverse_read_strand, pairs_are_opp_strand ):
    """Find the regions covered by this read

    """
    strand = get_strand( read, reverse_read_strand, pairs_are_opp_strand )

    # get the chromosome, correcting for alternate chrm names
    chrm = get_chrm( read, bam_obj )

    # we loop through each contig in the cigar string to deal
    # with junctions reads.
    # add 1 to the start because bam files are 0 based
    rv = []
    start = read.pos
    for contig_type, length in read.cigar:
        # if this is a match, add it 
        if contig_type == 0:
            yield ( chrm, strand, start, start + length - 1 )
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

    return

