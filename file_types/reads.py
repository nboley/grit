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

def read_pairs_are_on_same_strand( bam_obj, num_reads_to_check=100 ):
    # keep track of which fractiona re on the sam strand
    paired_cnts = {'no_mate': 0, 'same_strand': 1e-4, 'diff_strand': 1e-4}
    
    num_good_reads = 0
    num_observed_reads = 0
    for read in bam_obj:
        num_observed_reads += 1
        if read.is_paired and read.mate_is_unmapped:
            continue

        if not read.is_paired:
            paired_cnts['no_mate'] += 1        
        elif read.is_reverse != read.mate_is_reverse:
            paired_cnts['diff_strand'] += 1
        else:
            paired_cnts['same_strand'] += 1
            
        # keep collecting reads until we observe enough
        num_good_reads += 1
        if num_good_reads > num_reads_to_check:
            break
    
    # if the reads are single ended, then return True ( 
    #    because it doesnt really matter )
    if paired_cnts['no_mate'] == num_reads_to_check:
        return True
    
    if float(paired_cnts['same_strand'])/paired_cnts['diff_strand'] > 10:
        return True
    elif float(paired_cnts['diff_strand'])/paired_cnts['same_strand'] > 10:
        return False
    else:
        print >> sys.stderr, paired_cnts, num_observed_reads
        raise ValueError, "Reads appear to be a mix of reads on the same and different strands."

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

