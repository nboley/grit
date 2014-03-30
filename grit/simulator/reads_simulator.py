import sys
import os
import os.path
import numpy
import pickle
import pysam
import math
from random import random
from collections import defaultdict
import tempfile

DEFAULT_QUALITY_SCORE = 'r'
DEFAULT_BASE = 'A'
DEFAULT_FRAG_LENGTH = 150
DEFAULT_READ_LENGTH = 100
DEFAULT_NUM_FRAGS = 100
NUM_NORM_SDS = 4
FREQ_GTF_STRINGS = [ 'freq', 'frac' ]

# add slide dir to sys.path and import frag_len mod
#sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), ".." ))
import grit.frag_len as frag_len
from grit.files.gtf import load_gtf
from grit.files.reads import clean_chr_name

def fix_chr_name(x):
    return "chr" + clean_chr_name(x)

def get_transcript_sequence(transcript, fasta):
    """ get the mRNA sequence of the transcript from the gene seq
    """
    trans_seq = []
    for start, stop in transcript.exons:
        seq = fasta.fetch(fix_chr_name(transcript.chrm), start, stop+1)
        trans_seq.append( seq.upper() )
    
    trans_seq = "".join(trans_seq)
    return trans_seq


def get_cigar( transcript, start, stop ):
    """loop through introns within the read and add #N to the cigar for each 
    intron add #M for portions of read which map to exons
    """
    def calc_len(interval):
        return interval[1]-interval[0]+1
    cigar = []
    
    # find the exon index of the start
    genome_start = transcript.genome_pos(start)
    start_exon = next(i for i, (e_start, e_stop) in enumerate(transcript.exons)
                      if genome_start >= e_start and genome_start <= e_stop)
    genome_stop = transcript.genome_pos(stop-1)
    stop_exon = next(i for i, (e_start, e_stop) in enumerate(transcript.exons)
                     if genome_stop >= e_start and genome_stop <= e_stop)

    if start_exon == stop_exon:
        return "%iM" % (stop-start)

    tl = 0
    # add the first overlap match
    skipped_bases = sum(calc_len(e) for e in transcript.exons[:start_exon+1])
    cigar.append("%iM" % (skipped_bases-start))
    tl += skipped_bases-start
    
    # add the first overlap intron 
    cigar.append("%iN" % calc_len(transcript.introns[start_exon]))
    
    # add the internal exon and intron matches
    for i in xrange(start_exon+1, stop_exon):
        cigar.append("%iM" % calc_len(transcript.exons[i]))
        cigar.append("%iN" % calc_len(transcript.introns[i]))
        tl += calc_len(transcript.exons[i])

    # add the last overlap match
    skipped_bases = sum(e[1]-e[0]+1 for e in transcript.exons[:stop_exon])
    cigar.append("%iM" % (stop-skipped_bases))
    tl += stop - skipped_bases 
    
    assert tl == (stop-start)
    
    return "".join(cigar)

def build_sam_line( transcript, read_len, offset, read_identifier, quality_string ):
    """build a single ended SAM formatted line with given inforamtion
    
    """
    # set flag to indcate strandedness of read matching that of the transcript
    flag = 0
    if transcript.strand == '-': flag += 16
    
    # adjust start position to correct genomic position
    start = transcript.genome_pos(offset)
    # set cigar string corresponding to transcript and read offset
    cigar = get_cigar( transcript, offset, (offset + read_len) )
    # calculate insert size by difference of genomic offset and genomic offset+read_len
    insert_size = transcript.genome_pos(offset+read_len) - transcript.genome_pos(offset)
    # get slice of seq from transcript
    seq = ( transcript.seq[ offset : (offset + read_len) ] 
            if transcript.seq != None else '*' )
    # initialize sam lines with read identifiers and then add appropriate fields
    sam_line = '\t'.join( ( 
            read_identifier, str( flag ), fix_chr_name(transcript.chrm), 
            str(start+1),
            '255', cigar, "*", '0', str( insert_size ), seq, quality_string, 
            "NM:i:0", "NH:i:1" )  ) + "\n"

    return sam_line

def build_sam_lines( transcript, read_len, frag_len, offset, 
                     read_identifier, read_quals ):
    """build paired end SAM formatted lines with given information

    """
    # set ordered quals and reverse the qualities for the read on the negative strand
    ordered_quals = read_quals
    
    # determine whether read1 should be the 5' read or visa verses
    # and initialize attributes that are specific to a read number 
    # instead of 5' or 3' attribute
    if transcript.strand == '+':
        up_strm_read, dn_strm_read = (0, 1)
        flag = [ 99, 147 ]
        ordered_quals[1] = ordered_quals[1][::-1]
    else:
        up_strm_read, dn_strm_read = (1, 0)
        flag = [ 83, 163 ]
        ordered_quals[0] = ordered_quals[0][::-1]

    # get slice of seq from transcript
    seq = ['*', '*']
    if transcript.seq != None:
        seq[ up_strm_read ] = transcript.seq[offset:(offset + read_len)]
        seq[ dn_strm_read ] = transcript.seq[
            (offset + frag_len - read_len):(offset + frag_len)]
    
    # adjust five and three prime read start positions to correct genomic positions
    start = [ transcript.start, transcript.start ]
    start[ up_strm_read ] = transcript.genome_pos(offset)
    start[ dn_strm_read ] = transcript.genome_pos(offset + frag_len - read_len)
    
    # set cigar string for five and three prime reads
    cigar = [ None, None ]
    cigar[ up_strm_read ] = get_cigar( transcript, offset, (offset+read_len) )
    cigar[ dn_strm_read ] = get_cigar( 
        transcript, (offset+frag_len-read_len), (offset + frag_len))
    
    # calculate insert size by difference of the mapped start and end
    insert_size = (
        transcript.genome_pos(offset+read_len) - transcript.genome_pos(offset))
    insert_size = [ insert_size, insert_size ]
    insert_size[ dn_strm_read ] *= -1
    
    # initialize sam lines with read identifiers and then add appropriate fields
    sam_lines = [ read_identifier + '\t', read_identifier + '\t' ]    
    for i in (0,1):
        other_i = 0 if i else 1
        sam_lines[i] += '\t'.join( (
                str( flag[i] ), fix_chr_name(transcript.chrm), 
                str( start[i]+1 ),"255",
                cigar[i], "=", str( start[other_i] ), str( insert_size[i] ), 
                seq[i], ordered_quals[i], "NM:i:0", "NH:i:1" ) ) + "\n"

    return sam_lines

def write_fastq_lines( fp1, fp2, transcript, read_len, frag_len, offset, 
                       read_identifier ):
    """STUB for writing fastq lines to running through alignment pipeline

    """
    pass

def simulate_reads( genes, fl_dist, fasta, quals, num_frags, single_end, 
                    full_fragment, read_len, assay='RNAseq'):
    """write a SAM format file with the specified options

    """    
    # global variable that stores the current read number, we use this to 
    # generate a unique id for each read.
    global curr_read_index
    curr_read_index = 1
    
    def sample_fragment_length( fl_dist, transcript ):
        """Choose a random fragment length from fl_dist

        """
        # if the fl_dist is constant
        if isinstance( fl_dist, int ):
            assert fl_dist <= transcript.calc_length(), 'Transcript which ' + \
                'cannot contain a valid fragment was included in transcripts.'
            return fl_dist
        
        # Choose a valid fragment length from the distribution
        while True:
            fl_index = fl_dist.fl_density_cumsum.searchsorted( random() ) - 1
            fl = fl_index + fl_dist.fl_min
            
            # if fragment_length is valid return it
            if fl <= transcript.calc_length():
                return fl
        assert False
        
    def sample_read_offset( transcript, fl ):
        # calculate maximum offset
        max_offset = transcript.calc_length() - fl
        if assay == 'RNAseq':
            return int( random() * max_offset )
        elif assay == 'RAMPAGE':
            return 0
        elif assay == 'CAGE':
            return 0
        elif assay == 'PASseq':
            return max_offset
    
    def get_random_qual_score( read_len ):
        # if no quality score were provided
        if not quals:
            return DEFAULT_QUALITY_SCORE * read_len
        # else return quality string from input quality file 
        # scores are concatenated to match read_len if necessary
        else:
            qual_string = ''
            while len( qual_string ) < read_len:
                qual_string += str( quals[ int(random() * len(quals) ) ] )
            return qual_string[0:read_len]
    
    def get_random_read_pos( transcript ):
        while True:
            # find a valid fragment length
            fl = sample_fragment_length( fl_dist, transcript )
            if (fl >= read_len) or full_fragment: break

        # find a valid random read start position
        offset = sample_read_offset( transcript, fl )
        
        # get a unique string for this fragment
        global curr_read_index
        read_identifier = 'SIM:%015d:%s' % (curr_read_index, transcript.id)
        curr_read_index += 1
        
        return fl, offset, read_identifier

    def build_random_sam_line( transcript, read_len ):
        """build a random single ended sam line

        """
        fl, offset, read_identifier = get_random_read_pos( transcript )

        if full_fragment:
            read_len = fl
            
        # get a random quality scores
        if transcript.seq == None:
            read_qual = '*'
        else:
            read_qual = get_random_qual_score( read_len )

        # build the sam lines
        return build_sam_line( 
            transcript, read_len, offset, read_identifier, read_qual )

    def build_random_sam_lines( transcript, read_len ):
        """build random paired end sam lines

        """
        fl, offset, read_identifier = get_random_read_pos( transcript )

        # adjust read length so that paired end read covers the entire fragment
        if full_fragment:
            read_len = int( math.ceil( fl / float(2) ) )

        # get two random quality scores
        if transcript.seq == None:
            read_quals = ['*', '*']
        else:
            read_quals = [ get_random_qual_score( read_len ), 
                           get_random_qual_score( read_len ) ]

        sam_lines = build_sam_lines( 
            transcript, read_len, fl, offset, read_identifier, read_quals )
        
        return sam_lines
    
    def get_fl_min():
        if isinstance( fl_dist, int ):
            return fl_dist
        else:
            return fl_dist.fl_min
    
    def calc_scale_factor(t):
        if assay in ('RNAseq',):
            length = t.calc_length()
            if length < fl_dist.fl_min: return 0
            fl_min, fl_max = fl_dist.fl_min, min(length, fl_dist.fl_max)
            allowed_fl_lens = numpy.arange(fl_min, fl_max+1)
            weights = fl_dist.fl_density[
                fl_min-fl_dist.fl_min:fl_max-fl_dist.fl_min+1]
            mean_fl_len = float((allowed_fl_lens*weights).sum())
            return length - mean_fl_len
        elif assay in ('CAGE', 'RAMPAGE', 'PASseq'):
            return 1.0
    
    # initialize the transcript objects, and calculate their relative weights
    transcript_weights = []
    transcripts = []
    
    contig_lens = defaultdict(int)
    
    min_transcript_length = get_fl_min()
    for gene in genes:
        contig_lens[fix_chr_name(gene.chrm)] = max(
            gene.stop, contig_lens[fix_chr_name(gene.chrm)])
        for transcript in gene.transcripts:
            if fasta != None:
                transcript.seq = get_transcript_sequence(transcript, fasta)
            else:
                transcript.seq = None
            if transcript.fpkm != None:
                weight = transcript.fpkm*calc_scale_factor(transcript)
            elif transcript.frac != None:
                assert len(genes) == 1
                weight = transcript.frac
            else: assert False, "Transcript has neither an FPKM nor a frac"
            transcripts.append( transcript )
            transcript_weights.append( weight )
    
    #assert False
    assert len( transcripts ) > 0, "No valid trancripts."

    # normalize the transcript weights to be on 0,1
    transcript_weights = numpy.array(transcript_weights, dtype=float)
    transcript_weights = transcript_weights/transcript_weights.sum()
    transcript_weights_cumsum = transcript_weights.cumsum()

    # update the contig lens from the fasta file, if available 
    if fasta != None:
        for name, length in zip(fasta.references, fasta.lengths):
            contig_lens[fix_chr_name(name)] = max(
                length, contig_lens[name])

    # create the output directory
    bam_prefix = assay + ".sorted"
    
    with tempfile.NamedTemporaryFile( mode='w+' ) as sam_fp:
        # write out the header
        for contig, contig_len in contig_lens.iteritems():
            contig_len = 22422827
            data = ["@SQ", "SN:%s" % contig, "LN%i" % contig_len]
            sam_fp.write("\t".join(data) + "\n")
        
        while curr_read_index <= num_frags:
            # pick a transcript to randomly take a read from. Note that they 
            # should be chosen in proportion to the *expected number of reads*,
            # not their relative frequencies.
            transcript_index = \
                transcript_weights_cumsum.searchsorted( random(), side='left' )
            transcript = transcripts[ transcript_index ]

            if single_end:
                sam_line_s = build_random_sam_line( transcript, read_len )
            else:
                sam_line_s = build_random_sam_lines( transcript, read_len )
            sam_fp.writelines( sam_line_s )
    
        # create sorted bam file and index it
        sam_fp.flush()
        #sam_fp.seek(0)
        #print sam_fp.read()
        
        call = 'samtools view -bS {} | samtools sort - {}'
        os.system( call.format( sam_fp.name, bam_prefix ) )
        os.system( 'samtools index {}.bam'.format( bam_prefix ) )
        
    return
        
def build_objs( gtf_fp, fl_dist_const, 
                fl_dist_norm, full_fragment,
                read_len, fasta_fn, qual_fn ):

    genes = load_gtf( gtf_fp )
    gtf_fp.close()

    def build_normal_fl_dist( fl_mean, fl_sd ):
        fl_min = max( 0, fl_mean - (fl_sd * NUM_NORM_SDS) )
        fl_max = fl_mean + (fl_sd * NUM_NORM_SDS)
        fl_dist = frag_len.build_normal_density( fl_min, fl_max, fl_mean, fl_sd )
        return fl_dist

    if fl_dist_norm:
        fl_dist = build_normal_fl_dist( fl_dist_norm[0], fl_dist_norm[1] )
        assert fl_dist.fl_max > read_len or full_fragment, \
            'Invalid fragment length distribution and read length!!!'
    else:
        assert read_len < fl_dist_const or full_fragment, \
            'Invalid read length and constant fragment length!!!'
        fl_dist = fl_dist_const

    if fasta_fn:
        # create indexed fasta file handle object with pysam
        fasta = pysam.Fastafile( fasta_fn )
    else:
        fasta = None

    # if quals_fn is None, quals remains empty and reads will default to 
    # all base qualities of DEFAULT_BASE_QUALITY_SCORE
    quals = []
    if qual_fn:
        with open( quals_fn ) as quals_fp:
            for line in quals_fp:
                quals.append( line.strip() )
            quals = numpy.array( quals )

    return genes, fl_dist, fasta, quals

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(\
        description='Produce simulated reads in a perfecty aligned BAM file.' )
    
    # gtf is the only required argument
    parser.add_argument( 'gtf', type=file, \
                             help='GTF file from which to produce simulated reads ' + \
                             '(Note: Only the first trascript from this file will ' + \
                             'be simulated)' )

    parser.add_argument( 
        '--assay', choices=['RNAseq', 'RAMPAGE', 'CAGE', 'PASseq'],
        default='RNAseq', help='Which assay type to simulate from' )
    
    # fragment length distribution options
    parser.add_argument( '--fl-dist-const', type=int, default=DEFAULT_FRAG_LENGTH, \
                             help='Constant length fragments. (default: ' + \
                             '%(default)s)' )
    parser.add_argument( '--fl-dist-norm', \
                             help='Mean and standard deviation (format "mn:sd") ' + \
                             'used to create normally distributed fragment lengths.' )

    # files providing quality and sequnce information
    parser.add_argument( '--fasta', '-f', \
                             help='Fasta file from which to create reads ' + \
                             '(default: all sequences are "' + DEFAULT_BASE + \
                             '" * length of sequence)' )
    parser.add_argument( '--quality', '-q', \
                             help='Flat file containing one FASTQ quality score ' + \
                             'per line, created with get_quals.sh. (default: ' + \
                             'quality strings are  "' + str(DEFAULT_QUALITY_SCORE) + \
                             '" * length of sequence.)' )

    # type and number of fragments requested
    parser.add_argument( 
        '--num-frags', '-n', type=int, default=1000,
        help='Total number of fragments to create across all trascripts')
    parser.add_argument('--single-end', action='store_true', default=False, 
                        help='Produce single-end reads.' )    
    parser.add_argument('--paired-end', dest='single_end', action='store_false',
                        help='Produce paired-end reads. (default)' )    
    # XXX not sure if this works
    #parser.add_argument( 
    #    '--full-fragment', action='store_true', default=False, 
    #    help='Produce reads spanning the entire fragment.')
    
    parser.add_argument( '--read-len', '-r', type=int, default=DEFAULT_READ_LENGTH, \
                             help='Length of reads to produce in base pairs ' + \
                             '(default: %(default)s)' )

    # output options
    parser.add_argument( '--out_prefix', '-o', default='simulated_reads', \
                             help='Prefix for output FASTQ/BAM file ' + \
                             '(default: %(default)s)' )
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', \
                             help='Print status information.' )
    
    args = parser.parse_args()
    # set to false, but we may want to bring this option back
    args.full_fragment = False
    
    global VERBOSE
    VERBOSE = args.verbose

    if args.assay == 'CAGE':
        args.read_len = 28
        args.single_end = True
        
    # parse normal distribution argument
    if args.fl_dist_norm:
        try:
            mean, sd = args.fl_dist_norm.split( ':' )
            args.fl_dist_norm = [ int( mean ), int( sd ) ]
        except ValueError:
            args.fl_dist_norm = None
            print >> sys.stderr, \
              "WARNING: User input mean and sd are not formatted correctly.\n"+\
              "\tUsing default values.\n"

    return ( args.gtf, args.fl_dist_const, args.fl_dist_norm, 
             args.fasta, args.quality, args.num_frags, 
             args.single_end, args.full_fragment, 
             args.read_len, args.out_prefix, args.assay )

def main():
    ( gtf_fp, fl_dist_const, fl_dist_norm, fasta_fn, qual_fn, 
      num_frags, single_end, full_fragment, read_len, out_prefix, assay )\
        = parse_arguments()
    
    try: os.mkdir(out_prefix)
    except OSError: 
        ofname = os.path.join(out_prefix, assay + '.sorted.bam')
        if os.path.isfile(ofname):
            raise OSError, "File '%s' already exists" % ofname
    os.chdir(out_prefix)
    
    genes, fl_dist, fasta, quals = build_objs( 
        gtf_fp, fl_dist_const, 
        fl_dist_norm, full_fragment, read_len, 
        fasta_fn, qual_fn ) 
    
    """
    for gene in genes:
        for t in gene.transcripts:
            t.chrm = "chr" + t.chrm
            print t.build_gtf_lines(gene.id, {})
    assert False
    """
    simulate_reads( genes, fl_dist, fasta, quals, num_frags, single_end, 
                    full_fragment, read_len, assay=assay )

if __name__ == "__main__":
    main()
