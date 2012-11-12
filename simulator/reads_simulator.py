import sys
import os
import numpy
import pickle
import pysam
import math
from random import random

DEFAULT_QUALITY_SCORE = 'r'
DEFAULT_BASE = 'A'
DEFAULT_FRAG_LENGTH = 150
DEFAULT_READ_LENGTH = 100
DEFAULT_NUM_FRAGS = 100
NUM_NORM_SDS = 4
FREQ_GTF_STRINGS = [ 'freq', 'frac' ]

# add slide dir to sys.path and import frag_len mod
sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), "../sparsify/" ))
import frag_len

def calc_genomic_offset( transcript, offset ):
    """offset is a number of base pairs into the transcript
    
    loop through ordered exons until offset within transcript seq is reached
    """
    genomic_offset = 0

    prev_pos = 0
    for pos, exon_size in transcript.exons:
        # if offset does not extend through exon
        if offset < pos:
            genomic_offset += offset - prev_pos
            break
        
        # add exon and intron sizes to genomic_offset
        genomic_offset += exon_size
        # if pos is not in intron then this is the end of the transcript
        try:
            genomic_offset += transcript.introns[ pos ]
        except KeyError:
            break
        prev_pos = pos

    return genomic_offset

def get_cigar( transcript, start, stop ):
    """loop through introns within the read and add #N to the cigar for each intron
    add #M for portions of read which map to exons
    """
    cigar = ''
    prev_pos = start
    # loop though intron positions contined within this read 
    # adding skipped genomic position cigars
    for pos in sorted( set( range(start + 1, stop) ) & \
                           set( zip( *transcript.exons)[0] ) ):
        cigar += str( pos - prev_pos ) + "M" + \
            str( transcript.introns[ pos ] ) + "N"

        prev_pos = pos
    cigar += str( stop - prev_pos ) + "M"

    return cigar

def build_sam_line( transcript, read_len, offset, read_identifier, quality_string ):
    """build a single ended SAM formatted line with given inforamtion
    
    """
    # set flag to indcate strandedness of read matching that of the transcript
    flag = 0
    if transcript.strand == '-':
        flag += 16

    # adjust start position to correct genomic position
    start = transcript.start + calc_genomic_offset( transcript, offset )

    # set cigar string corresponding to transcript and read offset
    cigar = get_cigar( transcript, offset, (offset + read_len) )

    # calculate insert size by difference of genomic offset and genomic offset+read_len
    insert_size = calc_genomic_offset( transcript, (offset + read_len) ) - \
        calc_genomic_offset( transcript, offset )

    # get slice of seq from transcript
    seq = transcript.seq[ offset : (offset + read_len) ]

    # initialize sam lines with read identifiers and then add appropriate fields
    sam_line = '\t'.join( ( \
            read_identifier, str( flag ), transcript.chrm, str( start ), '255', \
                cigar, "*", '0', str( insert_size ), seq, quality_string, "NM:i:0", \
                "NH:i:1" )  ) + "\n"

    return sam_line

def build_sam_lines( transcript, read_len, frag_len, offset, \
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

    # adjust five and three prime read start positions to correct genomic positions
    start = [ transcript.start, transcript.start ]
    start[ up_strm_read ] += calc_genomic_offset( transcript, offset )
    start[ dn_strm_read ] += calc_genomic_offset( transcript, \
                                                 (offset + frag_len - read_len) )

    # set cigar string for five and three prime reads
    cigar = [ None, None ]
    cigar[ up_strm_read ] = get_cigar( transcript, offset, (offset + read_len) )
    cigar[ dn_strm_read ] = get_cigar( transcript, (offset + frag_len - read_len), \
                                      (offset + frag_len) )

    # calculate insert size by difference of genomic offset and genomic offset+frag_len
    insert_size = calc_genomic_offset( transcript, (offset + frag_len) ) - \
        calc_genomic_offset( transcript, offset )
    insert_size = [ insert_size, insert_size ]
    insert_size[ dn_strm_read ] *= -1

    # get slice of seq from transcript
    seq = [ None, None ]
    seq[ up_strm_read ] = transcript.seq[ offset : (offset + read_len) ]
    seq[ dn_strm_read ] = transcript.seq[ (offset + frag_len - read_len) : \
                                         (offset + frag_len) ]

    # initialize sam lines with read identifiers and then add appropriate fields
    sam_lines = [ read_identifier + '\t', read_identifier + '\t' ]    
    for i in (0,1):
        other_i = 0 if i else 1
        sam_lines[i] += '\t'.join( ( \
                str( flag[i] ), transcript.chrm, str( start[i] ), "255", cigar[i], \
                    "=", str( start[other_i] ), str( insert_size[i] ), seq[i], \
                    ordered_quals[i], "NM:i:0", "NH:i:1" ) ) + "\n"

    return sam_lines

def write_fastq_lines( fp1, fp2, transcript, read_len, frag_len, offset, \
                           read_identifier ):
    """STUB for writing fastq lines to running through alignment pipeline

    """
    pass

def simulate_reads( genes, fl_dist, fasta, quals, num_frags, single_end, \
                        full_fragment, read_len, out_prefix):
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
            assert fl_dist <= transcript.length, 'Transcript which ' + \
                'cannot contain a valid fragment was included in transcripts.'
            return fl_dist
        
        # Choose a valid fragment length from the distribution
        while True:
            fl_index = fl_dist.fl_density_cumsum.searchsorted( random() ) - 1
            fl = fl_index + fl_dist.fl_min
            
            # if fragment_length is valid return it
            if fl <= transcript.length:
                return fl
        assert False
            
    def sample_read_offset( transcript, fl ):
        # calculate maximum offset
        max_offset = transcript.length - fl + 1
        return int( random() * max_offset )
    
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
        read_identifier = 'SIM:%015d' % curr_read_index
        curr_read_index += 1
        
        return fl, offset, read_identifier

    def build_random_sam_line( transcript, read_len ):
        """build a random single ended sam line

        """
        fl, offset, read_identifier = get_random_read_pos( transcript )

        if full_fragment:
            read_len = fl
            
        # get a random quality scores
        read_qual = get_random_qual_score( read_len )

        # build the sam lines
        return build_sam_line( transcript, read_len, offset, read_identifier, read_qual )

    def build_random_sam_lines( transcript, read_len ):
        """build random paired end sam lines

        """
        fl, offset, read_identifier = get_random_read_pos( transcript )

        # adjust read length so that paired end read covers the entire fragment
        if full_fragment:
            read_len = int( math.ceil( fl / float(2) ) )

        # get two random quality scores
        read_quals = [ get_random_qual_score( read_len ), \
                           get_random_qual_score( read_len ) ]

        # build the sam lines
        return build_sam_lines( transcript, read_len, fl, \
                                    offset, read_identifier, read_quals )
    
    def get_fl_min():
        if isinstance( fl_dist, int ):
            return fl_dist
        else:
            return fl_dist.fl_min
    
    
    # initialize the transcript objects, and calculate their relative weights
    transcript_weights = []
    transcripts = []
    # count to assert there are basepairs across which to simulate
    total_bps = 0
    longest_trans = 0
    # sample a single transcript to get initial guess at shortest transcript
    shortest_trans = genes.values()[0].values()[0].get_length()

    min_transcript_length = get_fl_min()
    for gene in genes.values():
        for transcript in gene.values():
            # add the sequence information and intron/exon
            # structure to the transcript object
            transcript.process_transcript( fasta )

            # add the transcript and associated weight if there are valid fragments
            if transcript.length >= min_transcript_length:
                transcripts.append( transcript )
                # TODO:::for multiple genes this should be weighted across genes as well
                transcript_weights.append( \
                    transcript.freq * transcript.count_num_fragments( fl_dist ) )

            # update summary statistics
            if transcript.length > longest_trans:
                longest_trans = transcript.length
            if transcript.length < shortest_trans:
                shortest_trans = transcript.length
            total_bps += transcript.length
            
    assert total_bps > 0, "There are no bases across which to simulate reads!!!"
    assert len( transcripts ) > 0, "No valid trancripts."

    # normalize the transcript weights to be on 0,1
    tot_num_fragments = sum( transcript_weights )
    transcript_weights = [ num_fragments/tot_num_fragments \
                           for num_fragments in transcript_weights ]
    transcript_weights_cumsum = numpy.array( transcript_weights ).cumsum()

    if VERBOSE:
        avg_fl = fl_dist if isinstance( fl_dist, int ) else \
            (fl_dist.fl_min + fl_dist.fl_density_cumsum.searchsorted( 0.5 ) )
        fl_string = 'constant fragment length ' + str(avg_fl) + ' bps' if \
            isinstance( fl_dist, int ) else 'variable fragment length with median ' + \
            str(avg_fl)

        # determine approriate values for strings to be printed in summary
        ended_string = ' single-end ' if single_end else ' paired-end '

        read_string = 'full fragment ' if full_fragment else str(read_len) + ' bp '
        coverage = avg_fl * num_frags / float(total_bps) if full_fragment else \
            (2 * read_len * num_frags) / float(total_bps)

        # print summary inforamtion
        print '-' * 90
        print 'Simulating ' + str(num_frags) + ended_string + read_string + \
            'reads of ' + fl_string + ', producing approximately ' + \
            '{0:.3f}'.format( coverage ) + 'x coverage across all transcripts.'
        print '\tShortest transcript is ' + str(shortest_trans) + ' bps long.\n\t' + \
            'Longest transcript is ' + str(longest_trans) + ' bps long.'
        print '-' * 90

    sam_basename = out_prefix
    with open( sam_basename + ".sam", 'w' ) as sam_fp:
        while curr_read_index <= num_frags:
            # pick a transcript to randomly take a read from. Note that they should
            # be chosen in proportion to the *expected number of reads*, not their
            # relative frequencies.
            transcript_index = \
                transcript_weights_cumsum.searchsorted( random(), side='left' )
            transcript = transcripts[ transcript_index ]

            if single_end:
                sam_line_s = build_random_sam_line( transcript, read_len )
            else:
                sam_line_s = build_random_sam_lines( transcript, read_len )
            sam_fp.writelines( sam_line_s )
            
    # free the sequence information from the transcript objects
    # does not need to be done since all transcripts are processed together
    for transcript in transcripts:
        transcript.release()

    # create sorted bam file and index it
    sam_header_path = os.path.join( os.path.dirname(sys.argv[0]), 'sam_header.txt' )
    call = 'cat {0} {1}.sam | samtools view -bS - | samtools sort - {1}.sorted'
    os.system( call.format( sam_header_path, sam_basename ) )
    os.system( 'samtools index {0}.sorted.bam'.format( sam_basename ) )

    # make output directory and move sam and bam files into it
    os.mkdir( out_prefix )
    os.system( 'mv {0}.sam {0}.sorted.bam {0}.sorted.bam.bai {0}'.format( sam_basename ) )
    
    return

def call_slide( out_prefix, fl_dist_norm, fl_dist_const, gtf_fn ):
    """ Run slide on the simulated data and move the output to the output dir
    
    """
    # create paricular files and directory paths for options to slide
    slide_path = os.path.join( os.path.dirname(sys.argv[0]), '..', 'slide.py' )
    bam_path = os.path.join( out_prefix, out_prefix + '.sorted.bam' )
    fl_dist_norm_str = str(fl_dist_norm[0]) + ':' + str(fl_dist_norm[1]) \
        if fl_dist_norm else str(fl_dist_const) + ':' + '1'
    slide_out_dir = 'slide_' + out_prefix
    
    # call slide
    call = 'python {0} {1} {2} --fl_dist_norm {3} --out_prefix {4} --plot ' + \
        '--verbose --threads 1'
    os.system( call.format( slide_path, gtf_fn, bam_path, fl_dist_norm_str, \
                                slide_out_dir ) )
    os.system( 'mv {0} {1}'.format( slide_out_dir, out_prefix ) )

class Transcript( object ):
    def __init__( self, trans_name, data ):
        """Data contains chrm, strand, start, stop and freq (frac of major isoform)
        
        After process_transcript is run exons and introns are dicts 
        with key pos relative to the transcript seq (not genomic coordinates) and 
        value size of intron or exon
        """
        self.chrm = data[0]
        self.strand = data[1]
        self.freq = data[4]
        # if freq is 0 set it to 1/1000 to avoid division by 0
        if self.freq == 0:
            self.freq = 0.000001
        
        # a list of (start, stop) tuples
        self.exons = []
        # a list of intron positions relative to the transcript seq (key)
        # and corresonding length of intron (value)
        self.introns = {}
        # the sequence of the transcript from the reference
        self.seq = ''
    
    def add_exon( self, data ):
        assert data[0] == self.chrm
        assert data[1] == self.strand
        
        self.exons.append( ( data[2], data[3] ) )

    def get_length( self ):
        """Get transcript length before process_transcript has been run
        
        """
        length = 0
        for start, stop in self.exons:
            length += stop - start + 1
        return length

    def count_num_fragments( self, fl_dist ):
        # if a constant fragment length return the number of possible fragments
        # from this transcript
        if isinstance( fl_dist, int ):
            return self.get_length() - fl_dist + 1

        # else calculate the relative number of transcripts accounting for fl_density
        trans_len = self.get_length()
        num_fragments = 0.0
        for fl in xrange( fl_dist.fl_min, fl_dist.fl_max+1 ):
            num_fragments +=  (trans_len - fl + 1) * \
                fl_dist.fl_density[ fl - fl_dist.fl_min ]
        
        return num_fragments
    
    def process_transcript( self, fasta ):
        """produces transcript information from input exon values

        fasta is a pysam fasta file handle object accessed by fetch
        """
        self.start = min( zip( *self.exons )[0] )
        
        # used to assert no overlapping exons and calculate intron length
        prev_stop = -1
        # stores an exon end position relative to the seq and its length
        # key values to exons and introns are the same except for last exon key
        new_exons = []
        for start, stop in sorted( self.exons, key=lambda x:x[0] ):
            if start < prev_stop:
                raise ValueError, "Transcript " + self.name + \
                    " contains overlapping exons!!"

            # input intron position with respect to transcript seq and 
            # intron length into  intron_positions
            pos = len( self.seq )
            if pos:
                # if not first exon calculate intron size
                self.introns[ pos ] = start - prev_stop  - 1
            prev_stop = stop
            
            if fasta:
                # pysam uses 0-based indexing while gtf uses 1-based
                self.seq += fasta.fetch( self.chrm, start - 1, stop )
            else:
                # add appropriate number of DEFAULT_BASEs if fasta is not provided
                self.seq += DEFAULT_BASE * (stop - start + 1)
            new_exons.append( ( len( self.seq ), (stop - start + 1) ) )

        self.exons = new_exons
        # set exon length to the position of the last exon
        # corresponding to the lenght of the transcript
        self.length = self.exons[-1][0]
        
    def release( self ):
        self.seq = None
        self.exons = None
        self.introns = None
        
def parse_gtf_line( line  ):
    """Parse gtf line for either standard or custom gene object creation
    
    """
    data = line.split()
    # the type of element - ie exon, transcript, etc.
    gtf_type = data[2]
    
    # parse the meta data, and grab the gene name
    meta_data = dict( zip( data[8::2], data[9::2] ) )
    # remove ';' from value fields if ends in a ';'
    for key, val in meta_data.iteritems():
        if val.endswith( ';' ):
            meta_data[ key ] = val[:-1]
    gene_name = meta_data[ 'gene_id' ]
    trans_name = meta_data[ 'transcript_id' ]

    # set freq value from meta_data if available or default of 1/10,000
    freq=None
    for freq_string in FREQ_GTF_STRINGS:
        if freq_string in meta_data:
            freq = meta_data[ freq_string ]
            break
    if not freq:
        freq = 0.0001
        
    # fix meta_data to remove quotations
    if gene_name.startswith('"'):
        gene_name = gene_name[1:-1]
    if trans_name.startswith('"'):
        trans_name = trans_name[1:-1]
    if (not isinstance( freq, float )) and freq.startswith( '"' ):
        freq = float( freq[1:-1] )
    else:
        freq = float( freq )
    
    # fix the chromosome if necessary
    chrm = data[0]
    if chrm.startswith( 'chr' ):
        chrm = chrm[3:]

    # return names, type and data( chrm, strand, start, stop, freq )
    return gene_name, trans_name, gtf_type, \
        ( chrm, data[6], int(data[3]), int(data[4]), freq )

class Genes( dict ):
    """

    """
    def __init__( self, gtf_fp ):
        """Parse the gtf file at gtf_fp and return the genes and their unique exons.
        
        """
        self.filename = gtf_fp.name
        # TODO:::Add gene FPKM like value to determine 
        #        proportionof reads expected per gene
        
        # iterate through each line in the gtf, and organize exons
        # based upon the genes that they came from
        for line in gtf_fp:
            gene_name, trans_name, gtf_type, data = parse_gtf_line( line )

            if gtf_type == 'exon':
                try:
                    self[ gene_name ][ trans_name ].add_exon( data )
                except KeyError:
                    try:
                        self[ gene_name ][ trans_name ] = Transcript( trans_name, data )
                    except KeyError:
                        self[ gene_name ] = {}
                        self[ gene_name ][ trans_name ] = Transcript( trans_name, data )
                    self[ gene_name ][ trans_name ].add_exon( data )
                  
        gtf_fp.close()

def parse_custom_trans_string( custom_trans_string ):
    """Parse user input custom transcripts string
    
    user input must be of format:
        "exon_index,exon_index,...,exon_index:frequency|freq;exon_index..."
    user input is validated
    """
    transcripts = []
    for trans in custom_trans_string.split(";"):
        exons, freq = trans.split(":")

        try:
            freq = float( freq )
            exons = [ int(e) for e in exons.split(",") ]
            exons.sort()
        except ValueError:
            raise ValueError, "Custom gene transcript string is invalid!!!"
        transcripts.append( ( exons, freq ) )
    
    return transcripts

class CustomGene( dict ):
    """
    
    """
    def __init__( self, gtf_fp, custom_transcripts ):
        self.filename = gtf_fp.name
        exons = []

        for line in gtf_fp:
            gene_name, trans_name, gtf_type, data = parse_gtf_line( line )
            if gtf_type == 'exon':
                # set gene values from first exon found
                if not len( exons ):
                    self.chrm = data[0]
                    self.strand = data[1]
                    self.gene_name = gene_name
                # custom gene accepts a single gene so 
                # only get exons from the first gene in the gtf file
                if data[0] != self.chrm or \
                        data[1] != self.strand or \
                        gene_name != self.gene_name:
                    break
                # add exon start and stop to exons list
                exons.append( ( data[2], data[3] ) )
        
        gtf_fp.close()
        exons.sort()
        assert len( exons ), "No exons in first gene in " + self.filename

        # create data structure same as in Genes constructor
        self[ self.gene_name ] = {}
        for trans_num, (exon_indices, freq) in enumerate( custom_transcripts ):
            trans_name = "custom_trans_" + str( trans_num )
            self[ self.gene_name ][ trans_name ] = \
                Transcript( trans_name, (self.chrm, self.strand, 0, 0, freq) )
            
            # add exons to the transcript
            prev_stop = -1
            for start, stop in [ exons[i] for i in exon_indices ]:
                if prev_stop >= start:
                    raise ValueError, "A custom transcript contains overlapping exons!!!"
                prev_stop = stop
                
                data = (self.chrm, self.strand, start, stop, freq )
                self[ self.gene_name ][ trans_name ].add_exon( data )

def build_objs( gtf_fp, fl_dist_const, fl_dist_obj_fn, fl_dist_norm, full_fragment, \
                    read_len, fasta_fn, qual_fn, custom_trans ):
    if custom_trans:
        custom_trans = parse_custom_trans_string( custom_trans_string )
        genes = CustomGene( gtf_fp, custom_trans )
    else:
        genes = Genes( gtf_fp )
    gtf_fp.close()

    def build_normal_fl_dist( fl_mean, fl_sd ):
        fl_min = max( 0, fl_mean - (fl_sd * NUM_NORM_SDS) )
        fl_max = fl_mean + (fl_sd * NUM_NORM_SDS)
        fl_dist = frag_len.build_normal_density( fl_min, fl_max, fl_mean, fl_sd )
        return fl_dist

    if fl_dist_obj_fn:
        # if provided load fl_dist obj from file
        with open( fl_dist_fn ) as fl_dist_fp:
            fl_dist = pickle.load( fl_dist_fp )
        assert isinstance( fl_dist, frag_len.FlDist ), 'fl_dist_obj option must ' + \
            'be a pickeled FlDist object!!!'
        assert fl_dist.fl_max > read_len or full_fragment, \
            'Invalid fragment length distribution and read length!!!'
    elif fl_dist_norm:
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

    # fragment length distribution options
    parser.add_argument( '--fl_dist_const', type=int, default=DEFAULT_FRAG_LENGTH, \
                             help='Constant length fragments. (default: ' + \
                             '%(default)s)' )
    parser.add_argument( '--fl_dist_obj', \
                             help='Cached fragment length distribution object ' + \
                             'filename, created with get_frag_dist.py.' )
    parser.add_argument( '--fl_dist_norm', \
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
    parser.add_argument( '--num_frags', '-n', type=int, default=DEFAULT_NUM_FRAGS, \
                             help='Total number of fragments to create across all ' + \
                             'transcripts (default: %(default)s)' )
    parser.add_argument( '--single_end', action='store_true', default=False, \
                             help='Produce single-end reads.' )    
    parser.add_argument( '--paired_end', dest='single_end', action='store_false', \
                             help='Produce paired-end reads. (default)' )    
    parser.add_argument( '--full_fragment', action='store_true', default=False, \
                             help='Produce reads spanning the entire fragment. (If ' + \
                             'used in conjunction with paired_end option reads will ' + \
                             'cover approx. half the fragment each) Note: read_len ' + \
                             'option will be ignored if this option is set.' )
    parser.add_argument( '--read_len', '-r', type=int, default=DEFAULT_READ_LENGTH, \
                             help='Length of reads to produce in base pairs ' + \
                             '(default: %(default)s)' )

    # output options
    parser.add_argument( '--out_prefix', '-o', default='simulated_reads', \
                             help='Prefix for output FASTQ/BAM file ' + \
                             '(default: %(default)s)' )
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', \
                             help='Print status information.' )
    parser.add_argument( '--run_slide', action='store_true', default=False, \
                             help='Run slide on the simulated gtf file and the ' + \
                             'simulated bam file. (default: False) Note: Not ' + \
                             'compatible with fl_dist_obj option.' )

    # custom transcripts option
    parser.add_argument( '--custom_trans', \
                             help='User defined transcripts created from only the ' + \
                             'first gene in the gtf file in the format: ' + \
                             'exon_index,...,exon_index:frequency;exon_index,...')

    try:
        args = parser.parse_args()
    except IOError:
        raise IOError, 'Invalid GTF file provided!!!'

    global VERBOSE
    VERBOSE = args.verbose

    # parse normal distribution argument
    if args.fl_dist_norm:
        try:
            mean, sd = args.fl_dist_norm.split( ':' )
            args.fl_dist_norm = [ int( mean ), int( sd ) ]
        except ValueError:
            args.fl_dist_norm = None
            print "WARNING: User input mean and sd are not formatted correctly.\n" + \
                "\tUsing default values.\n"

    return args.gtf, args.fl_dist_const, args.fl_dist_obj, args.fl_dist_norm, \
        args.fasta, args.quality, args.num_frags, args.single_end, args.full_fragment, \
        args.read_len, args.out_prefix, args.custom_trans, args.run_slide

if __name__ == "__main__":
    gtf_fp, fl_dist_const, fl_dist_obj_fn, fl_dist_norm, fasta_fn, qual_fn, \
        num_frags, single_end, full_fragment, read_len, out_prefix, custom_trans, \
        run_slide = parse_arguments()
    
    genes, fl_dist, fasta, quals = build_objs( \
        gtf_fp, fl_dist_const, fl_dist_obj_fn, fl_dist_norm, full_fragment, read_len, \
            fasta_fn, qual_fn, custom_trans ) 

    """Test a single read for debugging
    transcript = genes.values()[0].values()[0]
    transcript.process_transcript( None )
    lines = build_sam_lines( transcript, 100, 700, 31, 'read1', ['r','r'] )
    for line in lines:
        print line,
    sys.exit()
    """

    simulate_reads( genes, fl_dist, fasta, quals, num_frags, single_end, full_fragment, 
                    read_len, out_prefix )
    
    if run_slide:
        call_slide( out_prefix, fl_dist_norm, fl_dist_const, gtf_fp.name )
