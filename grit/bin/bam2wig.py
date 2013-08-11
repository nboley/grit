import sys, os
import pysam
import numpy
import shutil
import subprocess
import tempfile
import time
from itertools import izip

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), ".." ) )
from grit.files.reads import iter_coverage_intervals_for_read, clean_chr_name, \
    get_strand, CAGEReads, RNAseqReads

import multiprocessing
import threading

BUFFER_SIZE = 50000000

class ProcessSafeOPStream( object ):
    def __init__( self, writeable_obj ):
        self.writeable_obj = writeable_obj
        self.lock = multiprocessing.Lock()
        self.name = self.writeable_obj.name
        return
    
    def write( self, data ):
        self.lock.acquire()
        self.writeable_obj.write( data )
        self.writeable_obj.flush()
        self.lock.release()
        return
    
    def close( self ):
        self.writeable_obj.close()

def update_buffer_array_from_rnaseq_read( buffer_array, 
                                          buffer_offset, 
                                          read, strand,
                                          reverse_read_strand,
                                          pairs_are_opp_strand ):
                                          
    """populate buffer with histogram of contiguous read regions

    """
    for start, stop in iter_coverage_intervals_for_read(read):
        rd_strand = get_strand( read, 
                                reverse_read_strand=reverse_read_strand, 
                                pairs_are_opp_strand=pairs_are_opp_strand )
        if rd_strand != strand: continue
        buffer_array[(start-buffer_offset):(stop-buffer_offset) + 1] += 1
    
    return

def update_buffer_array_from_rnaseq_read_generator( reads ):
    def update_buffer_array_from_read(buffer_array, buffer_offset,
                                      strand, reverse_read_strand, read):
        return update_buffer_array_from_rnaseq_read( 
            buffer_array, buffer_offset, 
            read, strand, reverse_read_strand,
            reads.pairs_are_opp_strand )
    
    return update_buffer_array_from_read

def update_buffer_array_from_polya_read(
        buffer_array, buffer_offset, strand, reverse_read_strand, read):
    rd_strand = '+' if read.is_reverse else '-'
    if reverse_read_strand:
        rd_strand = '-' if rd_strand == '+' else '+'
    
    # skip reads that dont match the filtering criterion
    if not read.is_read1: return
    if strand != rd_strand: return
    
    # determine which pos of the read corresponds to the 
    # poly(a) site
    if rd_strand == '+': pos = read.aend
    else: pos = read.pos
    
    # find the statmap posterior probabiliy, if available
    res = [ val for key, val in read.tags if key == 'XP' ]
    post_prb = 1.0 if len(res) == 0 else res[0]
    
    # update the array
    buffer_array[pos-buffer_offset] += post_prb
    
    return


def update_buffer_array_from_CAGE_read(
        buffer_array, buffer_offset, strand, reverse_read_strand, read):
    rd_strand = '+' if read.is_reverse else '-'
    if reverse_read_strand:
        rd_strand = '-' if rd_strand == '+' else '+'
    
    # skip reads that dont match the filtering criterion
    if strand != rd_strand: return
    
    # determine which pos of the read corresponds to the 
    # poly(a) site
    if rd_strand == '+': pos = read.pos
    else: pos = read.aend
    
    # find the statmap posterior probabiliy, if available
    res = [ val for key, val in read.tags if key == 'XP' ]
    try: post_prb = float(res[0])
    except Exception: post_prb = 1.0
    
    # update the array
    buffer_array[pos-buffer_offset] += post_prb
    
    return

def populate_cvg_array_for_contig( 
        merged_ofp, reads_fname, chrm, chrm_length, strand, 
        reverse_read_strand, update_buffer_array_from_read ):
    if VERBOSE: print "Starting ", chrm, strand
    
    # open the reads file - we pass in a filename and re-open to make 
    # this multi-process safe
    reads = pysam.Samfile( reads_fname, "rb" )
    
    # open a tempory file to write this to
    ofp = tempfile.NamedTemporaryFile(delete=False)
    
    # we make this two times the buffer size so that we can do a memmove rather 
    # then a new alloc when we write the buffer to disk. Note that the buffer 
    # size must be larger than the gap between reads for this to work. For a 
    # normal RNA experiment, 1MB should be more than enough, so we set it to 
    # 5MB to be extra safe.
    buffer_array = numpy.zeros( BUFFER_SIZE*2 )
    # stores after how many bases in the current contig ( chrm, strand ) the 
    # current buffer starts.
    buffer_offset = None
    for read in reads.fetch(chrm):
        # optimize the writing process by skipping regions that start before 
        # the first read. We subtract an additional megabase to account for 
        # weird potential read offsets
        if buffer_offset == None: 
            buffer_offset = max(0, read.pos-1e6)
        
        # if this read extends past the current buffer, then we need to write
        # it out to disk and move the unwritten portion tot he start of the 
        # buffer
        if read.pos > buffer_offset + BUFFER_SIZE:
            write_array_to_opstream( 
                ofp, buffer_array[:BUFFER_SIZE], 
                buffer_offset, chrm, chrm_length)
            
            # move the unwritten portion to the start of the buffer,
            # and zero out the end
            buffer_array[:BUFFER_SIZE] = buffer_array[BUFFER_SIZE:]
            buffer_array[BUFFER_SIZE:] = 0
            buffer_offset += BUFFER_SIZE

        update_buffer_array_from_read( 
            buffer_array, buffer_offset, strand, reverse_read_strand, read )
    
    #to make sure the rest of the buffer is stored on disk
    write_array_to_opstream( ofp, buffer_array[:BUFFER_SIZE], 
                             buffer_offset, chrm, chrm_length)
    reads.close()

    ofp.seek(0)
    merged_ofp.write( ofp.read() )
    ofp.close()

    if VERBOSE: print "Finished ", chrm, strand
    
    return


def write_array_to_opstream(ofp, buffer, buff_start, chrm, chrm_length ):
    """write buffer to disk, buff_start determines the start of buffer in 
       genomic coordinates.
    """
    chrm = clean_chr_name( chrm )
    
    prev_pos = 0
    prev_val = buffer[0]
    for pos, val in enumerate(buffer[1:]):
        # make sure this doesn't extend past the end of the chromosome
        # bedGraphs are 0-based, so use chrm_length-1
        if buff_start+pos+1 >= chrm_length:
            pos = chrm_length-buff_start-1
            break
        if val != prev_val:
            if prev_val > 1e-12:
                line = "chr%s\t%i\t%i\t%.2f" % (
                    chrm, buff_start+prev_pos, buff_start+pos+1, prev_val )
                ofp.write(line+"\n")
            prev_pos, prev_val = pos+1, val
    
    if prev_val > 1e-12:
        line = "chr%s\t%i\t%i\t%.2f" % (
            chrm, buff_start+prev_pos, buff_start+pos+1, prev_val )
        ofp.write(line+"\n")
    
    return


def build_chrm_sizes_file(reads):
    chrm_sizes_file = tempfile.NamedTemporaryFile(delete=True)
    # find the chrm names and their associated lengths
    chrm_lengths = zip(reads.references, reads.lengths)
    #write out the chromosomes and its corrosponding size to disk
    for chrm, chrm_length in chrm_lengths:
        chrm_sizes_file.write(chrm + "   " + str(chrm_length) +"\n")
    chrm_sizes_file.flush()
    
    return chrm_sizes_file

def generate_wiggle(reads, ofps, update_buffer_array_from_read, 
                    num_threads=1, reverse_read_strand=None ):
    all_args = []
    
    for chrm_length, chrm  in sorted(izip(reads.lengths, reads.references)):
        strands = ['+', '-'] if len(ofps) == 2 else [None,]
        for strand in strands:
            ofp = ofps[strand]
            all_args.append((ofp, reads.filename, chrm, chrm_length,
                             strand, reverse_read_strand,
                             update_buffer_array_from_read))
    
    if num_threads == 1:
        for args in reversed(all_args):
            populate_cvg_array_for_contig( *args )
    else:
        ps = [None]*num_threads
        while len( all_args ) > 0:
            for i in xrange(num_threads):
                if ps[i] == None or not ps[i].is_alive():
                    ps[i] = multiprocessing.Process( 
                        target=populate_cvg_array_for_contig, 
                        args=all_args.pop() )
                    ps[i].start()
                    break
            time.sleep( 0.1 )

        for p in ps:
            if p != None: p.join()
    
    for fp in ofps.values(): fp.close()
    reads.close()
    
    return

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(
        description='Get coverage bedgraphs from aligned reads.')
    parser.add_argument( '--mapped-reads-fname', required=True,
                         help='BAM or SAM file(s) containing the mapped reads.')
    parser.add_argument( '--out-fname-prefix', '-o', required=True, 
                         help='Output file(s) will be bigWig')
    parser.add_argument( '--assay', '-a', required=True, 
                         help='The assay type [(r)naseq, (c)age, (p)olya]')    
    parser.add_argument( '--bigwig', '-b', default=False, action='store_true', 
                         help='Build a bigwig instead of bedgraph.')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', 
                         help='Whether or not to print status information.')
    parser.add_argument( '--threads', '-t', default=1, type=int,
                         help='The number of threads to run.')
    
    parser.add_argument( '--reverse-read-strand', '-r', default=False, action='store_true',
                         help='Whether or not to reverse the strand of the read. default: False')
        
    args = parser.parse_args()
    global VERBOSE
    VERBOSE = args.verbose
    
    assay = {'c': 'cage', 'r': 'rnaseq', 'p': 'polya'}[args.assay.lower()[0]]
    if assay not in ('cage', 'rnaseq', 'polya'):
        raise ValueError, "Unrecongized assay (%s)" % args.assay
    
    return assay, args.mapped_reads_fname, args.out_fname_prefix, args.bigwig, \
        args.reverse_read_strand, args.threads

def build_bigwig_from_bedgraph(bedgraph_fp, chrm_sizes_file, op_fname):
    with tempfile.NamedTemporaryFile(delete=True) as sorted_ofp:
        if VERBOSE: print "Sorting ", bedgraph_fp.name
        subprocess.call( 
            ["sort -k1,1 -k2,2n " + bedgraph_fp.name,], 
            stdout=sorted_ofp, shell=True )
        sorted_ofp.flush()
        
        if VERBOSE: print "Building wig for", bedgraph_fp.name
        subprocess.check_call( [ "bedGraphToBigWig", 
                                 sorted_ofp.name, 
                                 chrm_sizes_file.name, 
                                 op_fname ] )
    return

def main():
    ( assay, reads_fname, op_prefix, build_bigwig, 
      reverse_read_strand, num_threads ) = parse_arguments()
    
    # initialize the assay specific options
    if assay == 'cage':
        reads = CAGEReads( reads_fname, "rb" )
        # the read strand reversal is done later, so set this to False
        reads.init(reverse_read_strand=False)
        update_buffer_array_from_read = update_buffer_array_from_CAGE_read
        stranded = True
    elif assay == 'polya':
        reads = pysam.Samfile( reads_fname, "rb" )
        update_buffer_array_from_read = update_buffer_array_from_polya_read
        stranded = True
    elif assay == 'rnaseq':
        reads = RNAseqReads( reads_fname, "rb" )
        # the read strand reversal is done later, so set this to False
        reads.init(reverse_read_strand=False)
        update_buffer_array_from_read = \
            update_buffer_array_from_rnaseq_read_generator(reads)
        stranded = reads.reads_are_stranded
    else:
        raise ValueError, "Unrecognized assay: '%s'" % assay
    
    # if we want to build a bigwig, make sure that the script is on the path
    if build_bigwig:
        try: 
            subprocess.check_call(["which", "bedGraphToBigWig"], stdout=None)
        except subprocess.CalledProcessError:
            raise ValueError, "bedGraphToBigWig does not exist on $PATH. " + \
                "You can still build a bedGraph by removing the --bigwig(-b) option."        
    
    # Open the output files
    if stranded:
        ofps = { '+' : ProcessSafeOPStream(
                open(op_prefix+".plus.bedgraph","w")), 
                 '-' : ProcessSafeOPStream(
                open(op_prefix+".minus.bedgraph", "w"))
               }
    else:
        ofps = { None: ProcessSafeOPStream(open(op_prefix+".bedgraph", "w")) }

    # write the bedgraph header information
    if not build_bigwig:
        for key, fp in ofps.iteritems():
            strand_str = "" if key == None else {'+': '.plus', '-': '.minus'}[key]
            fp.write( "track name=%s.%s type=bedGraph\n" \
                          % ( os.path.basename(op_prefix), strand_str ) )
    
    
    generate_wiggle( reads, ofps, update_buffer_array_from_read, num_threads,
                     reverse_read_strand=reverse_read_strand )
    
    # finally, if we are building a bigwig, build it, and then remove the bedgraph files
    if build_bigwig:
        # build the chrm sizes file.
        with build_chrm_sizes_file(reads) as chrm_sizes_file:        
            threads = []
            for strand, bedgraph_fp in ofps.iteritems():
                strand_str = "" if strand == None else ( 
                    {'+': '.plus', '-': '.minus'}[strand] )
                op_fname = op_prefix + strand_str + ".bw"

                t = threading.Thread( 
                    target=build_bigwig_from_bedgraph, 
                    args=(bedgraph_fp, chrm_sizes_file, op_fname) )
                t.start()
                threads.append( t )

            for t in threads:
                t.join()
        
        chrm_sizes_file.close()
    else:
        for fp in ofps.values(): fp.close()

if __name__ == "__main__":
    main()
