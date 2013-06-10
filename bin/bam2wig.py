import sys, os
import pysam
import numpy
import shutil
import subprocess
import tempfile
import time
from itertools import izip

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      ".." ) )
from grit.files.reads import iter_coverage_regions_for_read, clean_chr_name, \
    read_pairs_are_on_same_strand, get_strand

import multiprocessing

num_threads = 1
buffer_size = 5000000
bedGraphToBigWig_script = "/usr/local/bin/bedGraphToBigWig"


class ProcessSafeOPStream( object ):
    def __init__( self, writeable_obj ):
        self.writeable_obj = writeable_obj
        self.lock = multiprocessing.Lock()
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
                                          reads, 
                                          reverse_read_strand, 
                                          pairs_are_opp_strand ):
                                          
    """populate buffer with histogram of contiguous read regions

    """
    for chrm, rd_strand, start, stop in iter_coverage_regions_for_read(
            read, reads, reverse_read_strand, pairs_are_opp_strand):
        if rd_strand != strand: continue
        buffer_array[(start-buffer_offset):(stop-buffer_offset) + 1] += 1
    return

def update_buffer_array_from_rnaseq_read_generator( 
        reads, reverse_read_strand, pairs_are_opp_strand ):
    def update_buffer_array_from_read(buffer_array, buffer_offset, strand, read):
        return update_buffer_array_from_rnaseq_read( 
            buffer_array, buffer_offset, 
            read, strand, reads, 
            reverse_read_strand, pairs_are_opp_strand )
    
    return update_buffer_array_from_read

def update_buffer_array_from_polya_read(
        buffer_array, buffer_offset, strand, read):
    rd_strand = '+' if read.is_reverse else '-'
    
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

urvish_crap = """
    # check if the chromosome exists in the references
    ref_index = reads.references.index( chrm )
    if ref_index == -1:
        reads.close()
        return

    # try to get one read. If we cant, then skip this chromosome
    try:
        next( iter(reads.fetch(chrm)) )
    except StopIteration:
        reads.close()
        return


    #create file for output
    if strand:
        ofp = file(os.path.abspath("/tmp/{0}.{1}.bedGraph".format(clean_chr_name(chrm), strand)), "w")
    else:
        ofp = file(os.path.abspath("/tmp/{0}.bedGraph".format(clean_chr_name(chrm))), "w")




        # skip reads with the incorrect strand
        rd_strand = get_strand(read, reverse_read_strand, pairs_are_opp_strand )
        if strand != None and rd_strand != strand:                 
            continue            
        
        # skip reads that aren't mapped in pair
        if not read.is_proper_pair:
            continue
        

"""

def populate_cvg_array_for_contig( 
        merged_ofp, reads_fname, chrm, strand, update_buffer_array_from_read ):
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
    buffer_array = numpy.zeros( buffer_size*2 )
    # stores after how many bases in the current contig ( chrm, strand ) the 
    # current buffer starts.
    buffer_offset = None
    for read in reads.fetch(chrm):
        # optimize the writing process by skipping regions that start before 
        # the first read
        if buffer_offset == None: buffer_offset = read.pos
        
        # if this read extends past the current buffer, then we need to write
        # it out to disk and move the unwritten portion tot he start of the 
        # buffer
        if read.pos > buffer_offset + buffer_size:
            write_array_to_opstream( 
                ofp, buffer_array[:buffer_size], 
                buffer_offset, chrm)
            
            # move the unwritten portion to the start of the buffer,
            # and zero out the end
            buffer_array[:buffer_size] = buffer_array[buffer_size:]
            buffer_array[buffer_size:] = 0
            buffer_offset += buffer_size

        update_buffer_array_from_read( 
            buffer_array, buffer_offset, strand, read )
    
    #to make sure the rest of the buffer is stored on disk
    write_array_to_opstream( ofp, buffer_array[:buffer_size], 
                             buffer_offset, chrm)
    reads.close()

    ofp.seek(0)
    merged_ofp.write( ofp.read() )
    ofp.close()

    if VERBOSE: print "Finished ", chrm, strand
    
    return



def write_array_to_opstream(ofp, buffer, buff_start, chrm ):
    """write buffer to disk, buff_start determines the start of buffer in 
       genomic coordinates.
    """
    chrm = clean_chr_name( chrm )
    
    prev_pos = 0
    prev_val = buffer[0]
    for pos, val in enumerate(buffer[1:]):
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
    chrm_sizes_file = tempfile.NamedTemporaryFile(delete=False)
    chrm_names = reads.references
    chrm_lengths = zip(chrm_names, reads.lengths)
    #write out the chromosomes and its corrosponding size to disk
    for chrm, chrm_length in chrm_lengths:
        chrm_sizes_file.write(chrm + "   " + str(chrm_length) +"\n")
    
    return chrm_sizes_file

def generate_wiggle(reads_fname, op_prefix, assay, stranded=True ):
    reads = pysam.Samfile( reads_fname, "rb" )
    
    if assay == 'polya':
        update_buffer_array_from_read = update_buffer_array_from_polya_read
        stranded = True
    elif assay == 'rnaseq':
        update_buffer_array_from_read = \
            update_buffer_array_from_rnaseq_read_generator(reads, False, False)
    else:
        raise ValueError, "Unrecognized assay: '%s'" % assay
    
    if stranded:
        ofps = { '+' : ProcessSafeOPStream(
                open(op_prefix+".plus.bedgraph","w")), 
                 '-' : ProcessSafeOPStream(
                open(op_prefix+".minus.bedgraph", "w"))
               }
    else:
        ofps = { None: ProcessSafeOPStream(open(op_prefix+".bedgraph", "w")) }
    for key, fp in ofps.iteritems():
        strand_str = "" if key == None else {'+': '.plus', '-': '.minus'}[key]
        fp.write( "track name=%s.%s type=bedGraph\n" \
                      % ( os.path.basename(op_prefix), strand_str ) )
    
    all_args = []
    
    for chrm_length, chrm  in sorted(izip(reads.lengths, reads.references)):
        strands = ['+', '-'] if stranded else [None,]
        for strand in strands:
            ofp = ofps[strand]
            all_args.append((ofp, reads.filename, chrm, strand, 
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
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', 
                         help='Whether or not to print status information.')
    parser.add_argument( '--threads', '-t', default=10, type=int,
                         help='The number of threads to run.')

    parser.add_argument( '--reverse-read-strand', '-r', type=bool,
                         help='Whether or not to reverse the strand of the read. default: infer from junction reads')
    parser.add_argument('--stranded', type=bool,
                        help="true if reads are stranded")
    parser.add_argument('--reads-opp-strand', type=bool, 
                        help='True if reads are opp strand for each BAM file')
    parser.add_argument('--buffer-size', '-b', default=5000000, type=int,
                        help='The amount of memory(in base pairs) to use before flushing to disk')
    parser.add_argument('--force', '-f', default=False, action='store_true', 
                        help='force recreate of bigwig even if it exists')

    args = parser.parse_args()
    global buffer_size
    buffer_size = args.buffer_size
    global num_threads
    num_threads = args.threads    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.mapped_reads_fname, args.out_fname_prefix, \
        args.reverse_read_strand, args.reads_opp_strand, args.stranded

def main():
    reads_fnames, op_prefixes, reverse_read_strand, reads_opp_strand, stranded\
        = parse_arguments()
    
    generate_wiggle( reads_fnames, op_prefixes, "rnaseq" )

if __name__ == "__main__":
    main()
