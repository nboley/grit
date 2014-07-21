import sys, os
import pysam
import numpy
import shutil
import subprocess
import tempfile
import time
from itertools import izip

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), ".." ) )
from grit.files.reads import clean_chr_name, fix_chrm_name_for_ucsc, \
    CAGEReads, RAMPAGEReads, RNAseqReads, PolyAReads
from grit.lib.multiprocessing_utils import ProcessSafeOPStream

import multiprocessing
import threading

# if we choose the --ucsc option, then replace thsi function
# with fix_chrm_name_for_ucsc
def fix_chrm_name(name):
    return name

BUFFER_SIZE = 50000000

def populate_cvg_array_for_contig( 
        merged_ofp, reads, chrm, chrm_length, strand ):
    if VERBOSE: print "Starting ", chrm, strand
    
    # re-open the reads to make this multi-process safe
    reads.reload()
    
    # open a tempory file to write this to
    with tempfile.NamedTemporaryFile(delete=True) as ofp:
        # only find blocks of BUFFER_SIZE - to avoid running out of memory
        for block_index in xrange(int(chrm_length/BUFFER_SIZE)+1):
            buffer_array = reads.build_read_coverage_array( 
                chrm, strand, 
                block_index*BUFFER_SIZE, 
                (block_index+1)*BUFFER_SIZE )
            write_array_to_opstream( 
                ofp, buffer_array, block_index*BUFFER_SIZE, 
                chrm, chrm_length, strand)
        
        ofp.seek(0)
        merged_ofp.write( ofp.read() )
    
    if VERBOSE: print "Finished ", chrm, strand
    
    return


def write_array_to_opstream(ofp, buffer, buff_start, 
                            chrm, chrm_length, strand ):
    """write buffer to disk, buff_start determines the start of buffer in 
       genomic coordinates.
    """
    chrm = fix_chrm_name( clean_chr_name( chrm ) )
    
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
                write_val = -prev_val if strand == '-' else prev_val
                line = "%s\t%i\t%i\t%.2f" % (
                    chrm, buff_start+prev_pos, buff_start+pos+1, write_val )
                ofp.write(line+"\n")
            prev_pos, prev_val = pos+1, val
    
    if prev_val > 1e-12:
        write_val = -prev_val if strand == '-' else prev_val
        line = "%s\t%i\t%i\t%.2f" % (
            chrm, buff_start+prev_pos, buff_start+pos+1, write_val )
        ofp.write(line+"\n")
    
    return


def build_chrm_sizes_file(reads):
    chrm_sizes_file = tempfile.NamedTemporaryFile(delete=True)
    # find the chrm names and their associated lengths
    chrm_lengths = zip(reads.references, reads.lengths)
    #write out the chromosomes and its corrosponding size to disk
    for chrm, chrm_length in chrm_lengths:
        chrm_sizes_file.write(fix_chrm_name(chrm) + "   " + str(chrm_length) +"\n")
    chrm_sizes_file.flush()
    
    return chrm_sizes_file

def generate_wiggle(reads, ofps, num_threads=1, contig=None ):
    all_args = []
    for chrm_length, chrm  in sorted(izip(reads.lengths, reads.references)):
        strands = ['+', '-'] if len(ofps) == 2 else [None,]
        # skip regions not in the specified contig, if requested 
        if contig != None and clean_chr_name(chrm) != clean_chr_name(contig): 
            continue
        for strand in strands:
            ofp = ofps[strand]
            all_args.append((ofp, reads, chrm, chrm_length, strand ))
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
    
    return

def parse_arguments():
    allowed_assays = ['cage', 'rampage', 'rnaseq', 'polya']
    
    import argparse
    parser = argparse.ArgumentParser(
        description='Get coverage bedgraphs from aligned reads.')
    parser.add_argument( '--mapped-reads-fname', required=True,
                         help='BAM or SAM file(s) containing the mapped reads.')
    parser.add_argument( '--out-fname-prefix', '-o', 
                         help='Output file(s) will be bigWig')
    parser.add_argument( '--assay', '-a', required=True, 
                         choices=allowed_assays, help='The assay type')
    parser.add_argument( '--bigwig', '-b', default=False, action='store_true', 
                         help='Build a bigwig instead of bedgraph.')
    parser.add_argument( '--ucsc', default=False, action='store_true', 
                         help='Format the contig names to work with the UCSC genome browser.')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', 
                         help='Whether or not to print status information.')
    parser.add_argument( '--threads', '-t', default=1, type=int,
                         help='The number of threads to run.')
    
    parser.add_argument( '--region', 
        help='Only use the specified region ( currently only accepts a contig name ).')
    parser.add_argument( '--reverse-read-strand', '-r', default=False, action='store_true',
        help='Whether or not to reverse the strand of the read ( only sensible for RNAseq reads ). default: False')
    parser.add_argument( '--read-filter', default=None, choices=['1','2'],
        help='Filter paired end reads to only accept this read pair (ie uses the is_read1 pysam attribute)')

        
    args = parser.parse_args()
    global VERBOSE
    VERBOSE = args.verbose
    
    global fix_chrm_name
    if args.ucsc: fix_chrm_name = fix_chrm_name_for_ucsc
    
    assert args.read_filter in ( '1', '2', None )
    read_filter = int(args.read_filter) if args.read_filter != None else None
    
    if args.assay not in allowed_assays:
        raise ValueError, "Unrecongized assay (%s)" % args.assay
    
    region = args.region
    if region != None:
        if ':' in region or '-' in region:
            assert False, "Invalid contig name: %s" % region
    
    # if an output prefix isn't provided, then use the bam filename prefix
    if args.out_fname_prefix == None:
        fname_data = args.mapped_reads_fname.split(".")
        # remove bam and sorted suffixes
        while fname_data[-1] in ('bam', 'sorted'):
            del fname_data[-1]
        args.out_fname_prefix = ".".join(fname_data)
    
    return ( args.assay, args.mapped_reads_fname, args.out_fname_prefix, 
             args.bigwig, args.reverse_read_strand, read_filter, 
             args.region, args.threads )
        

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
      reverse_read_strand, read_filter, region, num_threads) = parse_arguments()
    
    # initialize the assay specific options
    if assay == 'cage':
        reads = CAGEReads( reads_fname, "rb" )
        reads.init(reverse_read_strand=False)
        stranded = True
        assert not reverse_read_strand
    elif assay == 'rampage':
        reads = RAMPAGEReads( reads_fname, "rb" )
        reads.init(reverse_read_strand=True)
        stranded = True
        assert not reverse_read_strand
    elif assay == 'polya':
        reads = PolyAReads( reads_fname, "rb" )
        reads.init(reverse_read_strand=reverse_read_strand, pairs_are_opp_strand=True)
        stranded = True
    elif assay == 'rnaseq':
        reads = RNAseqReads( reads_fname, "rb" )
        # the read strand reversal is done later, so set this to False
        reads.init(reverse_read_strand=reverse_read_strand)
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
    
    
    generate_wiggle( reads, ofps, num_threads, region )
    
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
    
    # close the reads files
    reads.close()


if __name__ == "__main__":
    main()
