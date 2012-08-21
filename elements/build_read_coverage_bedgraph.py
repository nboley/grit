# Copyright (c) 2011-2012 Nathan Boley

import sys, os
import pysam
import numpy

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../file_types/" ) )
from reads import iter_coverage_regions_for_read, clean_chr_name, \
    read_pairs_are_on_same_strand

import multiprocessing
from multiprocessing import Process
num_threads = 1

class threadsafeFile( file ):
    def __init__( self, *args ):
        self.lock = multiprocessing.Lock()
        file.__init__( self, *args )
    
    def write( self, data ):
        self.lock.acquire()
        file.write( self, data )
        file.flush( self )
        self.lock.release()

def populate_cvg_arrays( reads_fname, chrm, separate_read_pairs, 
                         reverse_read_strand, pairs_are_opp_strand ):
                     
    """Get all of the junctions represented in a reads object
    """

    reads = pysam.Samfile( reads_fname, "rb" )

    # get the refernece length. If there is no such reference, return
    ref_index = reads.references.index( chrm )
    if ref_index == -1:
        reads.close()
        return
    chrm_len = reads.lengths[ ref_index ]
    
    
    # try to get one read. If we cant, then skip this chromosome
    try:
        next( iter(reads.fetch(chrm)) )
    except StopIteration:
        reads.close()
        return

    read1_cov = { '+': numpy.zeros(chrm_len), '-': numpy.zeros(chrm_len) }
    if separate_read_pairs:
        read2_cov = { '+': numpy.zeros(chrm_len), '-': numpy.zeros(chrm_len) }
    else:
        read2_cov = None
    
    for i, read in enumerate(reads.fetch(chrm)):
        # skip reads that aren't mapped in pair
        if not read.is_proper_pair:
            continue
        
        for chrm, strand, start, stop in iter_coverage_regions_for_read( 
                read, reads, reverse_read_strand, pairs_are_opp_strand ):
            if separate_read_pairs and not read.is_read1:
                arr = read2_cov[ strand ]
            else:
                arr = read1_cov[ strand ]
            arr[start:stop+1] += 1.0
    
    reads.close()
    
    return read1_cov, read2_cov

def build_bedgraph_lines_from_array( ofp, array, chrm ):
    lines = []
    prev_pos = 0
    prev_val = array[0]
    for pos, val in enumerate(array[1:]):
        if val != prev_val:
            if prev_val > 1e-12:
                lines.append( "chr%s\t%i\t%i\t%.2f" % (
                        chrm, prev_pos, pos+1, prev_val) )
            prev_pos = pos+1
            prev_val = val
    
    if prev_val > 1e-12:
        lines.append( "chr%s\t%i\t%i\t%.2f" % (chrm, prev_pos, pos+2, prev_val) )
    ofp.write( "\n".join( lines ) + "\n" )
    return


def populate_wiggles_worker( reads_fname, chrm, rd1_ofps, rd2_ofps, 
                             reverse_read_strand, pairs_are_opp_strand, 
                             separate_read_pairs ):
    res = populate_cvg_arrays( reads_fname, chrm, 
                               separate_read_pairs, 
                               reverse_read_strand,
                               pairs_are_opp_strand )
    
    # if we couldn't get a single read, then continue
    if res == None: return
    else: read1_arr, read2_arr = res

    bedgraph_lines = build_bedgraph_lines_from_array( 
        rd1_ofps['+'], read1_arr['+'], clean_chr_name(chrm) )
    bedgraph_lines = build_bedgraph_lines_from_array( 
        rd1_ofps['-'], read1_arr['-'], clean_chr_name(chrm) )

    if rd2_ofps != None:
        bedgraph_lines = build_bedgraph_lines_from_array( 
            rd2_ofps['+'], read2_arr['+'], clean_chr_name(chrm) )
        bedgraph_lines = build_bedgraph_lines_from_array( 
            rd2_ofps['-'], read2_arr['-'], clean_chr_name(chrm) )
    
    return

def populate_wiggle( reads_fname, chrms, rd1_ofps, rd2_ofps, 
                     reverse_read_strand, pairs_are_opp_strand ):
    """Get all of the junctions represented in a reads object
    """
    separate_read_pairs = ( rd2_ofps != None )
    
    all_args = []
    for chrm in chrms:
        all_args.append( ( reads_fname, chrm, rd1_ofps, rd2_ofps, 
                           reverse_read_strand, pairs_are_opp_strand,
                           separate_read_pairs ) )

    if num_threads == 1:
        for args in all_args:
            populate_wiggles_worker( *args )
    else:
        from time import sleep
        
        ps = []
        for thread_id, args in zip( range(num_threads), all_args ):
            p = Process(target=populate_wiggles_worker, name=args[1], args=args)
            p.start()
            ps.append( p )
        
        del all_args[:num_threads]
        
        while len(all_args) > 0:
            sleep( 0.5 )
            for p_i, p in enumerate(ps):
                if len( all_args ) == 0:
                    break
                if not p.is_alive():
                    ps[p_i].join()
                    args = all_args.pop()
                    ps[p_i] = Process(target=populate_wiggles_worker, 
                                      name=args[1], args=args)
                    ps[p_i].start()
                
        for p in ps:
            p.join()
        
        assert len( all_args ) == 0
    
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
    parser.add_argument( '--threads', '-t', default=1, type=int,
        help='The numnber of threads to run.')
    args = parser.parse_args()
    
    global num_threads
    num_threads = args.threads
    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.mapped_reads_fname, args.chrm_sizes, args.out_fname_prefix, \
        args.reverse_read_strand, args.merge_read_ends

def main():
    reads_fname, chrm_sizes, op_prefix, reverse_read_strand, merge_read_ends \
        = parse_arguments()

    track_prefix = os.path.basename( op_prefix )
    
    # make sure that we can open the reads object
    reads = pysam.Samfile( reads_fname, "rb" )
    chrm_names = reads.references
    pairs_are_opp_strand = not read_pairs_are_on_same_strand( reads )
    reads.close()
    
    # write the wiggle to disk
    if merge_read_ends:
        rd1_fps = { 
            '+': threadsafeFile("{0}.{1}.bedGraph".format( op_prefix, "plus"), "w"),
            '-': threadsafeFile("{0}.{1}.bedGraph".format( op_prefix, "minus"), "w")
            }
        for key, val in rd1_fps.iteritems():
            strand_str = 'plus' if key == '+' else 'minus'
            val.write("track type=bedGraph name={0}_{1}\n".format(track_prefix, strand_str))

        rd2_fps = None
    else:
        rd1_fps = {
            '+': threadsafeFile("{0}.{1}.{2}.bedGraph".format( op_prefix, "rd1", "plus"), "w"),
            '-': threadsafeFile("{0}.{1}.{2}.bedGraph".format( op_prefix, "rd1", "minus"), "w")
            }
        for key, val in rd1_fps.iteritems():
            strand_str = 'plus' if key == '+' else 'minus'
            val.write("track type=bedGraph name={0}_rd1_{1}\n".format(track_prefix, strand_str))

        rd2_fps = {
            '+': threadsafeFile("{0}.{1}.{2}.bedGraph".format( op_prefix, "rd2", "plus"), "w"),
            '-': threadsafeFile("{0}.{1}.{2}.bedGraph".format( op_prefix, "rd2", "minus"), "w")
            }
        for key, val in rd2_fps.iteritems():
            strand_str = 'plus' if key == '+' else 'minus'
            val.write("track type=bedGraph name={0}_rd2_{1}\n".format(track_prefix, strand_str))

        
    
    # populate the wiggle from the bam file
    populate_wiggle( reads_fname, chrm_names, 
                     rd1_fps, rd2_fps,
                     reverse_read_strand, pairs_are_opp_strand )
    return

if __name__ == "__main__":
    main()
