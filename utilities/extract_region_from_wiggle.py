#!/usr/bin/python

# import python mods
import os 
import sys

sys.path.append( os.path.join(os.path.dirname(__file__), "..", "file_types") )
from wiggle import Wiggle, GenomicInterval

VERBOSE = False

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description='Merge counts of many wiggle files into one wig file.')
    
    parser.add_argument( 'filterRegion', \
        help="A filter region of the form chr:strand:start-stop. Use '.' for both strands.")
    parser.add_argument( 'wig', type=file, 
                             help='Wiggle file to extract region from.')
    parser.add_argument( 'chrmSizes', type=file, 
                             help='A chromosome sizes file.')
    parser.add_argument( 'outFnamePrefix', help='Output filename - +|-.wig will be appended.' )
    args = parser.parse_args()
    
    # set flag args
    global VERBOSE
    
    chrm, strand, locs = args.filterRegion.strip().split(":")
    assert strand in ".-+"
    start, stop = locs.split("-")
    chrm, start, stop = chrm.strip(), int(start), int(stop)
    filter_region = GenomicInterval( chrm, strand, start, stop )
    
    return args.wig, args.outFnamePrefix, args.chrmSizes, filter_region

def main():
    wig, out_fname_prefix, chrm_sizes_fp, filter_region = parse_arguments()
    
    op_wiggle = Wiggle( chrm_sizes_fp )
    op_wiggle.load_data_from_fp( wig )

    if VERBOSE:
        print 'Writing merged wiggle...'
    
    op_wiggle.write_wiggles( out_fname_prefix + ".+.wig",  
                             out_fname_prefix + ".-.wig",  
                             ignore_zeros=True, 
                             filter_region = filter_region )
    
if __name__ == '__main__':
    main()
