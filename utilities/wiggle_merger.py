#!/usr/bin/python

# Copyright (c) 2011-2012 Nathan Boley

# import python mods
import os 
import sys
import tempfile
import subprocess

sys.path.append( os.path.join(os.path.dirname(__file__), "..", "file_types") )
from wiggle import Wiggle, GenomicInterval

VERBOSE = False

def create_sorted_temp_file( fname ):
    return ofp

def merge_sorted_temp_files( fps ):
    if VERBOSE:
        print >> sys.stderr, "Merging %i bedgraphs" % len( fps )
    ofp = tempfile.NamedTemporaryFile("w+")
    ifnames = [ fp.name for fp in fps ]
    args = [ "unionBedGraphs", "-i" ]
    args.extend( ifnames )
    rc = subprocess.check_call( args, stdout=ofp, stderr=subprocess.PIPE )
    ofp.flush()
    ofp.seek(0)
    return ofp


def fast_merge_bedgraphs( fps ):
    if VERBOSE:
        print >> sys.stderr, "Sorting bedgraphs"
    # first, sort each bedgraph into a temporary file, and close the orig fps
    sorted_temp_fps = []
    calls = []
    for fp in fps:
        ofp = tempfile.NamedTemporaryFile("w+")
        sorted_temp_fps.append( ofp )
        call = subprocess.Popen( [ "sortBed", "-i", "%s" % fp.name ], 
                               stdout=ofp, stderr=subprocess.PIPE ) 
        calls.append( call )
        
    for call in calls:
        call.wait()
    
    for sorted_fp in sorted_temp_fps:
        sorted_fp.flush()
        sorted_fp.seek( 0 )
        
    for fp in fps:
        fp.close()
    
    # now, merge them together
    union_bed_fp = merge_sorted_temp_files( sorted_temp_fps )
    
    # close the sorted files 
    for fp in sorted_temp_fps:
        fp.close()
        
    if VERBOSE:
        print >> sys.stderr, "Summing bedgraphs"
    # finally, merge the lines together
    for line in union_bed_fp:
        data = line.split()
        data[3] = str( sum( map( float, data[3:] ) ) )
        del data[4:]
        print "\t".join( data )
    
    union_bed_fp.close()
    
    return

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Merge counts of many wiggle files into one wig file.')
    parser.add_argument( 'wigs', type=file, nargs='+', \
                             help='Wiggle files to merge into one wiggle.')
    parser.add_argument( '--chrm-sizes', type=file, # required=True, \
                             help='A chromosome sizes file.')
    parser.add_argument( '--out-fname-prefix', default="merged_wiggle", \
                             help='Output filenames will be PREFIX.plus.wig and ' \
                             + "PREFIX.minus.wig. Default: %(default)s")
    parser.add_argument( '--track-name-prefix', default="merged_wig", \
                             help='The output track names will be PREFIX_plus and PREFIX_minus.')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', \
                             help='Whether or not to print status information.')
    parser.add_argument( '--filter-region', '-f', \
                             help='A filter region of the form chr:start-stop.')
    args = parser.parse_args()
    
    # set flag args
    global VERBOSE
    VERBOSE = args.verbose
    
    if args.filter_region != None:
        chrm, locs = args.filter_region.strip().split(":")
        start, stop = locs.split("-")
        chrm, start, stop = chrm.strip(), int(start), int(stop)
        filter_region = GenomicInterval( chrm, '.', start, stop )
    else:
        filter_region = None
    
    return args.wigs, args.out_fname_prefix, args.chrm_sizes, \
        args.track_name_prefix, filter_region

def main():
    wiggles, out_fname_prefix, chrm_sizes_fp, track_name_prefix, filter_region \
        = parse_arguments()
    
    merged_wiggle = Wiggle( chrm_sizes_fp )

    for input_wiggle_fp in wiggles:
        if VERBOSE:
            print "Adding ", input_wiggle_fp.fname
        merged_wiggle.load_data_from_fp( input_wiggle_fp )

    if VERBOSE:
        print 'Writing merged wiggle...'
    
    ofn_template = out_fname_prefix + ".{strand}.bedGraph"
    merged_wiggle.write_wiggles( ofn_template.format(strand='plus'),  \
                                 ofn_template.format(strand='minus'), \
                                 ignore_zeros=True, \
                                 track_name_prefix=track_name_prefix, \
                                 filter_region = filter_region )
    
if __name__ == '__main__':
    main()
