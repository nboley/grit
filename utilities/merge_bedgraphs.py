#!/usr/bin/python

# Copyright (c) 2011-2012 Nathan Boley

# import python mods
import os 
import sys
import tempfile
import subprocess

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


def merge_bedgraphs( fps, final_output_fp ):
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
        final_output_fp.write( "\t".join( data ) + "\n" )
    
    union_bed_fp.close()
    
    return

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='Merge counts of many bedgraph files into one bedgraph file.')
    parser.add_argument( 'files', type=file, nargs='+', \
                             help='BedGraph files to merge into one wiggle.')
    
    parser.add_argument('--out-fname',help='Output filenames. Default: STDOUT.')
                     
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', \
                             help='Whether or not to print status information.')
    args = parser.parse_args()
    
    # set flag args
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.files, args.out_fname

def main():
    input_files, out_fname = parse_arguments()

    if out_fname == None:
        ofp = sys.stdout
    else:
        ofp = open( out_fname, "w" )
    
    if len( input_files ) == 1:
        import shutil
        shutil.copyfile( input_files[0].name, out_fname )
        print "Warning: only passed 1 input file. Making a copy.", input_files
    else:
        merge_bedgraphs( input_files, ofp )
    
    if out_fname != None:
        ofp.close()
    
    return
    
if __name__ == '__main__':
    main()
