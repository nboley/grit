import sys, os

import numpy
import pickle

import matplotlib.pyplot as plt

sys.path.append( os.path.join( os.path.dirname(__file__), "..", 'file_types' ) )
from wiggle import wiggle
from exons_file import parse_exons_file

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description='Plot the average coverage over all exons.')
    parser.add_argument( '--exons', '-e',type=file, \
        help='GTF file which contains exons.')
    parser.add_argument( '--chrm_sizes_fname', '-c', type=file, \
        help='File with chromosome names and sizes.')
    parser.add_argument( '--cvrg_wigs', '-w', type=file, nargs='*', \
        help='WIG file with RNAseq read coverage. Note: strand will be infered from filename (plus/minus).')
    
    parser.add_argument( '--out_fname', '-o', type=argparse.FileType( 'w' ), \
        help='Output file to dump exon counts.')
    parser.add_argument( '--cached_exon_counts', '-p', type=file,\
        help='Previously cached out_fname from this script.')
    args = parser.parse_args()
    
    if args.cached_exon_counts == None:
        if not all( item != None for item in [ args.exons, args.chrm_sizes_fname, args.cvrg_wigs ] ):
            parser.print_help()
            sys.exit()
    
    return args.exons, args.cvrg_wigs, args.chrm_sizes_fname, args.out_fname, args.cached_exon_counts

def main():
    exons_fp, cvrg_fps, chrm_sizes_fp, out_fp, pickled_fp = parse_arguments()
    
    if pickled_fp == None:
        all_exons = parse_exons_file( exons_fp )
        cvrg = wiggle( chrm_sizes_fp )
        chrm_sizes_fp.close()
        
        print 'building coverage...'
        for cvrg_fp in cvrg_fps:
            cvrg.load_data_from_fname( cvrg_fp )
            cvrg_fp.close()
        
        print 'converting coverage...'
        cumsum_cvrg = {}
        for key, chrm_cvrg in cvrg.iteritems():
            cumsum_cvrg[key] = numpy.cumsum( chrm_cvrg )
        
        print 'getting exon coveragees...'
        exon_counts = []
        keys = set( all_exons ).intersection( cumsum_cvrg )
        for key in sorted( keys ):
            for start, stop in all_exons[key]:
                exon_count = cumsum_cvrg[key][stop] - cumsum_cvrg[key][start-1]
                avg_cvrg = exon_count / (stop - start + 1)
                exon_counts.append( avg_cvrg )
        
        if out_fp != None:
            print 'dumping exons...'
            pickle.dump( exon_counts, out_fp )
            out_fp.close()
    
    else:
        exon_counts = pickle.load( pickled_fp )
    
    print min( exon_counts ), max( exon_counts )
    plt.hist( exon_counts, bins=101, range=(0,100) )
    plt.xlabel( 'Average Exon Coverage' )
    plt.ylabel( 'Num Exons' )
    plt.title( 'Exon Coverage (0-100 reads/bp)\n1 bp bins' )
    plt.show()
    
    return

if __name__ == '__main__':
    main()

