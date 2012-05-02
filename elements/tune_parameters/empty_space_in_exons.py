MIN_AVG_CVRG = 2.0

import sys, os

import numpy
import pickle

import matplotlib.pyplot as plt

sys.path.append( os.path.join( os.path.dirname(__file__), "..", 'file_types' ) )
from wiggle import wiggle
from exons_file import parse_exons_file

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description='Get transcript start sites(TSS) from a variety of sources.')
    parser.add_argument( '--exons', '-e',type=file, \
        help='GTF file which contains exons.')
    parser.add_argument( '--chrm_sizes_fname', '-c', type=file, \
        help='File with chromosome names and sizes.')
    parser.add_argument( '--cvrg_wigs', '-w', type=file, nargs='*', \
        help='WIG file with RNAseq read coverage. Note: strand will be infered from filename (plus/minus).')
    
    parser.add_argument( '--min_avg_cvrg', '-m', type=int, default=MIN_AVG_CVRG, \
        help='Exons with average coverage below this value will not be considered for empty regions. ' + \
            'default: %(default)f')
    
    parser.add_argument( '--out_fname', '-o', type=argparse.FileType( 'w' ), \
        help='Output file to dump exon counts.')
    parser.add_argument( '--cached_empty_regions', '-p', type=file,\
        help='Previously cached out_fname from this script.')
    args = parser.parse_args()
    
    if args.cached_empty_regions == None:
        if not all( item != None for item in [ args.exons, args.chrm_sizes_fname, args.cvrg_wigs ] ):
            parser.print_help()
            sys.exit()
    
    return args.exons, args.cvrg_wigs, args.chrm_sizes_fname, args.out_fname, args.cached_empty_regions

def main():
    exons_fp, cvrg_fps, chrm_sizes_fp, out_fp, pickled_fp = parse_arguments()
    
    if pickled_fp == None:
        print 'loading exons...'
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
        
        print 'calculating zero intervals...'
        cvrg.calc_zero_intervals()
        
        print 'getting largest empty regions in each exon...'
        exon_empty_regions = []
        keys = set( all_exons ).intersection( cvrg )
        for key in keys:
            for start, stop in all_exons[key]:
                if cumsum_cvrg[key][stop] - cumsum_cvrg[key][start-1] / (stop - start + 1) < MIN_AVG_CVRG:
                    continue
                
                max_empty_interval = 0
                # search for the first empty interval where empty_region_stop >= exon_start
                # Note empty regions are non-overlapping
                index = cvrg.zero_intervals[key][:,1].searchsorted( start )
                while index < len( cvrg.zero_intervals[key] ) and \
                        cvrg.zero_intervals[key][index][0] <= stop:
                    max_empty_interval = \
                        max( cvrg.zero_intervals[key][index][1] - cvrg.zero_intervals[key][index][0] + 1, \
                                 max_empty_interval )
                    index += 1
                
                exon_empty_regions.append( max_empty_interval )
        
        if out_fp != None:
            print 'dumping empty region sizes...'
            pickle.dump( exon_empty_regions, out_fp )
            out_fp.close()
    
    else:
        exon_empty_regions = pickle.load( pickled_fp )
    
    print min( exon_empty_regions ), max( exon_empty_regions )
    print len( [i for i in exon_empty_regions if i == 0] )
    print len( exon_empty_regions )
    plt.hist( exon_empty_regions, bins=500, range=(1,500), cumulative=True )
    plt.xlabel( 'Empty Region Size (1-500 bps)' )
    plt.ylabel( 'Num Exons (cumulative)' )
    plt.title( 'Empty Region Size Overlapping Flybase Exons\nMinimum Exon Coverage: {0:f}'.format(MIN_AVG_CVRG) )
    plt.show()
    
    return

if __name__ == '__main__':
    main()

