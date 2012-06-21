# Copyright (c) 2011-2012 Nathan Boley

import os
import sys
import numpy
import subprocess

sys.path.append( os.path.join( os.path.dirname(__file__), "../sparsify/" ) )
from sparsify_transcripts import build_gene_transcripts, build_reads_objs

from gene_models import parse_gff_line, GeneBoundaries
from reads import BinnedReads
from collections import defaultdict 

MINIMAL_VERBOSE = True



def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Determine valid transcripts and estimate frequencies.')
    parser.add_argument( 'gtf', type=file, \
                             help='GTF file processed for expression')
    parser.add_argument( 'bam_fns', nargs='+', metavar='bam',\
                             help='list of bam files to for which to produce expression')
    parser.add_argument( '--fl_dist', \
                             help='a pickled fl_dist object(default:generate fl_dist from input bam)')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', \
                             help='Whether or not to print status information.')
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    # we change to the output directory later, and these files need to opened in 
    # each sub-process for thread safety, so we get the absokute path while we can.
    bam_fns = [ os.path.abspath( bam_fn ) for bam_fn in args.bam_fns ]
    
    return args.gtf, bam_fns, args.fl_dist

def build_binned_reads_mat( ):
    pass


def main():
    # Get file objects from command line
    gtf_fp, bam_fns, fl_dist_fn = parse_arguments()
    
    # create gene object and close gtf_file
    genes = GeneBoundaries( gtf_fp )
    if MINIMAL_VERBOSE:
        print "Built gene objects from gtf file."
    # move file position back to beginning of file to be read for creating transcripts
    gtf_fp.seek(0)
    
    # create gene_transcripts dict
    gene_transcripts = build_gene_transcripts( gtf_fp, genes )
    gtf_fp.close()
    if MINIMAL_VERBOSE:
        print "Built transcripts objects from gtf file."
    
    # load the read objects, possible loading an already estimated fl dist
    reads_objs = build_reads_objs( bam_fns, fl_dist_fn )
    num_reads_dict = {}
    for bam_fn in sorted(bam_fns):
        res = subprocess.Popen(["samtools idxstats " + bam_fn,], shell=True, stdout=subprocess.PIPE )
        output = res.communicate()[0]
        # int(i.split()[2])
        num_reads_dict[ bam_fn ] = sum( int(i.split()[2])  for i in output.split("\n") if i != "" )
    
    for gene in genes.values():
        bin_types = set()
        
        # get dictionaries of the binned read cnts
        binned_read_cnt_dicts = []
        keys = sorted( sorted(reads_objs) )
        for reads_fn in keys:
            reads = reads_objs[ reads_fn ]
            br_cnts = BinnedReads( gene, reads, reads.read_group_mappings ).binned_reads
            bins = bin_types.update( item[2] for item in br_cnts.keys() )
            new_cnts = defaultdict( int )
            for key, val in br_cnts.iteritems():
                new_cnts[ key[2] ] += val
            binned_read_cnt_dicts.append( new_cnts )
        
        bin_types = sorted( bin_types )
        if len( bin_types ) == 0: continue
        
        # build the array
        cnts_array = []
        num_reads = [ num_reads_dict[key] for key in keys ]
        for bin_type in bin_types:
            cnts = [ cd[ bin_type ] for cd in binned_read_cnt_dicts ]
            # cnts_array.append( [ bin_type, ] )
            cnts_array.append( cnts)
        
        def calc_test_stat( cnts1, N1, cnts2, N2 ):
            return (abs((cnts1/N1) - (cnts2/N2))).sum()
        
        def sample_from_multinomials( cnts1, N1, cnts2, N2, N):
            p1 = cnts1/N1

            p2 = cnts2/N2
            
            p = (cnts1 + cnts2)/( N1 + N2 )
            
            NULL_stats = []
            for i in xrange( N ):
                s_1 = numpy.array( numpy.random.multinomial( N1, p, 1 ), dtype=float )
                s_2 = numpy.array( numpy.random.multinomial( N2, p, 1 ), dtype=float )
                NULL_stats.append( calc_test_stat( s_1, N1, s_2, N2 ) )
            
            #import pylab as P
            #n, bins, patches = P.hist(NULL_stats, 50, normed=1, histtype='stepfilled')
            #P.show()
            
            return numpy.array( NULL_stats )
        
        cnts_array = numpy.array( cnts_array, dtype=float )
        if cnts_array.sum() < 30: 
            continue
        
        
        test_stat = calc_test_stat( cnts_array[:,0], num_reads[0], 
                                    cnts_array[:,1], num_reads[1] )
        NULL_stats = sample_from_multinomials( cnts_array[:,0], num_reads[0], 
                                               cnts_array[:,1], num_reads[1], 1000 )
        NULL_mean = NULL_stats.mean()
        NULL_sd = NULL_stats.std()
        z_score = ( test_stat - NULL_mean )/NULL_sd
        print "%s\t%.2f\t%.2e\t%.2e\t%.2e" % ( gene.name.ljust(20), z_score, NULL_mean, test_stat, NULL_sd )
        
        # sys.exit()
            
if __name__ == "__main__":
    main()
