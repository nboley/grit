# Copyright (c) 2011-2012 Nathan Boley

import sys
import os
import numpy
import re

from collections import defaultdict
from pysam import Fastafile

VERBOSE = False
MIN_CANONICAL_ENTROPY_SCORE = 1.5
MIN_NON_CANONICAL_ENTROPY_SCORE = 1.7
MIN_NUM_SAMPLES = 2
MIN_SCORE_RATIO = 0.01
MAX_INTRON_SIZE = 50000
THRESHHOLD_BY_UNIQ = False

sys.path.append(os.path.join(os.path.dirname(__file__), "../../", "file_types"))
import junctions_file
from junctions_file import parse_jn_gffs, Junction, GenomicInterval, get_jn_type


def merge_jns( gff_fnames, num_threads=1 ):
    all_jns = parse_jn_gffs( gff_fnames, min( num_threads, len( gff_fnames ) ) )
    
    grpd_jns = defaultdict( lambda: defaultdict( list ) )
    for sample_id, jns in enumerate(all_jns):
        for jn in jns:
            grpd_jns[ jn.region ][ sample_id ].append( jn )
        
    return grpd_jns

def calc_entopy_score( cnts ):
    cnts = numpy.array( cnts, dtype=float )
    if cnts.sum() < 1e-6: return 0
    
    pi_s = (cnts+1)/cnts.sum()
    return -(pi_s*numpy.log2(pi_s)).sum()


def iter_filtered_jns( grpd_jns, fasta_obj=None ):    
    # find the junction type for each intron
    jns_types = None
    if fasta_obj != None:
        jns_types = {}
        for jn_region in grpd_jns.iterkeys():
            jns_types[jn_region] = get_jn_type( 
                region.chr, region.start, region.stop, fasta_obj, region.strand)
    
    if VERBOSE: 
        print >> sys.stderr, "Finished building the junction types object."
    
    # find the junction count at each splice site
    splice_site_counts = defaultdict( int )
    for jn_region, src_grpd_jns in grpd_jns.iteritems():
        cnt_sum = 0
        uniq_cnt_sum = 0
        for jns in src_grpd_jns.itervalues():
            for jn in jns:
                cnt_sum += jn.cnt
                uniq_cnt_sum += jn.uniq_cnt
        
        splice_site_counts[ (jn.strand, jn.chrm, jn.start) ] = \
            max( splice_site_counts[ (jn.strand, jn.chrm, jn.start) ], cnt_sum )
        splice_site_counts[ (jn.strand, jn.chrm, jn.stop) ] = \
            max( splice_site_counts[ (jn.strand, jn.chrm, jn.stop) ], cnt_sum )
    
    if VERBOSE: 
        print >> sys.stderr, "Finished building the splice site counts objects."
    
    def test_entropy_score( jn_region, src_grpd_jns, jn_type ):
        if jn_type == 'canonical': 
            min_entropy_score = MIN_CANONICAL_ENTROPY_SCORE
        else:
            min_entropy_score = MIN_NON_CANONICAL_ENTROPY_SCORE
                
        num_pass_thresh = 0
        for jns in src_grpd_jns.itervalues():
            if THRESHHOLD_BY_UNIQ:
                scores = [ jn.uniq_cnt for jn in jns ]
            else:
                scores = [ jn.cnt for jn in jns ]
            
            if calc_entopy_score( scores ) >= min_entropy_score:
                num_pass_thresh += 1
        
        if num_pass_thresh < MIN_NUM_SAMPLES:
            return False
        
        intron_size = jn_region.stop - jn_region.start + 1
        if intron_size > MAX_INTRON_SIZE:
            return False
        
        return True
    
    def test_non_canonical( jn ):
        if jns_types == None:
            return True
        
        if jns_types[ jn ] == 'canonical_wrong_strand':
            new_strand = '-' if jn.strand == '+' else '+'
            new_jn_region = GenomicInterval( jn.chr, new_strand, jn.start, jn.stop )
            if new_jn_region in jns_types:
                return False
        
        return True

    def test_jn_ratio( jn, cnt ):
        max_val = max( splice_site_counts[ (jn.strand, jn.chr, jn.start) ], 
                       splice_site_counts[ (jn.strand, jn.chr, jn.stop) ]  )
        return float( cnt )/max_val > MIN_SCORE_RATIO
    

    for jn_region, src_grpd_jns in grpd_jns.iteritems():
        cnt = sum( sum( jn.cnt for jn in jns  ) 
                   for jns in src_grpd_jns.values())
        
        uniq_cnt = sum( sum( jn.uniq_cnt for jn in jns  ) 
                        for jns in src_grpd_jns.values())
        
        jn_type = jns_types[jn_region] if jns_types != None else None

        num_passed_thresh = 0
        for jns in src_grpd_jns.itervalues():
            if test_entropy_score( jn_region, src_grpd_jns, jn_type ) \
                    and test_non_canonical( jn_region ) \
                    and test_jn_ratio( jn_region, cnt ):
                num_passed_thresh += 1
            if num_passed_thresh >= MIN_NUM_SAMPLES: break
        
        if num_passed_thresh >= MIN_NUM_SAMPLES:
            yield Junction(jn_region,jn_type=jn_type,cnt=cnt,uniq_cnt=uniq_cnt)
    
    return

def iter_merged_jns( grpd_jns, fasta_obj=None ):
    for jn_region, src_grpd_jns in grpd_jns.iteritems():
        cnt = sum( sum( jn.cnt for jn in jns  ) 
                   for jns in src_grpd_jns.values())
        
        uniq_cnt = sum( sum( jn.uniq_cnt for jn in jns  ) 
                        for jns in src_grpd_jns.values())
        
        if fasta_obj != None:
            jn_type = get_jn_type( 
                jn_region.chr, jn_region.start, jn_region.stop, 
                fasta_obj, jn_region.strand)
        else:
            jn_type = None
        
        yield Junction(jn_region,jn_type=jn_type,cnt=cnt,uniq_cnt=uniq_cnt)
    
    return


def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(description=\
       'Parse junctions from provided bam file and write to std out.')
    parser.add_argument( 'jns_fns', nargs="+", \
                             help="GFF files to merge jns from." )
    parser.add_argument( '--fasta', '-s', default=None,\
       help='Fasta file used to determine the consensus intron type (fasta ' + \
           'file expects to be indexed: use samtools faidx).')    
    parser.add_argument( '--filter', '-f', default=False, action='store_true', \
       help='Whether or not to print status information.' )

    
    
    
    
    parser.add_argument( '--maximum-intron-size', type=int, \
                             default=MAX_INTRON_SIZE, \
       help='The maximum allowable intron size.' )

    parser.add_argument( '--threads', '-t', type=int, default=1,
       help='The number of threads to use.' ) 
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
       help='Whether or not to print status information.' )
    
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    junctions_file.VERBOSE = VERBOSE

    return args.jns_fns, args.fasta, args.filter, args.threads, args.maximum_intron_size


def main():
    jns_fns, fasta_fn, do_filter, num_threads, maximum_intron_size \
        = parse_arguments()
    fasta_obj = None if fasta_fn == None else Fastafile( fasta_fn )

    grpd_jns = merge_jns( jns_fns, num_threads )
    if VERBOSE: print >> sys.stderr, "Finished merging junctions."

    if do_filter:
        jns_iter = sorted( iter_filtered_jns( grpd_jns ) )
    else:
        jns_iter = sorted( iter_merged_jns( grpd_jns ) )
    
    if VERBOSE: print >> sys.stderr, "Finished filtering junctions."
        
    for grp_id, jn in enumerate( jns_iter ):
        sys.stdout.write( jn.build_gff_line(
                group_id = 'group_id "%i";' % (grp_id+1)) + "\n")
    
    return

if __name__ == "__main__":
    #import cProfile
    #cProfile.run( "main()" )
    main()
