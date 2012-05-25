import sys
import os
import numpy
import re

from collections import defaultdict
from pysam import Fastafile

VERBOSE = False
MIN_NUM_OVERLAP_BASES_DEFAULT = 3
MAX_INTRON_SIZE = 50000

sys.path.append( os.path.join( os.path.dirname(__file__), "../../", "file_types" ) )
from junctions_file import parse_junctions_file, build_jn_line, GenomicInterval, get_jn_type

def merge_jns( gff_fnames ):
    all_jns = defaultdict( lambda: defaultdict( lambda: defaultdict(int) ) )
    for sample_id, gff_fname in enumerate(gff_fnames):
        with open( gff_fname ) as gff_fp:
            jns = parse_junctions_file(gff_fp)
            if VERBOSE:
                print >> sys.stderr, "Finished parsing", gff_fp.name
            for jn, cnt, grp in jns.iter_jns_and_cnts_and_grps():
                all_jns[jn][sample_id][grp] += cnt
    
    return all_jns

def iter_filtered_jns( junctions, is_valid=lambda x: True, fasta_obj=None ):    
    jns_types = None
    if fasta_obj != None:
        jns_types = {}
        for jn, grp_data in junctions.iteritems():
            jns_types[jn] = get_jn_type( region.chr, region.start, region.stop, \
                                           fasta_obj, region.strand )
    
    def test_non_canonical( jn ):
        if jns_types == None:
            return True
        
        if jns_types[ jn ] == 'canonical_wrong_strand':
            new_strand = '-' if jn.strand == '+' else '+'
            new_jn = GenomicInterval( jn.chr, new_strand, jn.start, jn.stop )
            if new_jn in jns_types:
                return False
        
        return True

    for jn, grp_data in junctions.iteritems():
        if is_valid( jn, grp_data ) and test_non_canonical( jn ):
            cnt = sum( sum(entry.values()) for entry in grp_data.values() )
            yield jn, cnt, jns_types[jn] if jns_types != None else None
    
    return

def is_valid_jn_factory( min_entropy_score, \
                         min_num_success_samples, \
                         maximum_intron_size ):
    def calc_entopy_score( cnts ):
        cnts = numpy.array( cnts, dtype=float )
        if cnts.sum() < 1e-6: return 0
        
        pi_s = cnts/cnts.sum()
        return -(pi_s*numpy.log2(pi_s)).sum()
    
    
    def f( intron, grp_data ):
        num_pass_thresh = 0
        for key, val in grp_data.iteritems():
            if calc_entopy_score( val.values()) >= min_entropy_score:
                num_pass_thresh += 1
        
        if num_pass_thresh < min_num_success_samples:
            return False
        
        intron_size = intron[3] - intron[2] + 1
        if intron_size > maximum_intron_size:
            return False
        
        return True
    return f


def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(description=\
       'Parse junctions from provided bam file and write to std out.')
    parser.add_argument( 'jns_fns', nargs="+", help="GFF files to merge jns from." )
    parser.add_argument( '--fasta', '-s', default=None,\
       help='Fasta file used to determine the consensus intron type (fasta file ' + \
           'expects to be indexed: use samtools faidx).')    
    parser.add_argument( '--filter', '-f', default=False, action='store_true', \
       help='Whether or not to print status information.' )
    parser.add_argument( '--min-num-overlap-bases', type=int, \
                             default=MIN_NUM_OVERLAP_BASES_DEFAULT, \
       help='The minimum number of different start bases we must observe for a valid junction.' )
    parser.add_argument( '--maximum-intron-size', type=int, \
                             default=MAX_INTRON_SIZE, \
       help='The maximum allowable intron size.' )
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', \
       help='Whether or not to print status information.' )
    
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose

    return args.jns_fns, args.fasta, args.filter, \
                             args.min_num_overlap_bases, args.maximum_intron_size


def main():
    jns_fns, fasta_fn, do_filter, min_num_overlap_bases, maximum_intron_size \
        = parse_arguments()
    fasta_obj = None if fasta_fn == None else Fastafile( fasta_fn )
    
    jns = merge_jns( jns_fns )
    
    if do_filter==True:
        is_valid_jn = is_valid_jn_factory( 2.0, 2, maximum_intron_size )
    else:
        is_valid_jn = ( lambda intron, grp_data: True )
    
    for grp_id, (region, cnt, junction_type) in enumerate(sorted( \
            iter_filtered_jns( jns, is_valid_jn ) )):
        sys.stdout.write( build_jn_line( \
                region, grp_id, cnt, fasta_obj = fasta_obj, intron_type=junction_type) + "\n")
    
    return

if __name__ == "__main__":
    main()

