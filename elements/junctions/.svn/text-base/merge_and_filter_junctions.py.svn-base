import sys
import os

from collections import defaultdict
from pysam import Fastafile

VERBOSE = False
MIN_NUM_OVERLAP_BASES_DEFAULT = 3
MAX_INTRON_SIZE = 50000

sys.path.append( os.path.join( os.path.dirname(__file__), "../../", "file_types" ) )
from junctions_file import parse_junctions_file, build_jn_line

def merge_jns( gff_fnames ):
    all_jns = defaultdict( lambda: defaultdict(int) )
    for gff_fname in gff_fnames:
        with open( gff_fname ) as gff_fp:
            jns = parse_junctions_file(gff_fp)
            if VERBOSE:
                print >> sys.stderr, "Finished parsing", gff_fp.name
            for jn, cnt, grp in jns.iter_jns_and_cnts_and_grps():
                all_jns[ jn ][ grp ] += cnt
    
    return all_jns

def iter_filtered_jns( junctions, is_valid=lambda x: True ):
    for jn, grp_data in junctions.iteritems():
        if is_valid( jn, grp_data ):
            yield jn, sum( grp_data.values() )
    return

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

def is_valid_jn_factory( min_num_overlap_bases, maximum_intron_size ):
    def f( intron, grp_data ):
        if len( grp_data ) < min_num_overlap_bases:
            return False
        
        intron_size = intron[3] - intron[2] + 1
        if intron_size > maximum_intron_size:
            return False
        
        return True
    return f

def main():
    jns_fns, fasta_fn, do_filter, min_num_overlap_bases, maximum_intron_size \
        = parse_arguments()
    fasta_obj = None if fasta_fn == None else Fastafile( fasta_fn )
    
    jns = merge_jns( jns_fns )
    
    if do_filter==True:
        is_valid_jn = is_valid_jn_factory( min_num_overlap_bases, maximum_intron_size )
    else:
        is_valid_jn = ( lambda intron, grp_data: True )
    
    for grp_id, (region, cnt) in enumerate(sorted( \
            iter_filtered_jns( jns, is_valid_jn ) )):
        sys.stdout.write( build_jn_line( \
                region, grp_id, cnt, fasta_obj = fasta_obj) + "\n")
    
    return

if __name__ == "__main__":
    main()

