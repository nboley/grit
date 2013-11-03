import sys

import sys, os
sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), ".." ) )

from grit.files.junctions import load_junctions_in_bam
from grit.files.bed import create_bed_line
from grit.files.reads import RNAseqReads, MergedReads, fix_chrm_name_for_ucsc
from grit.lib.multiprocessing_utils import ProcessSafeOPStream

def log_statment(statement):
    print >> sys.stderr, statement

VERBOSE = False
NTHREADS = 1
FIX_CHRM_NAMES_FOR_UCSC = False

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(
        description = 'Extract all junctions in a bam file' )
    parser.add_argument(
        'bam', type=argparse.FileType('rb'), nargs='+',
        help='Indexed BAM file(s) containing RNASeq to extract jns from' )
    parser.add_argument( '--rnaseq-read-type', nargs='+',
        choices=["forward", "backward"],
        help='Whether or not the first RNAseq read in a pair needs to be reversed to be on the correct strand.')
    parser.add_argument(
        '--ofname', '-o', type=argparse.FileType('w'), default=sys.stdout,
        help='Output stream (default: stdout)')
    parser.add_argument( '--region', help='Region to find junctions in (ie: chr4:+:100000-200000)')
    parser.add_argument( '--ucsc', default=False, action='store_true',
        help='Try to format contig names in the ucsc format (typically by prepending a chr).')    
    parser.add_argument(
        '--nthreads', '-t', type=int, default=1,
        help='Number of threads (default: %(default)d)')
    parser.add_argument( 
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    if len(args.bam) != len(args.rnaseq_read_type):
        raise ValueError, "--rnaseq-read-type must be the same length as the bam argument"
    
    global VERBOSE
    VERBOSE = args.verbose
    global NTHREADS
    NTHREADS = args.nthreads
    global FIX_CHRM_NAMES_FOR_UCSC
    if args.ucsc: FIX_CHRM_NAMES_FOR_UCSC = True
    
    reverse_rnaseq_strand = ( 
        True if args.rnaseq_read_type == 'backward' else False )
    reads =  MergedReads([
        RNAseqReads(bam.name).init(reverse_read_strand=reverse_rnaseq_strand)
        for bam in args.bam ])
    
    return reads, args.region, args.ofname

def main():
    reads, region, ofp = parse_arguments()
    if region != None:
        chrm, strand, pos = region.split(":")
        start, stop = [int(x) for x in pos.split("-")]
        region = [(chrm, strand, start, stop),]
    jns = load_junctions_in_bam(reads, regions=region, nthreads=NTHREADS)
    print >> ofp, "track name=junctions useScore=1"
    for (contig, strand), contig_jns in jns.iteritems():
        if FIX_CHRM_NAMES_FOR_UCSC: contig = fix_chrm_name_for_ucsc(contig)
        for (start, stop), cnt in contig_jns:
            print >> ofp, create_bed_line(
                contig, strand, start, stop, name='intron', score=cnt)
                            

if __name__ == '__main__':
    main()
