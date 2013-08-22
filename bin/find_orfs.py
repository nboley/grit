#!/usr/bin/python
import sys, os
sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), ".." ) )

import grit.proteomics.ORF
from grit.files.gtf import load_gtf
from grit.lib.multiprocessing_utils import ProcessSafeOPStream

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(
        description = 'Find open reading frames(ORF) in the input GTF file '
        'and output gtf file with annotated reading frames.' )
    parser.add_argument(
        'gtf', type=file,
        help='GTF file to search for ORFs.' )
    parser.add_argument(
        'fasta',
        help='Fasta file with reference sequence.' )
    parser.add_argument(
        '--min-aas', '-m', type=int,
        help='Number of amino acids to require for an open read frame. ' +
        '(default: {0:d})'.format( grit.proteomics.ORF.MIN_AAS_PER_ORF ))

    parser.add_argument(
        '--only-longest-orf', default=False, action='store_true',
        help='If this is set, only report the longest ORF per transcript. ' )
    parser.add_argument(
        '--dont-include-stop-codon', default=False, action='store_true',
        help='If this is set, don\'t include the stop codon in the CDS region.')
    
    parser.add_argument(
        '--threads', '-t', type=int, default=1,
        help='Number of threads with which to find ORFs. ' +
        '(default: %(default)d)')
    
    parser.add_argument(
        '--output-filename', '-o',
        help='Output file. (default: stdout)')
    parser.add_argument(
        '--fasta-output-filename', '-f',
        help='Write protein sequences to --fasta-output-filename, if set.')

    parser.add_argument( 
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    if args.min_aas != None: 
        grit.proteomics.ORF.MIN_AAS_PER_ORF = args.min_aas
    
    # create default if no prefix provided or if same as gtf filename
    if args.output_filename == None:
        gtf_ofp = ProcessSafeOPStream( sys.stdout )
    else:
        gtf_ofp = ProcessSafeOPStream( open( args.output_filename, 'w' ) )
        
    fa_ofp =  ProcessSafeOPStream( open(args.fasta_output_filename, 'w') ) if \
        args.fasta_output_filename != None else None
    
    # set flag args
    grit.proteomics.ORF.VERBOSE = args.verbose
    grit.proteomics.ORF.ONLY_USE_LONGEST_ORF = args.only_longest_orf
    grit.proteomics.ORF.INCLUDE_STOP_CODON = not args.dont_include_stop_codon
    
    return args.gtf, args.fasta, args.threads, gtf_ofp, fa_ofp

def main():
    gtf_fp, fasta_fn, threads, gtf_ofp, fa_ofp = parse_arguments()
    genes = load_gtf( gtf_fp.name )
    grit.proteomics.ORF.find_all_orfs(genes, fasta_fn, gtf_ofp, fa_ofp, threads)
    
    return

if __name__ == '__main__':
    main()
