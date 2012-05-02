#!/usr/bin/python

# import python mods
import os 
import sys

# declare constants
VERBOSE = False

BLAST_PRGM = 'blastp'
BLAST_DB = 'Microbial/7227'

E_MAX = 1e-10

# import biological mods
from Bio.Blast import NCBIWWW, NCBIXML

def run_blast( fasta_string, out_fp, saved_blast_fp, save_fp ):
    """ Either run blastp or load saved xml file
    """
    if saved_blast_fp != None:
        result_handle = saved_blast_fp
    else:
        print 'Running qblast remotely this may take some time.'
        result_handle = NCBIWWW.qblast( \
            BLAST_PRGM, BLAST_DB, fasta_string, expect=E_MAX )
        
        # if a save file is provided save the xml results
        if save_fp:
            save_fp.write( result_handle.read() )
            # seek back to read the file for parsing
            result_handle.seek(0)
            save_fp.close()
    
    # parse results for useful information. For available info see:
    #     http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc89
    blast_records = NCBIXML.parse( result_handle )
    
    # output description of hits for each query
    for blast_record in blast_records:
        out_fp.write( blast_record.query + '\n' )
        for description in blast_record.descriptions:
            dsc_line = '\t'.join( \
                (description.title, str(description.e), str(description.score)))
            out_fp.write( '\t' + dsc_line + '\n' )
        
        out_fp.write( '\n' )
    
    # result handle must be closed after the records are read
    result_handle.close()
    out_fp.close()
    
    return

def build_objects( fasta_fp ):
    # slurp whole file into a string
    # fix this to return a list of fasta record strings
    fasta_string = fasta_fp.read()
    
    return fasta_string

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(description='BlastP the records in the input fasta protein file.')
    parser.add_argument( '--fasta', '-f', type=file, \
        help='Fasta file of proteins to blast.')
    parser.add_argument( '--saved_blast', '-p', type=file, \
        help='XML from previously run blastp.')
    parser.add_argument( '--save_file', '-s',\
        help='Filename in which to save the XML blastp results. ' + \
                             '(default: do not save XML)' )
    global E_MAX
    parser.add_argument( '--e_val', '-e', type=float,\
        help='Maximum blast e-value hits to accept. default: {0:e}'.format( \
            E_MAX ))
    
    parser.add_argument( '--out_fname', '-o', \
        help='Output filename (default: stdout)')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', \
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    if not any( (args.fasta, args.saved_blast) ):
        print parser.print_help()
        sys.exit()
    
    # set global flag args
    global VERBOSE
    VERBOSE = args.verbose
    
    if args.e_val: E_MAX = args.e_val
    
    out_fp = open( args.out_fname, "w" ) if args.out_fname else sys.stdout
    save_fp = open( args.save_file, "w" ) \
        if args.save_file != None and args.saved_blast == None else None
    
    return args.fasta, out_fp, args.saved_blast, save_fp

def main():
    fasta_fp, out_fp, saved_blast_fp, save_fp = parse_arguments()
    fasta_string = build_objects( fasta_fp ) if saved_blast_fp == None else None
    
    run_blast( fasta_string, out_fp, saved_blast_fp, save_fp )
    
    return

if __name__ == '__main__':
    main()
