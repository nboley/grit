import os, sys

sys.path.append( os.path.join(os.path.dirname(__file__), 
                              "../file_types/fast_gtf_parser") )
from gtf import load_gtf
VERBOSE = False
import psycopg2

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Load GTF file into database.')
    parser.add_argument( 'gtf', type=file, help='GTF file to load.')
        
    parser.add_argument( '--host', default='localhost', help='Database host. ' )
    parser.add_argument( '--db-name', default='rnaseq_data', 
                         help='Database name. ' )
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
                             help='Whether or not to print status information.')
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    conn = psycopg2.connect("dbname=%s host=%s" % ( args.db_name, args.host) )
    return args.gtf, conn

def add_gene_to_db( gene, conn ):
    for trans in gene.transcripts:
        print trans.chrm, trans.strand
        print trans.exons
        print trans.us_exons, trans.cds_exons, trans.ds_exons
        print trans
    sys.exit()

def main():
    gtf_fp, conn = parse_arguments()
    genes = load_gtf( gtf_fp.name )
    for gene in genes:
        add_gene_to_db( gene, conn )


if __name__ == '__main__':
    main()
