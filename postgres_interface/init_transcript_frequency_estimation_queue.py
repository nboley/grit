import os, sys
import pickle

VERBOSE = False
import psycopg2

from load_gene_from_db import load_gene_from_db

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='Estimate transcript frequencies from DB.')

    # TODO add code to determine this automatically
    parser.add_argument( '--annotation-name', required=True, type=int,
                         help='Annotation for which to estimate transcript frequencies.')

    parser.add_argument( '--bam-fname', required=True, type=file,
                         help='File containing the reads to estimate transcript frequencies from.')
    parser.add_argument( '--fl-dists', type=file,
                         help='File containing the estimated fragment length dist. ( default: BAM_FN.fl_dist )')
    
    parser.add_argument( '--host', default='localhost', help='Database host.' )
    parser.add_argument( '--db-name', default='rnaseq_data', 
                         help='Database name. ' )
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
                             help='Print status information.')

    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    # load the fl dists
    fl_dists_fp = open( args.bam_fname.name + ".fldist" ) \
        if args.fl_dists == None else args.fl_dists
    
    conn = psycopg2.connect("dbname=%s host=%s" % ( args.db_name, args.host) )
    return conn, args.annotation_name, args.bam_fname, fl_dists_fp

def main():
    conn, ann_id, bam_fp, fl_dists_fp = parse_arguments()
    # add the bam, bam.bai, and fl_dist object to S3
    # 
    # XXX for now, just add the actual filename
    
    # 
    # add the genes and file to the db table
    # 
    cursor = conn.cursor()
    query = "INSERT INTO gene_expression_queue " \
          + "    SELECT id as gene, " \
          + "           %i as annotation, " % ann_id \
          + "         '%s' as bam_fn " % bam_fp.name \
          + "    FROM genes WHERE annotation = %i;" % ann_id
    cursor.execute( query )
    conn.commit()

    return

if __name__ == '__main__':
    main()
