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
    parser.add_argument( '--bam-fname', required=True, type=str,
                         help='File containing the reads to estimate transcript frequencies from.')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
                             help='Print status information.')

    parser.add_argument( '--db-name', default='rnaseq', 
                         help='Database to insert the data into. ' )
    parser.add_argument( '--db-host', 
                         help='Database host. default: socket connection' )
    parser.add_argument( '--db-user', 
                         help='Database connection user. Default: unix user' )
    parser.add_argument( '--db-pass', help='DB connection password.' )

    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose

    conn_str = "dbname=%s" % args.db_name
    conn = psycopg2.connect(conn_str)
    
    return conn, args.annotation_name, args.bam_fname

def main():
    conn, ann_id, bam_fname = parse_arguments()
    # add the bam, bam.bai, and fl_dist object to S3
    # 
    # XXX for now, just add the actual filename
    if not bam_fname.startswith("s3" ):
        bam_fname = os.path.abspath( bam_fname )
    # 
    # add the genes and file to the db table
    # 
    cursor = conn.cursor()
    print "Inserting %s %s" % ( ann_id, bam_fname )
    query = "INSERT INTO gene_expression_queue " \
          + "    SELECT id as gene, " \
          + "         '%s' as bam_fn " % bam_fname \
          + "    FROM annotations.genes WHERE annotation = %i;" % ann_id
    cursor.execute( query )
    conn.commit()
    print "FINISHED Inserting %s %s" % ( ann_id, bam_fname )
    conn.close()
    
    return

if __name__ == '__main__':
    main()
