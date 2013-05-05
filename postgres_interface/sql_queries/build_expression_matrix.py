import os, sys

import numpy
import psycopg2

from itertools import izip

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='Estimate transcript frequencies from DB.')

    parser.add_argument( '--host', default='localhost', help='Database host.' )
    parser.add_argument( '--db-name', default='rnaseq_data', 
                         help='Database name. ' )

    args = parser.parse_args()
        
    conn = psycopg2.connect("dbname=%s host=%s user=nboley" % (args.db_name, args.host))
    
    return conn

def get_rpkm_for_fn( cursor, fn ):
    query = """
SELECT rpkm FROM (
SELECT transcript, COALESCE( rpkm, 0 ) as rpkm
  FROM (   SELECT transcript, rpkm 
            FROM transcript_expression
           WHERE reads_fn = '%s' ) as foo
    RIGHT JOIN annotations.transcripts 
    ON ( foo.transcript = annotations.transcripts.id  ) 
ORDER BY transcript
) foo2;
""" % fn
    cursor.execute( query )
    return numpy.array( cursor.fetchall() )[:,0]

def main():
    conn = parse_arguments()
    cursor = conn.cursor()

    cursor.execute( "select distinct reads_fn from transcript_expression;" )
    fnames = sorted( x[0] for x in cursor.fetchall() )
    
    cursor.execute( "select id from annotations.transcripts order by id;" )
    transcript_names = [ x[0] for x in cursor.fetchall() ]
    
    data = []
    for fname in fnames:
        print >> sys.stderr, fname
        data.append( get_rpkm_for_fn( cursor, fname ) )
    
    print "transcript\t" + "\t".join( "_".join(fname.split("/")[3:5]) for fname in fnames )
    for transcript, values in izip( transcript_names, numpy.vstack( data ).T):
        print transcript + "\t" + "\t".join( "%e" % x for x in values )
    
if __name__ == '__main__':
    main()
