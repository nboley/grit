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

    # TODO add code to determine this automatically
    parser.add_argument( '--annotation-name', required=True, help='Database host. ' )
        
    parser.add_argument( '--host', default='localhost', help='Database host. ' )
    parser.add_argument( '--db-name', default='rnaseq_data', 
                         help='Database name. ' )
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
                             help='Whether or not to print status information.')
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    conn = psycopg2.connect("dbname=%s host=%s" % ( args.db_name, args.host) )

    return args.gtf, args.annotation_name, conn

def add_exons_to_db( transcript, conn ):
    cursor = conn.cursor()
    for (start, stop) in transcript.exons:
        query = "INSERT INTO exons " \
            + "( transcript, contig, strand, location )" \
            + "VALUES ( '%s', '(\"%s\", 1)', '%s', '[%s, %s]')" % \
            ( transcript.id, transcript.chrm, transcript.strand,
              start, stop )
        cursor.execute( query )
    cursor.close()
    

def add_transcript_to_db( gene, annotation_key, transcript, conn ):
    cursor = conn.cursor()
    query = "INSERT INTO transcripts " \
        + "( id, gene, annotation, contig, strand, location )" \
        + "VALUES ( '%s', '%s', %s, '(\"%s\", 1)', '%s', '[%s, %s]')" % \
        ( transcript.id, gene.id, annotation_key, transcript.chrm, 
          transcript.strand, transcript.start, transcript.stop )
    cursor.execute( query )
    cursor.close()
    
    add_exons_to_db( transcript, conn )
    
    return 

def add_gene_to_db( gene, conn, annotation_key ):
    cursor = conn.cursor()
    
    #add the gene entry
    query = "INSERT INTO genes ( id, annotation, contig, strand, location ) " \
        + "VALUES ( '%s', %s, '(\"%s\", 1)', '%s', '[%s, %s]')" % \
        ( gene.id, annotation_key, 
          gene.chrm, gene.strand, gene.start, gene.stop )
    cursor.execute( query )
    cursor.close()
    
    for trans in gene.transcripts:
        add_transcript_to_db( gene, annotation_key, trans, conn )

def add_annotation_to_db( conn, name, description='NULL' ):
    cursor = conn.cursor()
    cursor.execute("INSERT INTO annotations ( name, description )"
                   + "VALUES ( %s, %s ) RETURNING id", (name, description) )
    rv = cursor.fetchone()[0]
    cursor.close()
    return int(rv)

def main():
    gtf_fp, ann_name, conn = parse_arguments()
    
    # load the genes
    genes = load_gtf( gtf_fp.name )
    
    # add the annotation, and return the pkey
    annotation_key = add_annotation_to_db( conn, ann_name )
    
    # add all of the genes to the DB
    for i, gene in enumerate(genes):
        print "(%i/%i)    Processing %s" % ( i+1, len(genes), gene.id )
        add_gene_to_db( gene, conn, annotation_key )
    
    conn.commit()
    conn.close()

if __name__ == '__main__':
    main()
