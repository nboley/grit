import os, sys

sys.path.append( os.path.join(os.path.dirname(__file__), 
                              "../file_types/fast_gtf_parser") )
from gtf import load_gtf
VERBOSE = False
import psycopg2

def add_exons_to_db( transcript, conn ):
    cursor = conn.cursor()
    for (start, stop) in transcript.exons:
        query = "INSERT INTO annotations.exons " \
            + "( transcript, location )" \
            + "VALUES ( '%s', '[%s, %s]')" % \
            ( transcript.id, start, stop )
        cursor.execute( query )
    cursor.close()
    return

def add_transcript_regions_to_db( transcript, conn ):
    cursor = conn.cursor()
    # add the coding sequence
    if transcript.is_protein_coding:
        start = transcript.relative_pos( transcript.cds_region[0] )
        stop = transcript.relative_pos( transcript.cds_region[1] )
        query = "INSERT into annotations.transcript_regions " \
              + "VALUES ( '%s', '[%s,%s]', '%s' )" \
              % (transcript.id, start, stop, 'CDS')
        cursor.execute( query )
    
    cursor.close()
    return

def add_transcript_to_db( gene_id, (chrm, assembly_id), transcript, conn ):
    assert chrm == transcript.chrm
    cursor = conn.cursor()
    query = "INSERT INTO annotations.transcripts " \
        + "( id, gene, contig, strand, location )" \
        + "VALUES ( '%s', %s, '(\"%s\", %s)', '%s', '[%s, %s]')" % \
        ( transcript.id, gene_id, chrm, assembly_id,
          transcript.strand, transcript.start, transcript.stop )
    cursor.execute( query )
    cursor.close()
    
    add_exons_to_db( transcript, conn )
    
    add_transcript_regions_to_db( transcript, conn )
    
    return 

def add_gene_to_db( gene, annotation_key, assembly_id, conn ):
    cursor = conn.cursor()
    
    #add the gene entry
    query_template = "INSERT INTO annotations.genes " \
        + "( name, annotation, contig, strand, location ) " \
        + "VALUES ( '%s', %s, '(\"%s\", %i)', '%s', '[%s, %s]') " \
        + "RETURNING id;"
    query = query_template % ( gene.id, annotation_key, gene.chrm, 
                               assembly_id, gene.strand, gene.start, gene.stop )
    cursor.execute( query )
    gene_id = int(cursor.fetchone()[0])
    cursor.close()
    
    for trans in gene.transcripts:
        add_transcript_to_db( gene_id, (gene.chrm, assembly_id), trans, conn )

def add_annotation_to_db( conn, name, description='NULL' ):
    cursor = conn.cursor()
    cursor.execute("INSERT INTO annotations.annotations ( name, description )"
                   + "VALUES ( %s, %s ) RETURNING id", (name, description) )
    rv = cursor.fetchone()[0]
    cursor.close()
    return int(rv)

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Load GTF file into database.')
    parser.add_argument( 'gtf', type=file, help='GTF file to load.')

    # TODO add code to determine this automatically
    parser.add_argument( '--annotation-name', required=True, 
           help='Name of the annotation we are inserting. ' )
    parser.add_argument( '--assembly-id', required=True, 
           help='ID of the assembly this annotation is based upon.' )

    parser.add_argument( '--host', default='localhost', help='Database host. ' )
    parser.add_argument( '--db-name', default='rnaseq_data', 
                         help='Database to insert the data into. ' )
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
                             help='Whether or not to print status information.')
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    conn = psycopg2.connect("dbname=%s host=%s" % ( args.db_name, args.host) )

    return args.gtf, args.annotation_name, int(args.assembly_id), conn

def main():
    gtf_fp, ann_name, assembly_id, conn = parse_arguments()
    
    # add the annotation, and return the pkey
    annotation_key = add_annotation_to_db( conn, ann_name )
    
    # load the genes
    genes = load_gtf( gtf_fp.name )
    
    # add all of the genes to the DB
    for i, gene in enumerate(genes):
        print "(%i/%i)    Processing %s" % ( i+1, len(genes), gene.id )
        add_gene_to_db( gene, annotation_key, assembly_id, conn )
    
    conn.commit()
    conn.close()

if __name__ == '__main__':
    main()
