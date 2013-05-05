import os, sys
import psycopg2

sys.path.append( os.path.join(os.path.dirname(__file__), "../file_types/") )
from gtf_file import Gene, Transcript

VERBOSE = False

def parse_integer_range( interval ):
    start, stop = interval[1:-1].split(",")
    return int(start), int(stop)-bool(interval[-1]==')')

def parse_contig_composite( data ):
    name, ann_id = data[1:-1].split(",")
    return ( name, int(ann_id) )

def load_gene_object( gene_id, conn ):
    cursor = conn.cursor()
    query = """SELECT contig, strand, location 
               FROM annotations.genes WHERE id = '%s';""" % gene_id
    cursor.execute( query )
    res = cursor.fetchall()
    if len(res) == 0: 
        raise ValueError, "Unknown gene: %s" % gene_id
    assert len(res) == 1
    
    chrm, strand, region = res[0]
    start, stop = parse_integer_range(region)
    cursor.close()
    
    return Gene( gene_id, parse_contig_composite(chrm), 
                 strand, start, stop, [] )

def load_transcript_ids( gene, conn ):
    cursor = conn.cursor()
    query = "SELECT id FROM annotations.transcripts WHERE gene = %s;" % gene.id
    cursor.execute( query )
    rv = [ x[0] for x in cursor.fetchall() ]
    cursor.close()
    return rv

def load_transcript( gene, trans_id, conn ):
    cursor = conn.cursor()
    # load the exons
    query = "SELECT location FROM annotations.exons " \
          + "WHERE transcript = '%s';" % trans_id
    cursor.execute( query )
    exons = [ parse_integer_range(x[0]) for x in cursor.fetchall() ]
    
    trans = Transcript( trans_id, gene.chrm, gene.strand, exons, None, 
                        gene_id=gene.id )
    cursor.close()
    return trans
    
    # try and load the coding sequence
    query = "SELECT region FROM annotations.transcript_regions " \
          + "WHERE transcript = '%s' and type='CDS';" % trans_id
    cursor.execute( query )
    res = cursor.fetchone()
    if len( res ) > 0:
        cds_region = parse_integer_range( res[0] )
        genome_cds_region = ( trans.genome_pos(cds_region[0]), 
                              trans.genome_pos(cds_region[1]) ) 
        trans.add_cds_region( genome_cds_region )
    
    return trans

def load_transcripts( gene, conn ):
    assert len( gene.transcripts ) == 0
    transcript_ids = load_transcript_ids( gene, conn )
    for t_id in transcript_ids:
        gene.transcripts.append( load_transcript( gene, t_id, conn ) )
    
    return

def load_gene_from_db( gene_id, conn):
    gene = load_gene_object( gene_id, conn )
    load_transcripts( gene, conn )
    return gene

#print load_gene_from_db( 'RAB4B', 1 )
