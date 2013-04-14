import os, sys
import numpy
import psycopg2
import pickle

"""
CREATE OR REPLACE VIEW single_exon_genes AS
  SELECT genes.id AS gene, genes.annotation as annotation
    FROM genes, transcripts, exons
   WHERE genes.id = transcripts.gene
     AND genes.annotation = transcripts.annotation
     AND transcripts.id = exons.transcript
GROUP BY genes.id, genes.annotation
  HAVING count(*) = 1;

CREATE VIEW single_exon_gene_locations AS
SELECT exons.contig, exons.strand, exons.location 
  FROM single_exon_genes, transcripts, exons
 WHERE single_exon_genes.gene = transcripts.gene
   AND single_exon_genes.annotation = transcripts.annotation
   AND transcripts.id = exons.transcript;
"""

from load_gene_from_db import load_gene_from_db

sys.path.append( os.path.join(os.path.dirname(__file__), "../sparsify/") )
from frag_len import FlDist

sys.path.append( os.path.join(os.path.dirname(__file__), "../file_types/") )
from reads import Reads

import pysam

def build_fl_dist(bam_fn, conn, nsamples=25000, min_exon_size=600):
    reads = pysam.Samfile( bam_fn, "rb" )
    reads = Reads( bam_fn, "rb" )
    
    cursor = conn.cursor()
    query = """SELECT * FROM annotations.single_exon_genes"""
    cursor.execute( query )
        
    frag_lens = []
    for gene_id in cursor.fetchall():
        gene = load_gene_from_db( gene_id, conn )
        exon = gene.transcripts[0].exons[0]
        chrm = 'chr%s' % gene.chrm[0]
        if exon[1] - exon[0] < min_exon_size: continue
        try:
            for r1, r2 in reads.iter_paired_reads(
                    chrm, gene.strand, exon[0], exon[1], min_read_len=10, 
                    ignore_partial_alignments=False):
                f_start = min( r1.pos, r2.pos )
                f_stop = max( r1.aend, r2.aend )
                frag_lens.append( f_stop - f_start + 1 )
        except ValueError, inst:
            print "ValueError:, ", inst
            continue
        
        if len( frag_lens ) > nsamples:
            break

    cursor.close()
    
    # trim the top and bottom 1%
    frag_lens.sort()
    n_frags = len( frag_lens )
    frag_lens = frag_lens[int(n_frags*0.01):-int(n_frags*0.01)]

    # aggregate the fragment lengths, and build an flDist object
    min_fl, max_fl = frag_lens[0], frag_lens[-1]
    frag_dist = numpy.zeros( max_fl - min_fl + 1 )
    for frag_len in frag_lens:
        frag_dist[ frag_len-min_fl ] += 1
    fl_dist = FlDist( int(min_fl), int(max_fl), frag_dist )
    
    with open( bam_fn + ".fldist", "w" ) as ofp:
        pickle.dump( { 'mean': fl_dist }, ofp )

def main():
    bam_fn = sys.argv[1]
    conn = psycopg2.connect("dbname=%s host=%s user=nboley" 
                            % ('rnaseq', 'localhost'))
    build_fl_dist(bam_fn, conn)
    
if __name__ == '__main__':
    main()
