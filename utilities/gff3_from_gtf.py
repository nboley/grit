import os, sys
from collections import defaultdict

sys.path.append( os.path.join( os.path.dirname( __file__ ),
                               "..", "file_types" ) )
from gtf_file import create_gff3_line
sys.path.append( os.path.join( os.path.dirname( __file__ ),
                               "..", "file_types", "fast_gtf_parser" ) )
from gtf import load_gtf

def main():
    genes = load_gtf( sys.argv[1] )
    for gene in genes:
        gene_lines = []
        # build the gene line
        gene_lines.append( create_gff3_line( 
                gene.chrm, gene.strand, gene.start, gene.stop,
                'gene', gene.id, gene.id ) )
        
        # build the transcript lines, and group the exons
        exons = defaultdict( list )
        for t in gene.transcripts:
            gene_lines.append( create_gff3_line( 
                t.chrm, t.strand, t.start, t.stop,
                'mRNA', t.id, t.id, parents=[gene.id,] ) )
            for exon in t.exons:
                exons[exon].append( t.id )
        
        # build the exon lines
        for i, ((start, stop), parents) in enumerate(exons.iteritems()):
            exon_id = gene.id + "_E%i" % i
            gene_lines.append( create_gff3_line( 
                gene.chrm, gene.strand, start, stop,
                'exon', exon_id, exon_id, parents=parents ) )
            
        print "\n".join( gene_lines )
    
    return

main()
