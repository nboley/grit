# Copyright (c) 2011-2012 Nathan Boley

import os, sys

sys.path.append( os.path.join( os.path.dirname( __file__), 
                               "../file_types/" ) )
from gtf_file import load_gtf

sys.path.append( os.path.join( os.path.dirname( __file__), "../file_types/" ) )
from bed import create_bed_line

                            
def main( gtf_fname ):
    genes = load_gtf( gtf_fname )
    
    single_exon_genes = set()
    tss_exons = set()
    internal_exons = set()
    introns = set()
    tes_exons = set()
    
    for gene in genes:
        for t in gene.transcripts:
            if len( t.exons ) == 1:
                single_exon_genes.add( 
                    ( t.chrm, t.strand, t.exons[0][0], t.exons[0][1] ) )
            else:
                tss_exon = t.exons[0] if t.strand == '+' else t.exons[-1]
                tss_exons.add( (t.chrm, t.strand, tss_exon[0], tss_exon[1]) )
                internal_exons.update( (t.chrm, t.strand, x1, x2) 
                                       for x1, x2 in t.exons[1:-1] )
                introns.update( (t.chrm, t.strand, x1, x2) for x1, x2 in t.introns )
                tes_exon = t.exons[-1] if t.strand == '+' else t.exons[0] 
                tes_exons.add( (t.chrm, t.strand, tes_exon[0], tes_exon[1]) )
    
    print 'track name="extracted_elements" visibility=2 itemRgb="On"'
    data = [ single_exon_genes, tss_exons, internal_exons, introns, tes_exons ]
    labels = [ "single_exon_gene", "tss_exon", "internal_exon", 
               "intron", "tes_exon" ]
    colors = [ '000,000,000', '140,195,59', '000,000,000',
               '200,200,200', '255,51,255' ]
    for exons, label, color in zip( data, labels, colors ):
        for exon in exons:
            print create_bed_line( 
                exon[0], exon[1], exon[2], exon[3]+1, 
                name=label, color=color,
                use_thick_lines=bool( label != 'intron' ) )
    
    return

if __name__ == '__main__':
    main( sys.argv[1] )
