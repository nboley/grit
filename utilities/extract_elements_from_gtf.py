import os, sys

sys.path.append( os.path.join( os.path.dirname( __file__), 
                               "../file_types/fast_gtf_parser/" ) )
from gtf import load_gtf

sys.path.append( os.path.join( os.path.dirname( __file__), "../file_types/" ) )
from gtf_file import create_gff_line, create_gtf_line, GenomicInterval

def build_transcript_bnds_line( tr, gene_id ):
    start = tr.exons[0][0]
    stop = tr.exons[-1][1]
    region = GenomicInterval( tr.chrm, tr.strand, start, stop )
    return create_gtf_line( region, gene_id, tr.id, 
                            {'type': 'transcript'}, score='.', feature='exon'  )

def build_t_s_lines( gene, tss=True ):
    t_ss = set()
    for tr in gene[-1]:
        if ( tr.strand == '+' and tss ) or ( tr.strand == '-' and not tss ):
            coord = tr.exons[0][0]
        else:
            coord = tr.exons[-1][1]
        t_ss.add( coord )

    if tss:
        label = 'TSS'
    else:
        label = 'TES'
    
    lines = []
    for t_s in sorted(t_ss):
        region = GenomicInterval( tr.chrm, tr.strand, t_s, t_s )
        lines.append( create_gff_line( 
                region, 'gene_id "%s"; type "%s";' % (gene[0], label),
                score='.', feature='exon'  ) )
    
    return "\n".join( lines )

def build_gene_bnds_line( gene ):
    # find the gene bounds
    gene_start = 1e50
    gene_stop = 0
    for tr in gene[-1]:
        gene_start=min( gene_start, tr.exons[0][0] )
        gene_stop=max( gene_start, tr.exons[-1][1] )
    
    region = GenomicInterval( gene[1], gene[2], gene_start, gene_stop )
    return create_gff_line( region, gene[0], score='.', feature='exon' )
                            

def main( element_type, gtf_fname ):
    extract_gene_bnds = False
    extract_trans_bnds = False
    extract_tss_bnds = False
    extract_tes_bnds = False

    element_types = ("gene", "transcript", "TSS", "TES")
    if element_type not in element_types:
        raise ValueError, "Uncrecognized element type '%s' ( expecting %s) " % (
            element_type, ", ".join( element_types ) )
    
    if element_type == 'gene':
        extract_gene_bnds = True
    elif element_type == 'transcript':
        extract_trans_bnds = True
    elif element_type == 'TSS':
        extract_tss_bnds = True
    elif element_type == 'TES':
        extract_tes_bnds = True
    
    genes = load_gtf( gtf_fname )
    for gene in genes:
        if extract_gene_bnds:
            print build_gene_bnds_line( gene )
        if extract_tss_bnds:
            print build_t_s_lines( gene, tss=True )
        if extract_tes_bnds:
            print build_t_s_lines( gene, tss=False )
        
        for trans in gene[-1]:
            if extract_trans_bnds:
                print build_transcript_bnds_line( trans, gene[0] )
    
    return

if __name__ == '__main__':
    main( sys.argv[1], sys.argv[2] )
