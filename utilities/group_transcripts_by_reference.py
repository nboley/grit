import sys, os
from itertools import izip, combinations
from collections import defaultdict

sys.path.append( os.path.join( os.path.dirname( __file__), 
                               "../file_types/fast_gtf_parser/" ) )
from gtf import load_gtfs

# bijective hexavijezimal

def get_introns_by_gene_name( genes ):
    intron_to_gene_name = {}
    trans_to_ref_name = {}
    for ( gene_name, chr, strand, start, stop, transcripts ) in genes:
        for trans in transcripts:
            trans_to_ref_name[ tuple(trans[1][1:-1]) ] = trans[0]
            for intron in izip( trans[1][1:-2:2], trans[1][2:-1:2] ):
                intron_to_gene_name[ intron ] = gene_name
    
    return intron_to_gene_name, trans_to_ref_name

def find_associated_genes( trans, intron_to_gene_name ):
    genes = set()
    for intron in izip( trans[1][1:-2:2], trans[1][2:-1:2] ):
        if intron in intron_to_gene_name:
            genes.add( intron_to_gene_name[ intron ] )
    return tuple(sorted(genes))

def get_good_gene_merges( fname ):
    good_gene_merges = set()
    with open( fname ) as fp:
        for line in fp:
            gene_names = tuple(sorted(line.split()))
            good_gene_merges.add( gene_names )
            # add all of the combinations
            for i in xrange( 2, len(gene_names) ):
                for new_gene_names in combinations(gene_names, i):
                    good_gene_merges.add( tuple(sorted(new_gene_names)) )
    
    return good_gene_merges

class HexavijezimalCntr( object ):
    def __init__( self ):
        self.value = 0
    def __str__( self ):
        return str( self.value )
    def increment( self ):
        self.value += 1

def main( ann_fname, ref_fname, good_gene_merges_fname ):
    NOVEL_CNTR = 1
    # even if the transcript match is identical, we may have alternate 5' and 3'
    # ends which means multiple 'identical' transcripts. So keep a counter. If 
    # the counter is at 0 *for that ref trans id*, keep the ref name. Otherwise,
    # increment it by one. 
    IDENT_CNTR = defaultdict( HexavijezimalCntr )
    # similar for gene merges, but index by the tpule of gene names
    GENE_MERGE_CNTR = defaultdict( HexavijezimalCntr )
    # for novel splices, againb, index by gene name
    ALTERNATE_SPLICE_CNTR = defaultdict( HexavijezimalCntr )
    
    good_gene_merges = get_good_gene_merges( good_gene_merges_fname )
    novel_genes, ref_genes =  load_gtfs( ( ann_fname, ref_fname ), 2 )
    intron_to_gene_name, trans_to_ref_name = get_introns_by_gene_name(ref_genes)

    trans_id_to_names = {}
    for ( gene_name, chr, strand, start, stop, transcripts ) in novel_genes:
        for trans in transcripts:
            assoc_genes = find_associated_genes(trans, intron_to_gene_name)
            # if this is a completely novel transcript
            if len( assoc_genes ) == 0:
                trans_id_to_names[trans[0]] = ( 
                    gene_name, "NB%s" % str(NOVEL_CNTR).zfill( 8 ) )
                NOVEL_CNTR += 1
            # if there is only 1 corresponding ref gene
            elif len( assoc_genes ) == 1:
                # check to see if this is an exact match
                if tuple(trans[1][1:-1]) in trans_to_ref_name:
                    i_trans_name = trans_to_ref_name[tuple(trans[1][1:-1])]
                    if IDENT_CNTR[ i_trans_name ].value > 0:
                        new_i_trans_name = i_trans_name + ".%s" \
                            % (IDENT_CNTR[ i_trans_name ])
                    else:
                        new_i_trans_name = i_trans_name
                    
                    IDENT_CNTR[ i_trans_name ].increment()
                    
                    trans_id_to_names[ trans[0] ] = ( 
                        assoc_genes[0], new_i_trans_name )
                # if it's not, it must be a partial match
                else:
                    ALTERNATE_SPLICE_CNTR[ assoc_genes[0] ].increment()
                    i_trans_name = "%s.%s" % ( 
                        assoc_genes[0], ALTERNATE_SPLICE_CNTR[ assoc_genes[0] ])
                    trans_id_to_names[ trans[0] ] = ( 
                        assoc_genes[0], i_trans_name )
            # if this is a gene merge
            elif len( assoc_genes ) > 1:
                # make sure that it's white listed
                if assoc_genes not in good_gene_merges:
                    continue
                i_gene_name = "/".join( assoc_genes )
                ALTERNATE_SPLICE_CNTR[ i_gene_name ].increment()
                i_trans_name = "%s.%s" % ( 
                    i_gene_name, ALTERNATE_SPLICE_CNTR[ i_gene_name ] )
            
            print trans_id_to_names[ trans[0] ]
        raw_input()

if __name__ == "__main__":
    main( sys.argv[1], sys.argv[2], sys.argv[3] )
