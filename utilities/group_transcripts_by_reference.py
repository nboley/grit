import sys, os
from itertools import izip, combinations, chain
from collections import defaultdict
import re

sys.path.append( os.path.join( os.path.dirname( __file__), 
                               "../file_types/fast_gtf_parser/" ) )
from gtf import load_gtfs, Transcript

sys.path.append( os.path.join( os.path.dirname( __file__), "../file_types/" ) )
from gtf_file import create_gff_line, parse_gff_line, GenomicInterval

# bijective hexavijezimal

def get_introns_by_gene_name( genes ):
    intron_to_gene_name = defaultdict( lambda: defaultdict(set) )
    cds_bnds_to_gene_name = defaultdict( lambda: defaultdict(set) )
    trans_to_ref_name = defaultdict( lambda: defaultdict(set) )
    
    for ( gene_name, chrm, strand, start, stop, transcripts ) in genes:
        for trans in transcripts:
            trans_to_ref_name[(chrm, strand)][trans.IB_key()].add(
                (trans.id, gene_name ))
            if trans.is_protein_coding:
                cds_bnds_to_gene_name[(chrm, strand)][trans.cds_region].add(
                    gene_name)
            for bndry in trans[1][1:-1]:
                intron_to_gene_name[(chrm, strand)][bndry].add(gene_name)
    
    return intron_to_gene_name, trans_to_ref_name, cds_bnds_to_gene_name

def find_associated_genes( trans, intron_to_gene_name ):
    # genes which are associated with every bndry
    fully_assoc_genes = set()
    # genes associated with *any* bndry
    genes = set()
    
    # initialize the gene sets with genes from the first boundary
    first_bndry = trans[1][1]
    if first_bndry in intron_to_gene_name:
        genes.update( intron_to_gene_name[ first_bndry ] )
        fully_assoc_genes.update( intron_to_gene_name[ first_bndry ] )
    
    # for a gene to be fully associated it has to be associated with *every* 
    # boundary, so we take the intersection. to be associated, it has to be 
    # associated with any, so we take the union
    for bndry in trans[1][2:-1]:
        if bndry in intron_to_gene_name:
            assoc_genes = intron_to_gene_name[ bndry ] 
            genes.update( assoc_genes )
            if len( assoc_genes ) > 0:
                fully_assoc_genes.intersection_update(
                    intron_to_gene_name[ bndry ])
    
    return tuple(sorted(fully_assoc_genes)), tuple(sorted(genes))

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

"""
def get_blacklist_regions( fname, ref_genes ):
    gene_regions = {}
    for gene_name, chrm, strand, start, stop, transcripts in ref_genes:
        gene_regions[gene_name] = ( chrm, strand, start, stop )
    
    blacklist_regions = defaultdict( set )
    with open( fname ) as fp:
        for line in fp:
            gene_name = line.strip()
            try:
                chrm, strand, start, stop = gene_regions[ gene_name ]
            except KeyError:
                print >> sys.stderr, \
                    "WARNING: Can't find black list gene '%s'" % gene_name
                continue
            
            blacklist_regions[(chrm, strand)].add( (start, stop ) )
    
    for key in blacklist_regions:
        blacklist_regions[key] = sorted( blacklist_regions[key] )
    
    return blacklist_regions
"""

def get_blacklist_genes( fname ):
    with open( fname ) as fp:
        return set( fp.read().strip().split() )

class HexavijezimalCntr( object ):
    def __init__( self ):
        self.value = 0
    
    _char_mapping = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.lower()
    
    def __str__( self ):
        return_chars = []
        assert self.value > 0
        tmp_val = self.value-1
        
        while tmp_val > 0:
            remainder = (tmp_val-1)%26
            return_chars.append( self._char_mapping[remainder] )
            tmp_val =  (tmp_val + 1 - remainder)/26
        
        return "".join( reversed(return_chars) )
    
    def increment( self ):
        self.value += 1

"""
def is_blacklisted( region, bl_regions ):
    start, stop = region
    for bl_start, bl_stop in bl_regions:
        if not( stop < bl_start or start > bl_stop ):
            return True
    
    return False
"""

def main( ann_fname, ref_fname, good_gene_merges_fname, 
          bl_regions_fname, tss_exons_fname ):
    ann_ofp = open( "filtered.renamed.transcripts.gtf", "w" )

    novel_gene_names_set = set()
    NOVEL_GENE_CNTR = 0
    NOVEL_TRANS_CNTR = defaultdict( HexavijezimalCntr )
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
    # bl_regions = get_blacklist_regions( bl_regions_fname, ref_genes )
    blacklist_genes = get_blacklist_genes( bl_regions_fname )
    intron_to_gene_name, trans_to_ref_name, cds_bnds_to_gene_name \
        = get_introns_by_gene_name(ref_genes)
    
    # keep track of the genes that we've seen exactly
    observed_reference_trans = set()
    
    trans_id_to_names = {}
    for ( gene_name, chrm, strand, start, stop, transcripts ) in novel_genes:
        num_trans = 0

        # check to see if this overlaps a black list region
        # if is_blacklisted( (start, stop), bl_regions[(chrm, strand)] ):
        #    continue
        
        for trans in transcripts:
            fully_assoc_genes, assoc_genes = find_associated_genes(
                trans, intron_to_gene_name[(chrm, strand)])
            
            # if this intersects with any blacklist gene, then
            # skip this transcript
            if any( x in blacklist_genes for x in assoc_genes ):
                continue
            
            # if this is a completely novel transcript
            if len( assoc_genes ) == 0:
                if gene_name not in novel_gene_names_set:
                    NOVEL_GENE_CNTR += 1
                    novel_gene_names_set.add( gene_name )
                
                new_gene_name = "mgn%s" % str(NOVEL_GENE_CNTR).zfill( 5 )
                if new_gene_name not in NOVEL_TRANS_CNTR:
                    NOVEL_TRANS_CNTR[ new_gene_name ].increment()
                NOVEL_TRANS_CNTR[ new_gene_name ].increment()
                
                new_trans_name = new_gene_name + '.' + \
                    str( NOVEL_TRANS_CNTR[ new_gene_name ] )
                trans_id_to_names[trans.id] = ( new_gene_name, new_trans_name )
            # if there is only 1 corresponding ref gene
            elif len(fully_assoc_genes) > 0 or len( assoc_genes ) == 1:
                # check to see if this is an exact match
                if trans.IB_key() in trans_to_ref_name[(chrm, strand)]:
                    i_assoc_transripts = trans_to_ref_name[
                        (chrm, strand)][trans.IB_key()]
                    
                    i_trans_name, best_gene_name = \
                        next(iter(i_assoc_transripts))
                    
                    observed_reference_trans.update( 
                        x[0] for x in i_assoc_transripts )
                    
                    if IDENT_CNTR[ i_trans_name ].value > 0:
                        new_i_trans_name = i_trans_name + ".%s" \
                            % (IDENT_CNTR[ i_trans_name ])
                    else:
                        IDENT_CNTR[ i_trans_name ].increment()
                        new_i_trans_name = i_trans_name
                    
                    IDENT_CNTR[ i_trans_name ].increment()
                    
                    trans_id_to_names[ trans.id ] = ( 
                        best_gene_name, new_i_trans_name )
                # if it's not, it must be a partial match
                else:
                    # try to find the best gene name
                    # first search in CDS matched genes
                    cds_gene_names = cds_bnds_to_gene_name[
                        (chrm, strand)][trans.cds_region]
                    if len( cds_gene_names ) > 0:
                        fully_assoc_cds_gene_names = cds_gene_names.intersection(
                            fully_assoc_genes)
                        if len( fully_assoc_cds_gene_names ) > 0:
                            best_gene_name = next(iter(fully_assoc_cds_gene_names))
                        else:
                            best_gene_name = next(iter(cds_gene_names))
                    # if there are no matching CDS genes, then prefer fully 
                    # matching genes, 
                    else:
                        if len( fully_assoc_genes ) > 0:
                            best_gene_name = fully_assoc_genes[0]
                        else:
                            best_gene_name = assoc_genes[0]
                
                    if best_gene_name not in ALTERNATE_SPLICE_CNTR:
                        ALTERNATE_SPLICE_CNTR[ best_gene_name ].increment()
                    ALTERNATE_SPLICE_CNTR[ best_gene_name ].increment()
                    i_trans_name = "%s.%s" % ( 
                        best_gene_name, ALTERNATE_SPLICE_CNTR[ best_gene_name ])
                    trans_id_to_names[ trans.id ] = ( 
                        best_gene_name, i_trans_name )
            # if this is a gene merge
            elif len( assoc_genes ) > 1:
                assert len( fully_assoc_genes ) == 0
                # make sure that it's white listed
                if assoc_genes not in good_gene_merges:
                    continue
                i_gene_name = "/".join( assoc_genes )
                ALTERNATE_SPLICE_CNTR[ i_gene_name ].increment()
                i_trans_name = "%s.%s" % ( 
                    i_gene_name, ALTERNATE_SPLICE_CNTR[ i_gene_name ] )
                trans_id_to_names[ trans.id ] = ( i_gene_name, i_trans_name )

            new_gene_name, new_trans_name = trans_id_to_names[ trans.id ]
            
            trans.id = new_trans_name
            gtf_lines = trans.build_gtf_lines( new_gene_name,
                {'gene_type': 'gene', 'transcript_type': 'mRNA'} )
            
            ann_ofp.write(  gtf_lines + "\n" )
        
    ann_ofp.close()
    
    tss_exons_mapping = defaultdict( set )
    with open( tss_exons_fname ) as fp:
        for line in fp:
            data = parse_gff_line( line )
            if data == None: continue
            rgn = data.region
            if data.region.strand == '+':
                tss_exons_mapping[(rgn.chr, rgn.strand, rgn.stop) ].add(
                    rgn.start )
            else:
                tss_exons_mapping[(rgn.chr, rgn.strand, rgn.start)].add(
                    rgn.stop )
    
    def fix_5p_bnd( tr ):
        if tr.strand == '+':
            key = (tr.chrm, tr.strand, tr.exons[0][1])
            bnd = tr.exons[0][0]
        else:
            key = (tr.chrm, tr.strand, tr.exons[-1][0])
            bnd = tr.exons[-1][1]
        
        # if this is a novel bnd, return the original transcript
        if key not in tss_exons_mapping:
            return tr
        
        # otherwise, get the previous boundaries
        if not tr.is_protein_coding:
            bnds = tr.exon_bnds
        else:
            bnds = list( chain( chain(*tr.us_exons), 
                                chain(*tr.cds_exons), 
                                chain(*tr.ds_exons) ) )
        if tr.strand == '+':
            other_bnds = [i for i in tss_exons_mapping[key] 
                              if i <= bnd and bnd - i < 500 ]
            if len( other_bnds ) == 0: return tr
            new_distal = min( other_bnds )
            if len(tr.us_exons) > 0:
                bnds[0] = new_distal
            else:
                bnds.insert(0, bnds[1]-1) 
                bnds.insert(0, new_distal )
        else:
            other_bnds = [ i for i in tss_exons_mapping[key] 
                              if i >= bnd and i - bnd < 500  ]
            if len( other_bnds ) == 0: return tr
            new_distal = max( other_bnds )
            if len(tr.ds_exons) > 0:
                bnds[-1] = new_distal
            else:
                bnds.append(bnds[-1]+1) 
                bnds.append(new_distal)
                
        return Transcript( tr.id, tr.chrm, tr.strand, bnds, tr.cds_region )

    unobserved_trans_ids = set()
    for item in trans_to_ref_name.values():
        for i_item in item.values():
            for trans_id, gene_id in i_item:
                if trans_id not in observed_reference_trans:
                    unobserved_trans_ids.add( trans_id )
    
    fb_ofp = open( "unobserved_reference.gtf", "w" )
    for gene in ref_genes:
        for trans in gene[-1]:
            if trans.id in unobserved_trans_ids:
                trans = fix_5p_bnd( trans )
                lines = trans.build_gtf_lines(gene[0], 
                    { 'gene_type': 'gene', 'transcript_type': 'mRNA' } )
                fb_ofp.write( lines + "\n" )
    
    fb_ofp.close()
    
    return

if __name__ == "__main__":
    main( sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5] )
