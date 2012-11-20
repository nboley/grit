# Copyright (c) 2011-2012 Nathan Boley

import sys, os
from itertools import izip, combinations, chain
from collections import defaultdict, namedtuple
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
            equiv_key = trans.IB_key()
            trans_to_ref_name[(chrm, strand)][equiv_key].add(
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

def get_good_gene_merges( fp ):
    good_gene_merges = set()
    for line in fp:
        gene_names = tuple(sorted(line.split()))
        good_gene_merges.add( gene_names )
        # add all of the combinations
        for i in xrange( 2, len(gene_names) ):
            for new_gene_names in combinations(gene_names, i):
                good_gene_merges.add( tuple(sorted(new_gene_names)) )
    
    return good_gene_merges

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

class MutableInt(object):
    def __init__( self, value=0 ):
        self.value = value
    def increment( self ):
        self.value += 1
    def __str__( self ):
        return str( self.value )
    
Cntrs = namedtuple( "Cntrs", (
        'NOVEL_GENE', 'NOVEL_TRANS', 'IDENT', 'GENE_MERGE', 'ALTERNATE_SPLICE'))
def init_counters():
    NOVEL_GENE = MutableInt(0)
    NOVEL_TRANS = defaultdict( HexavijezimalCntr )
    # even if the transcript match is identical, we may have alternate 5' and 3'
    # ends which means multiple 'identical' transcripts. So keep a counter. If 
    # the counter is at 0 *for that ref trans id*, keep the ref name. Otherwise,
    # increment it by one. 
    IDENT = defaultdict( HexavijezimalCntr )
    # similar for gene merges, but index by the tpule of gene names
    GENE_MERGE = defaultdict( HexavijezimalCntr )
    # for novel splices, againb, index by gene name
    ALTERNATE_SPLICE = defaultdict( HexavijezimalCntr )
    
    return Cntrs( NOVEL_GENE, NOVEL_TRANS, IDENT, GENE_MERGE, ALTERNATE_SPLICE)

def build_gtf_lines_for_gene( gene, cntrs, 
                              intron_to_gene_name, trans_to_ref_name, 
                              cds_bnds_to_gene_name, observed_reference_trans,
                              novel_gene_names_set,
                              good_gene_merges, blacklist_genes):
    lines = []
    num_trans = 0
    trans_id_to_names = {}
    
    for trans in gene.transcripts:
        fully_assoc_genes, assoc_genes = find_associated_genes(
            trans, intron_to_gene_name[(gene.chrm, gene.strand)])

        # if this intersects with any blacklist gene, then
        # skip this transcript
        if blacklist_genes != None \
                and any( x in blacklist_genes for x in assoc_genes ):
            continue

        # if this is a completely novel transcript
        if len( assoc_genes ) == 0:
            if gene.id not in novel_gene_names_set:
                cntrs.NOVEL_GENE.increment()
                novel_gene_names_set.add( gene.id )

            new_gene_name = "mgn%s" % str(cntrs.NOVEL_GENE).zfill( 5 )
            if new_gene_name not in cntrs.NOVEL_TRANS:
                cntrs.NOVEL_TRANS[ new_gene_name ].increment()
            cntrs.NOVEL_TRANS[ new_gene_name ].increment()

            new_trans_name = new_gene_name + '.' + \
                str( cntrs.NOVEL_TRANS[ new_gene_name ] )
            trans_id_to_names[trans.id] = ( new_gene_name, new_trans_name )
        # if there is only 1 corresponding ref gene
        elif len(fully_assoc_genes) > 0 or len( assoc_genes ) == 1:
            # check to see if this is an exact match
            if trans.IB_key() in trans_to_ref_name[(gene.chrm, gene.strand)]:
                i_assoc_transripts = trans_to_ref_name[
                    (gene.chrm, gene.strand)][trans.IB_key()]

                i_trans_name, best_gene_name = \
                    next(iter(i_assoc_transripts))

                observed_reference_trans.update( 
                    x[0] for x in i_assoc_transripts )

                if cntrs.IDENT[ i_trans_name ].value > 0:
                    new_i_trans_name = i_trans_name + ".%s" \
                        % (cntrs.IDENT[ i_trans_name ])
                else:
                    cntrs.IDENT[ i_trans_name ].increment()
                    new_i_trans_name = i_trans_name

                cntrs.IDENT[ i_trans_name ].increment()

                trans_id_to_names[ trans.id ] = ( 
                    best_gene_name, new_i_trans_name )
            # if it's not, it must be a partial match
            else:
                # try to find the best gene name
                # first search in CDS matched genes
                cds_gene_names = cds_bnds_to_gene_name[
                    (gene.chrm, gene.strand)][trans.cds_region]
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

                if best_gene_name not in cntrs.ALTERNATE_SPLICE:
                    cntrs.ALTERNATE_SPLICE[ best_gene_name ].increment()
                cntrs.ALTERNATE_SPLICE[ best_gene_name ].increment()
                i_trans_name = "%s.%s" % ( 
                    best_gene_name, cntrs.ALTERNATE_SPLICE[ best_gene_name ])
                trans_id_to_names[ trans.id ] = ( 
                    best_gene_name, i_trans_name )
        # if this is a gene merge
        elif len( assoc_genes ) > 1:
            assert len( fully_assoc_genes ) == 0
            # make sure that it's white listed
            if good_gene_merges != None \
                    and assoc_genes not in good_gene_merges:
                continue
            i_gene_name = "/".join( assoc_genes )
            cntrs.ALTERNATE_SPLICE[ i_gene_name ].increment()
            i_trans_name = "%s.%s" % ( 
                i_gene_name, cntrs.ALTERNATE_SPLICE[ i_gene_name ] )
            trans_id_to_names[ trans.id ] = ( i_gene_name, i_trans_name )

        new_gene_name, new_trans_name = trans_id_to_names[ trans.id ]

        trans.id = new_trans_name
        gtf_lines = trans.build_gtf_lines( new_gene_name,
            {'gene_type': 'gene', 'transcript_type': 'mRNA'}, source='grit' )

        lines.append( gtf_lines )
    
    return lines 

def filter_and_rename_transcripts( ref_genes, novel_genes, ann_ofp,
                                   good_gene_merges=None, blacklist_genes=None):
    cntrs = init_counters()

    novel_gene_names_set = set()    
    intron_to_gene_name, trans_to_ref_name, cds_bnds_to_gene_name \
        = get_introns_by_gene_name(ref_genes)
    
    # keep track of the genes that we've seen exactly
    observed_reference_trans = set()
    
    for gene in novel_genes:
        lines = build_gtf_lines_for_gene( 
            gene, cntrs, intron_to_gene_name, trans_to_ref_name, 
            cds_bnds_to_gene_name, observed_reference_trans, 
            novel_gene_names_set,
            good_gene_merges, blacklist_genes )
        ann_ofp.write( "\n".join( lines ) + '\n' )
    
    return observed_reference_trans

def build_tss_exons_mapping( tss_exons_fp ):
    tss_exons_mapping = defaultdict( set )
    for line in tss_exons_fp:
        data = parse_gff_line( line )
        if data == None: continue
        rgn = data.region
        if data.region.strand == '+':
            tss_exons_mapping[(rgn.chr, rgn.strand, rgn.stop) ].add(
                rgn.start )
        else:
            tss_exons_mapping[(rgn.chr, rgn.strand, rgn.start)].add(
                rgn.stop )
    
    return tss_exons_mapping

def fix_5p_bnd( tr, tss_exons_mapping ):
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
        if tr.is_protein_coding and len(tr.us_exons) > 0:
            bnds[0] = new_distal
        else:
            bnds.insert(0, bnds[1]-1) 
            bnds.insert(0, new_distal )
    else:
        other_bnds = [ i for i in tss_exons_mapping[key] 
                          if i >= bnd and i - bnd < 500  ]
        if len( other_bnds ) == 0: return tr
        new_distal = max( other_bnds )
        if tr.is_protein_coding and len(tr.ds_exons) > 0:
            bnds[-1] = new_distal
        else:
            bnds.append(bnds[-1]+1) 
            bnds.append(new_distal)

    return Transcript( tr.id, tr.chrm, tr.strand, bnds, tr.cds_region )

def parse_arguments():
    global MIN_AAS_PER_ORF
    global OUTPUT_PROTEINS
    global VERBOSE
    import argparse
    
    parser = argparse.ArgumentParser(
        description = 'Rename a transcript file with respect to a reference, '\
        + 'and perform reference guided filtering.')
    
    parser.add_argument(
        'annotation', type=file,
        help='GTF containing the annotation to be renamed.' )
    parser.add_argument(
        'reference', type=file,
        help='GTF file containing the reference annotation.' )

    parser.add_argument(
        '--output-file', 
        help='File to write renamed and filtered transcripts to.' )
    
    parser.add_argument(
        '--good-gene-merges', type=file,
        help='A file containing allowed gene merges WRT the reference. ' \
            + ' default: Allow all gene merges.' )
    parser.add_argument(
        '--blacklist-genes', type=file,
        help='Remove any transcripts in the annotation that overlap these genes.')
    parser.add_argument(
        '--tss-exons', type=file,
        help='A gff containing tss exons. When merging back in missed' \
            + ' transcripts, if this file contains a tss exon with the same' \
            + ' internal boundary then we will use that instead of the ref.')
        
    args = parser.parse_args()
    
    ofp = open( args.output_file, "w") \
        if args.output_file != None else sys.stdout
    
    return args.annotation, args.reference, ofp, \
        args.good_gene_merges, args.blacklist_genes, \
        args.tss_exons

def main():
    ann_fp, ref_fp, ofp, \
        ggm_fp, blg_fp, tss_exons_fp = parse_arguments()
        
    ggm = get_good_gene_merges( ggm_fp ) if ggm_fp != None else None
    bl_genes = get_blacklist_genes(blg_fp) if blg_fp != None else None
    novel_genes, ref_genes =  load_gtfs( ( ann_fp.name, ref_fp.name ), 2 )

    observed_reference_trans = filter_and_rename_transcripts( 
        ref_genes, novel_genes, ofp, ggm, bl_genes )
    
    tss_exons_mapping = build_tss_exons_mapping( tss_exons_fp ) \
        if tss_exons_fp != None else set()
    
    unobserved_trans_ids = set()
    for gene in ref_genes:
        for trans in gene.transcripts:
            if trans.id not in observed_reference_trans:
                unobserved_trans_ids.add( trans.id )
    
    """
    for gene in ref_genes:
        for trans in gene[-1]:
            if trans.id in unobserved_trans_ids:
                trans = fix_5p_bnd( trans, tss_exons_mapping )
                lines = trans.build_gtf_lines(gene[0], 
                    { 'gene_type': 'gene', 'transcript_type': 'mRNA' },
                                              source="reference")
                ofp.write( lines + "\n" )
    """
    
    ofp.close()
    
    return

if __name__ == "__main__":
    main()
