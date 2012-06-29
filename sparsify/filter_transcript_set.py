# Copyright (c) 2011-2012 Nathan Boley

import sys, os
from collections import defaultdict
from itertools import chain

sys.path.append( os.path.join( os.path.dirname(__file__), 
                               '..', 'file_types', 'fast_gtf_parser' ) )
from gtf import load_gtf, Transcript

VERBOSE = False

def merge_distal_bndries( transcripts ):
    def make_mapping( max_link_size, is_downstream ):
        bnd_index = -1 if is_downstream else 0
        
        # find all of the distal bnds
        bnds = []
        for trans in transcripts:
            bnds.append( trans.exon_bnds[bnd_index] )
        
        # cluster them using single linkage
        bnds = sorted( bnds, reverse=is_downstream )
        grpd_bnds = [ [bnds[0],], ]
        for bnd in bnds[1:]:
            if abs(grpd_bnds[-1][-1] - bnd) > max_link_size:
                grpd_bnds.append( [bnd,] )
            else:
                grpd_bnds[-1].append( bnd )
        
        # make a mapping from bnds to corrected bnds
        rv = {}
        for bnd_grp in grpd_bnds:
            for bnd in bnd_grp:
                rv[bnd] = bnd_grp[0]
        
        return rv
    
    us_mapping = make_mapping( 25, False )
    ds_mapping = make_mapping( 25, True )
    
    new_transcripts = set()
    for tr in transcripts:
        if tr.is_protein_coding:
            new_bnds = list ( chain( chain(*tr.us_exons),  
                                     chain(*tr.cds_exons),     
                                     chain(*tr.ds_exons) ) )
        else:
            new_bnds = tr.exon_bnds
        
        new_bnds[0] = us_mapping[ new_bnds[0] ]
        new_bnds[-1] = ds_mapping[ new_bnds[-1] ]
        
        new_transcripts.add( Transcript( 
                tr.id, tr.chrm, tr.strand, new_bnds, tr.cds_region )  )
    
    return sorted(new_transcripts)

def split_transcripts_by_CDS( gene, fl_transcripts_set ):    
    # build the unique set of elements
    uniq_fp_utrs = set()
    uniq_cds = set()
    uniq_tp_utrs = set()
    
    for transcript in gene[-1]:
        if not transcript.is_protein_coding: continue
        uniq_fp_utrs.add( transcript.fp_utr_exons )
        uniq_cds.add( transcript.cds_exons )
        uniq_tp_utrs.add( transcript.tp_utr_exons )
    
    final_transcripts = []
    trans_to_iter_through = []
    # add non coding transcripts
    for transcript in gene[-1]:
        if not transcript.is_protein_coding:
            final_transcripts.append( transcript )
        else:
            trans_to_iter_through.append( transcript )

    # add transcripts that have full length support
    unadded_trans = []
    for transcript in trans_to_iter_through:
        if transcript.introns in fl_transcripts_set:
            final_transcripts.append( transcript )
            uniq_fp_utrs.discard( transcript.fp_utr_exons )
            uniq_cds.discard( transcript.cds_exons )
            uniq_tp_utrs.discard( transcript.tp_utr_exons )
        else:
            unadded_trans.append( transcript )
    trans_to_iter_through = unadded_trans
    
    # remove 3' utr's that share a coding exon with another coding seqeunce
    cds_exons_set = set()
    for transcript in trans_to_iter_through:
        cds_exons_set.update( transcript.cds_exons )
    for tp_utr in list(uniq_tp_utrs):
        for tp_utr_exon in tp_utr:
            num_utr_exons = 0
            if tp_utr_exon in cds_exons_set:
                num_utr_exons += 1
            if num_utr_exons >= 1:
                uniq_tp_utrs.remove( tp_utr )
                break

    # remove 5' utr's that share a coding exon with another coding seqeunce
    for fp_utr in list(uniq_fp_utrs):
        for fp_utr_exon in fp_utr:
            num_utr_exons = 0
            if fp_utr_exon in cds_exons_set:
                num_utr_exons += 1
            if num_utr_exons >= 1:
                uniq_fp_utrs.remove( fp_utr )
                break
    
    def find_transcripts_that_match_elements( 
            trans_iter, fp_flag=True, cds_flag=True, tp_flag=True ):
        unadded_trans = []
        # add transcripts with unique everything
        for tr in trans_iter:
            if ( (not fp_flag) or tr.fp_utr_exons in uniq_fp_utrs ) \
                    and ( (not tp_flag) or tr.tp_utr_exons in uniq_tp_utrs ) \
                    and ( (not cds_flag) or tr.cds_exons in uniq_cds ):
                final_transcripts.append( tr )
                if fp_flag:  uniq_fp_utrs.remove( tr.fp_utr_exons )
                if cds_flag:  uniq_cds.remove( tr.cds_exons )
                if tp_flag: uniq_tp_utrs.remove( tr.tp_utr_exons )
            else:
                unadded_trans.append( tr )
        
        return unadded_trans
    
    trans_to_iter_through = find_transcripts_that_match_elements( 
        trans_to_iter_through, True, True, True )
    trans_to_iter_through = find_transcripts_that_match_elements( 
        trans_to_iter_through, True, True, False )
    trans_to_iter_through = find_transcripts_that_match_elements( 
        trans_to_iter_through, False, True, True )
    trans_to_iter_through = find_transcripts_that_match_elements( 
        trans_to_iter_through, False, True, False )
    trans_to_iter_through = find_transcripts_that_match_elements( 
        trans_to_iter_through, True, False, False )
    trans_to_iter_through = find_transcripts_that_match_elements( 
        trans_to_iter_through, False, False, True )
        
    return final_transcripts

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(\
        description='Find exons from wiggle and junctions files.')

    parser.add_argument( 'novel_genes', type=file, \
        help='GTF with transcripts to be filtered.')
    parser.add_argument( 'ref_genes', type=file, \
        help='GTF with reference transcripts.')
    parser.add_argument( 'cdnas', type=file, \
        help='GTF with full length cDNAs.')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
        help='Whether or not to print status information.')

    args = parser.parse_args()
    
    # set flag args
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.novel_genes, args.ref_genes, args.cdnas

def main():
    novel_genes_fp, ref_genes_fp, cdnas_fp = parse_arguments()
    
    # load the proteomics gtf file
    # transcripts should have: 5' UTR, CDS, 3'UTR
    novel_genes = load_gtf( novel_genes_fp )
        
    # build a set of the full length transcripts
    ref_genes = load_gtf( ref_genes_fp )
    
    # find the full length cDNAs
    cdna_genes = load_gtf( cdnas_fp )
    
    full_len_transcripts = defaultdict( set )
    for gene in chain( ref_genes, cdna_genes ):
        for tr in gene[-1]:
            full_len_transcripts[(tr.chrm, tr.strand)].add( tr.introns )
    
    for gene in novel_genes:
        gene[-1] = merge_distal_bndries(gene[-1])
        
        filtered_transcripts = split_transcripts_by_CDS( 
            gene, full_len_transcripts[(gene[1], gene[2])] )
        for tr in filtered_transcripts:
            print tr.build_gtf_lines( gene[0], {} )
    
main()
