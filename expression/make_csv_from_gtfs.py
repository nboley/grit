# Copyright (c) 2011-2012 Nathan Boley

import os, sys
from collections import defaultdict

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "..", "file_types" ) )
from gtf_file import parse_gff_line

def get_exons_and_genes_in_gff(gff_fname):
    fp = open( gff_fname )
    gene_expression_scores = {}
    exon_expression_scores = {}
    jn_expression_scores = defaultdict(int)
    for line in fp:
        data = parse_gff_line( line )
        if data == None: 
            continue
        if data.feature == 'gene':
            gene_expression_scores[data.group] = data.score
        elif data.feature == 'exon':
            exon_expression_scores[data.region] = data.score
        elif data.feature == 'intron':
            jn_expression_scores[data.region] = data.score
        else:
            raise ValueError, "UNRECOGNIZED FEATURE" + line.strip()
    
    fp.close()
    sample_name = ".".join(os.path.basename(gff_fname).split(".")[:-2])
    return sample_name, gene_expression_scores, \
        exon_expression_scores, jn_expression_scores

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(\
        description="Combine expression scores into csv's. Currently only uses"\
            + "exon, gene, and intron element types. If an output filename is"\
            + "not given, then that data type is not collected.")
    
    parser.add_argument( 'expression_gff', nargs="+", type=str, \
        help='GFF file with expression values.')

    parser.add_argument( '--gene-output-filename', \
        help='Filename in which the gene expression data should be written.')
    parser.add_argument( '--exon-output-filename', \
        help='Filename in which the exon expression data should be written.')
    parser.add_argument( '--intron-output-filename', \
        help='Filename in which the intron expression data should be written.')
    
    args = parser.parse_args()
    return args.expression_gff, args.gene_output_filename, \
        args.exon_output_filename, args.intron_output_filename

def main():
    exp_gff_fnames, g_opfname, e_opfname, i_opfname = parse_arguments()

    sample_names = []
    gene_names = set()
    ge_scores = []
    exon_regions = set()
    ee_scores = []
    jn_regions = set()
    jn_scores = []
    for fname in exp_gff_fnames:
        sample_name, gene_expression_scores, exon_expression_scores, \
            jn_expression_scores =  get_exons_and_genes_in_gff( fname )
        sample_names.append( sample_name )
        ge_scores.append( gene_expression_scores )
        gene_names.update( gene_expression_scores.iterkeys() )
        ee_scores.append( exon_expression_scores )
        exon_regions.update( exon_expression_scores.iterkeys() )
        jn_scores.append( jn_expression_scores )
        jn_regions.update( jn_expression_scores.iterkeys() )
    
    # write out the gene names 
    if g_opfname != None:
        ofp = open( g_opfname, "w" )
        ofp.write( "gene_name," + ",".join( sample_names ) + "\n" )
        for gene_name in sorted( gene_names ):
            op_line = [ '"%s"' % gene_name,]
            for x in ge_scores:
                op_line.append( "%e" % x[gene_name] ) 
            ofp.write( ",".join( op_line ) + "\n" )
        ofp.close()

    # write out the exon scores
    if e_opfname != None:
        ofp = open( e_opfname, "w" )
        ofp.write( "gene_name," + ",".join( sample_names ) + "\n" )
        for exon_region in sorted( exon_regions ):
            exon_region_str = "%s:%s:%i-%i" % exon_region
            op_line = [exon_region_str,]
            for x in ee_scores:
                op_line.append( "%e" % x[exon_region] ) 
            ofp.write( ",".join( op_line ) + "\n" )
        ofp.close()

    # write out the jn scores
    if i_opfname != None:
        ofp = open( i_opfname, "w" )
        ofp.write( "gene_name," + ",".join( sample_names ) + "\n" )
        for jn_region in sorted( jn_regions ):
            jn_region_str = "%s:%s:%i-%i" % jn_region
            op_line = [jn_region_str,]
            for x in jn_scores:
                op_line.append( "%i" % int(x[jn_region]) ) 
            ofp.write( ",".join( op_line ) + "\n" )
        ofp.close()


main()
