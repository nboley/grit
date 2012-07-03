# Copyright (c) 2011-2012 Nathan Boley

import sys, os

from collections import defaultdict
from itertools import izip, chain

import numpy
from numpy import mean
from numpy.linalg import lstsq
from scipy.optimize import nnls

sys.path.append( os.path.join( os.path.dirname(__file__), 
                               '..', 'file_types', 'fast_gtf_parser' ) )
from gtf import load_gtf

sys.path.append(os.path.join(os.path.dirname(__file__), "..", 'file_types'))
from wiggle import Wiggle
from gtf_file import create_gff_line, GenomicInterval

def group_exons_by_gene( genes ):
    grpd_exons = defaultdict(set)
    tss_exons = set()
    tes_exons = set()
    for gene_name, chrm, strand, gene_start, gene_stop, transcripts in genes:
        for tran_name, transcript in transcripts:
            key = (gene_name, chrm, strand, gene_start, gene_stop)
            grpd_exons[ key ].update( 
                izip( transcript[:-1:2], transcript[1::2] ) )
            upstream_exon = (transcript[0], transcript[1])
            downstream_exon = (transcript[-1], transcript[-2])
            if strand == '+':
                tss_exons.add( upstream_exon )
                tes_exons.add( downstream_exon )
            else:
                tss_exons.add( downstream_exon )
                tes_exons.add( upstream_exon )
    
    return grpd_exons, tss_exons, tes_exons

def estimate_exon_expression( exons, read_cov, num_mapped_bases ):
    if num_mapped_bases == 0:
        return [0]*len( exons )
    
    starts = set( exon[0] for exon in exons )
    bndries = sorted(set(chain(*exons)))
    
    # build a mapping from pseudo cnt to average base coverage
    pseudo_cnts = {}
    for start, stop in izip(bndries[:-1], bndries[1:]):
        if stop in starts:
            stop -= 1
        pseudo_cnts[start] = 1000*read_cov[start:stop+1].mean()
    
    # build the design matrix. That is, for each exon, find which pseudo
    # exons it overlaps
    design_matrix = []
    for exon_start, exon_stop in exons:
        design_matrix.append( [] )
        for start in bndries[:-1]:
            if start < exon_start or start >= exon_stop:
                design_matrix[-1].append( 0 )
            else:
                design_matrix[-1].append( 1 )
    
    X = numpy.matrix( design_matrix )
    
    y = []
    for bndry in bndries[:-1]:
        y.append( pseudo_cnts[bndry] )
    y = numpy.array( y ).T
    
    # print y
    # print X
    # y = x*Beta
    # choose beta to minimize (y - X*Beta)^2
    # print nnls( X.T, y )
    coefs = nnls( X.T, y )[0]
    # print coefs
    # coefs = lstsq( X.T, y )[0]
    # turn the estiamtes into rpkm
    return [ int(x)/num_mapped_bases for x in coefs ]

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(\
        description='Estimate expression scores.')

    parser.add_argument( 'transcripts_gtf', type=file, \
        help='File with transcripts to estiamte expression for.')
    parser.add_argument( 'chrm_sizes_fname', type=file, \
        help='File with chromosome names and sizes.')
    parser.add_argument( 'rnaseq_wigs', type=file, nargs="+", \
        help='Read coverage wiggle files ( bedgraphs parse the fastest ).')
    
    parser.add_argument( '--cage-wigs', type=file, nargs='+', \
        help='wig files with cage reads, to further quantify TSS exons.')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
        help='Whether or not to print status information.')
    
    args = parser.parse_args()

    
    # set flag args
    global VERBOSE
    VERBOSE = args.verbose

    
    return args.transcripts_gtf.name, args.rnaseq_wigs, args.chrm_sizes_fname, args.cage_wigs

def main():
    transcripts_gtf_fname, wigs, chrm_sizes_fp, cage_wigs = parse_arguments()

    transcripts = load_gtf( transcripts_gtf_fname )
    grpd_exons, tss_exons, tes_exons = group_exons_by_gene( transcripts )

    read_cov = Wiggle( chrm_sizes_fp, wigs )    
    num_mapped_bases = sum( float(x.sum())/1e6 for x in read_cov.values() )
    
    for (gene_name, chrm, strand, gene_start, gene_stop), exons \
            in grpd_exons.iteritems():
        exon_num = 0
        exp_scores = estimate_exon_expression( 
            sorted(exons), read_cov[(chrm, strand)], num_mapped_bases )
        
        region = GenomicInterval( chrm, strand, gene_start, gene_stop )
        print create_gff_line( region, gene_name, 
                               mean(exp_scores), feature='gene', source='grit' )

        for (start, stop), score in izip( exons, exp_scores ):
            exon_num += 1
            region = GenomicInterval( chrm, strand, start, stop )
            print create_gff_line( region, "%s_%i" % (gene_name, exon_num), 
                                   score, feature='exon', source='grit' )
    
    return
    
    
    if cage_wigs != None:
        cage_cov = Wiggle( chrm_sizes_fp, cage_wigs )
    
    
    pass

main()
