import sys, os

import numpy
import pickle

import matplotlib.pyplot as plt

from collections import defaultdict

sys.path.append( os.path.join( os.path.dirname(__file__), "..", "..", "sparsify" ) )
from sparsify_transcripts import get_raw_transcripts

sys.path.append( os.path.join( os.path.dirname(__file__), "..", "..", "file_types" ) )
from wiggle import Wiggle

sys.path.append( os.path.join( os.path.dirname(__file__), "..", "distal_exons" ) )
from distal_exons import get_exon_cvrg_scores, get_exons_w_quantls

# should be in ['region', 'total', 'avg', 'quantl']
EXON_FILTER_TYPE = 'quantl'
SMOOTHING_WINDOW_SIZE = 20
SMOOTH_CVRG_FOR_QUANTL = True
REGION_FILTER_LEN = 25

# should be in [ 'ROC', 'cumsum', 'rev_cumsum' ]
PLOT_TYPE = 'ROC'
SHOW_ROC_TABLE = False

SHOW_EMPTY_GENES = False

def make_roc_plot_data( all_scores ):
    # calc ROC curve data
    true_tss_exon_scores = numpy.array( sorted(all_scores[0]) )
    
    # get all negatives with internal and tes exons
    false_tss_exon_scores = all_scores[1]
    false_tss_exon_scores.extend( all_scores[2] )
    false_tss_exon_scores = numpy.array( sorted( false_tss_exon_scores ) )
    
    total_positives = len( true_tss_exon_scores )
    total_negative = len( false_tss_exon_scores )
    
    num_empty_tss_exons = true_tss_exon_scores.searchsorted( 0, side='right' )
    print 'Number of empty tss exons: ' + str( num_empty_tss_exons )
    
    x_coords = []
    y_coords = []
    cutoff_val = 0.0
    while cutoff_val <= 1.00:
        FN = true_tss_exon_scores.searchsorted( cutoff_val, side='right' )
        TN = false_tss_exon_scores.searchsorted( cutoff_val, side='right' )
        TP = total_positives - FN
        FP = total_negative - TN
        
        sn = float( TP ) / total_positives
        sp = float( TN ) / total_negative
        if SHOW_ROC_TABLE:
            print str(1.0-sp), str(sn), cutoff_val
        
        x_coords.append( 1.0 - sp )
        y_coords.append( sn )
        
        cutoff_val += 0.001
    
    return x_coords, y_coords

def plot_scores( all_scores ):
    if PLOT_TYPE in ( 'cumsum', 'rev_cumsum' ):
        # flip cumsum
        if PLOT_TYPE == 'rev_cumsum':
            new_all_scores = [[],[],[]]
            for i, scores in enumerate(all_scores):
                for score in scores:
                    new_all_scores[i].append( 1.0-score )
            all_scores = new_all_scores
        
        plt.hist( all_scores, bins=1001, range=(0,1), histtype='step', 
                  cumulative=True )
        # plt.axis( [0,1,0,5000] )
        plt.xlabel( 'TSS soverage score' )
        plt.ylabel( 'Num Exons' )
        plt.legend( ['TSS', 'Internal', 'TES'] )
        plt.title( 'TSS coverage (TSS, internal and TES exons)' )
        plt.show()
    
    else:
        assert PLOT_TYPE == 'ROC'
        x_coords, y_coords = make_roc_plot_data( all_scores )
        
        plt.plot( x_coords, y_coords )
        plt.xlabel( 'FPR (1 - Specificity)' )
        plt.ylabel( 'TRP or Sensitivity' )
        plt.title( 'ROC curve for TSS selection cutoff_val' )
        plt.show()
    
    return

def get_exon_scores( gtf_fp, tss_cvrg_fps, chrm_sizes_fp, out_fp ):
    print 'parsing transcripts file...'
    all_trans = get_raw_transcripts( gtf_fp )
    
    print 'building tss coverage...'
    cvrg = Wiggle( chrm_sizes_fp, tss_cvrg_fps )
    chrm_sizes_fp.close()
    
    print 'separating exons...'
    # create dict with exon coordinates of tss, internal and tes exons
    gene_exons = defaultdict( lambda: [set(),set(),set()] )
    for gene_name, transcripts in all_trans.iteritems():
        for trans_name, exons in transcripts.iteritems():
            # skip single exon genes
            if len(exons) == 1: continue
            
            for i, exon in enumerate(exons):
                if (i == 0 and exon.strand == '+') or \
                        (i == len(exons)-1 and exon.strand == '-'):
                    gene_exons[ ( exon.chr, exon.strand, gene_name ) ][0].add(
                        (exon.start, exon.stop) )
                elif (i == len(exons)-1 and exon.strand == '+') or \
                        (i == 0 and exon.strand == '-'):
                    gene_exons[ ( exon.chr, exon.strand, gene_name ) ][2].add(
                        (exon.start, exon.stop) )
                else:
                    gene_exons[ ( exon.chr, exon.strand, gene_name ) ][1].add(
                        (exon.start, exon.stop) )
    
    print 'getting exon scores...'
    # get scores by using EXON_FILTER_TYPE
    tss_exon_scores = []
    int_exon_scores = []
    tes_exon_scores = []
    zero_tss_cvrg_genes = []
    
    cvrg_for_quantl = cvrg.get_smoothed_copy( SMOOTHING_WINDOW_SIZE ) \
        if EXON_FILTER_TYPE == 'quantl' and SMOOTH_CVRG_FOR_QUANTL \
        else cage_read_cvrg
    
    for (chrm, strand, gene_name), exons in gene_exons.iteritems():
        # get all the exons together for each cluster
        all_gene_exons = sorted([
                exon for grp_exons in exons for exon in grp_exons])
        
        if EXON_FILTER_TYPE in [ 'quantl' ]:
            try:
                exons_w_quantls = get_exons_w_quantls( 
                    all_gene_exons, cvrg_for_quantl[(chrm, strand)] )
            except ValueError:
                zero_tss_cvrg_genes.append( gene_name )
                continue
            
            for exon, quantl in exons_w_quantls:
                if exon in exons[0]:
                    tss_exon_scores.append( quantl )
                elif exon in exons[1]:
                    int_exon_scores.append( quantl )
                else:
                    assert exon in exons[2]
                    tes_exon_scores.append( quantl )
        
        else:
            assert EXON_FILTER_TYPE in [ 'region', 'total', 'avg' ]
            all_gene_exons_w_score, max_gene_score, total_gene_score = \
                get_exon_cvrg_scores( all_gene_exons, cvrg[(chrm, strand)], 
                                      EXON_FILTER_TYPE, REGION_FILTER_LEN )
            if total_gene_score == 0:
                zero_tss_cvrg_genes.append( gene_name )
                continue
            
            for exon, score in all_gene_exons_w_score:
                if exon in exons[0]:
                    tss_exon_scores.append( float( score ) / max_gene_score )
                elif exon in exons[1]:
                    int_exon_scores.append( float( score ) / max_gene_score )
                else:
                    assert exon in exons[2]
                    tes_exon_scores.append( float( score ) / max_gene_score )
    
    if SHOW_EMPTY_GENES:
        print zero_tss_cvrg_genes
    print 'Number of genes with no coverage: ' + str(len(zero_tss_cvrg_genes))
    
    # pack scores
    all_scores = (tss_exon_scores, int_exon_scores, tes_exon_scores)
    if out_fp != None:
        print 'dumping exons...'
        pickle.dump( all_scores, out_fp )
        out_fp.close()
    
    return all_scores

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(
        description='Get transcript start sites(TSS) from a variety of sources.')
    parser.add_argument(
        '--transcripts', '-t',type=file,
        help='GTF file which contains transcripts.')
    parser.add_argument(
        '--chrm_sizes_fname', '-c', type=file,
        help='File with chromosome names and sizes.')
    parser.add_argument(
        '--tss_cvrg_wigs', '-w', type=file, nargs='*',
        help='WIG file with RNAseq read coverage. ' +
        'Note: strand will be infered from filename (plus/minus).')
    
    parser.add_argument( '--out_fname', '-o', type=argparse.FileType( 'w' ), \
        help='Output file to dump tss and internal exon cvrg.')
    parser.add_argument( '--cached_exon_counts', '-p', type=file,\
        help='Previously cached out_fname from this script.')
    args = parser.parse_args()
    
    if args.cached_exon_counts == None:
        if not all( item != None for item in [ 
                args.transcripts, args.chrm_sizes_fname, args.tss_cvrg_wigs ] ):
            parser.print_help()
            sys.exit()
    
    return ( args.transcripts, args.tss_cvrg_wigs, args.chrm_sizes_fname, 
             args.out_fname, args.cached_exon_counts )

def main():
    gtf_fp, tss_cvrg_fps, chrm_sizes_fp, out_fp, pickled_fp = parse_arguments()
    
    if pickled_fp == None:
        all_scores = get_exon_scores( gtf_fp, tss_cvrg_fps,
                                      chrm_sizes_fp, out_fp )
    else:
        all_scores = pickle.load( pickled_fp )
    
    plot_scores( all_scores )
    
    return

if __name__ == '__main__':
    main()

