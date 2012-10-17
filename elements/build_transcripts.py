# Copyright (c) 2011-2012 Nathan Boley

MAX_NUM_TRANSCRIPTS = 20
MIN_INTRON_CNT_FRAC = 0.001
VERBOSE = False
MIN_VERBOSE = False

import sys
import os
import numpy

from collections import namedtuple
CategorizedExons = namedtuple( "CategorizedExons", ["TSS", "TES", "internal"] )

import igraph
from cluster_exons import \
    find_overlapping_exons, find_jn_connected_exons, find_paths


sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "./exons/" ) )
from build_genelets import cluster_exons, cluster_overlapping_exons

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../sparsify/" ) )
import transcripts as transcripts_module

sys.path.append( os.path.join(os.path.dirname(__file__), "../file_types") )
from exons_file import parse_exons_files
from junctions_file import parse_jn_gffs, Junctions

from itertools import product, izip, chain

from collections import defaultdict
import multiprocessing
import Queue

class multi_safe_file( file ):
    def __init__( self, *args ):
        args = list( args )
        args.insert( 0, self )
        file.__init__( *args )
        self.lock = multiprocessing.Lock()

    def write( self, line ):
        self.lock.acquire()
        file.write( self, line + "\n" )
        self.flush()
        self.lock.release()

def build_gtf_line( gene_name, chr_name, gene_strand, trans_name, exon_num, \
                        start, stop ):
    gtf_line = list()
    
    if not chr_name.startswith( 'chr' ):
        chr_name = 'chr' + chr_name
    gtf_line.append( chr_name )
    gtf_line.append( 'build_transcripts' )
    gtf_line.append( 'exon' )
    gtf_line.append( str(start))
    gtf_line.append( str(stop))
    gtf_line.append( '.' )
    gtf_line.append( gene_strand )
    gtf_line.append( '.' )
    gtf_line.append( 'gene_id "' + gene_name + \
                         '"; transcript_id "' + trans_name + \
                         '"; exon_number "' + str(exon_num) + '";' )
    return '\t'.join( gtf_line )

def build_gtf_lines( gene_name, chrm, gene_strand, transcript_name, exons ):
    lines = []
    for exon_num, ( start, stop ) in enumerate( exons ):
        lines.append( build_gtf_line(gene_name, chrm, gene_strand, \
                                        transcript_name, exon_num, start, stop))
        
    return lines

def cluster_exons( tss_exons, internal_exons, tes_exons, jns, strand ):
    assert isinstance( tss_exons, set )
    assert isinstance( internal_exons, set )
    assert isinstance( tes_exons, set )
    
    all_exons = sorted( chain(tss_exons, internal_exons, tes_exons) )
    
    genes_graph = igraph.Graph()
    genes_graph.add_vertices( xrange(len(all_exons)) )
    genes_graph.add_edges( find_overlapping_exons(all_exons) )
    genes_graph.add_edges( find_jn_connected_exons(all_exons, jns, strand ) )
    
    genes = genes_graph.components()
    for gene in genes:
        exons = [ all_exons[exon_i] for exon_i in gene ]
        yield CategorizedExons( tss_exons.intersection( exons ),
                                tes_exons.intersection( exons ),
                                internal_exons.intersection( exons ) )
    
    return

def build_transcripts( exons, jns, strand ):
    import networkx as nx
    # build a directed graph, with edges leading from exon to exon via junctions
    all_exons = sorted(chain(exons.TSS, exons.internal, exons.TES))
    graph = nx.DiGraph()
    graph.add_nodes_from( exons.TSS )
    graph.add_nodes_from( exons.internal )
    graph.add_nodes_from( exons.TES )
    
    edges = find_jn_connected_exons(all_exons, jns, strand, use_names=True )
    graph.add_edges_from( edges  )
    
    transcripts = []
    for tss in exons.TSS:
        for tes in exons.TES:
            for i,transcript in enumerate(nx.all_simple_paths(graph, tss, tes)):
                if i > MAX_NUM_TRANSCRIPTS: 
                    print >> sys.stderr, "TOO COMPLEX"
                    print >> sys.stderr, exons
                    return transcripts # []
                transcripts.append( transcript )
    
    return transcripts

def build_genes(all_internal_exons, all_tss_exons, all_tes_exons, jns, ofp):
    keys = set(chain(all_tss_exons.iterkeys(), 
                     all_internal_exons.iterkeys(), 
                     all_tes_exons.iterkeys()))
    gene_id = 0
    for (chrm, strand) in keys:
        tss_exons = set( map(tuple, all_tss_exons[(chrm, strand)].tolist()))
        internal_exons = set(map(tuple, all_internal_exons[(chrm, strand)].tolist()))
        tes_exons = set(map(tuple, all_tes_exons[(chrm, strand)].tolist()))
        for gene_exons in cluster_exons( 
                tss_exons, internal_exons, tes_exons, jns[(chrm, strand)], strand):
            transcripts = build_transcripts( 
                gene_exons, jns[(chrm, strand)], strand )
            gene_name = "GENE%i" % gene_id
            gene_id += 1
            for trans_id, exons in enumerate(transcripts):
                trans_name = gene_name + "_TRANS%i" % trans_id
                for line in build_gtf_lines( 
                    gene_name, chrm, strand, trans_name, exons):
                    ofp.write( line + "\n" )
                


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Determine valid transcripts.')
    parser.add_argument( '--exons', type=file, required=True, nargs="+", \
        help='GFF file(s) contaning exons.')
    parser.add_argument( '--tss', type=file, required=True, nargs="+", \
        help='Gff file(s) containing valid transcription stop site exons. ' )
    parser.add_argument( '--tes', type=file, required=True, nargs="+", \
        help='Gff file(s) containing valid transcription stop site exons.')
    parser.add_argument( '--junctions', type=file, required=True, nargs="+", \
        help='Gff file(s) containing valid introns.')
    parser.add_argument( '--single-exon-genes', type=file, nargs="*", 
        help='Single exon genes file.')    

    parser.add_argument( '--threads', '-t', type=int , default=1, \
                             help='Number of threads spawn for multithreading.')
    parser.add_argument( '--out-fname', '-o', default='',\
                             help='Output file name. (default: stdout)')
    parser.add_argument( '--log-fname', '-l',\
                             help='Output log file name. (default: sterr)')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
                             help='Whether or not to print status information.')
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    transcripts_module.VERBOSE = VERBOSE
    
    out_fp = open( args.out_fname, "w" ) if args.out_fname else sys.stdout
    
    log_fp = multi_safe_file( args.log_fname, 'w' ) \
        if args.log_fname else sys.stderr
    
    return args.exons, args.tss, args.tes, args.junctions, \
        args.single_exon_genes, args.threads, out_fp, log_fp

def load_jns( jn_gff_fps, min_cnt_ratio=0.1 ):
    all_jns = parse_jn_gffs( jn_gff_fps )
    donor_scores = defaultdict(lambda: defaultdict(int))
    acceptor_scores = defaultdict(lambda: defaultdict(int))
    for jn in all_jns:
        donor_scores[(jn.chrm, jn.strand)][jn.start] \
            = max( jn.cnt, donor_scores[(jn.chrm, jn.strand)][jn.start] )
        acceptor_scores[(jn.chrm, jn.strand)][jn.stop] \
            = max( jn.cnt, acceptor_scores[(jn.chrm, jn.strand)][jn.stop] )

    jns = defaultdict(list)    
    for jn in all_jns:
        donor_score = donor_scores[(jn.chrm, jn.strand)][jn.start]
        if float(jn.cnt)/donor_score < min_cnt_ratio:
            continue
        acceptor_score = acceptor_scores[(jn.chrm, jn.strand)][jn.stop]
        if float(jn.cnt)/acceptor_score < min_cnt_ratio:
            continue

        jns[(jn.chrm, jn.strand)].append( (jn.start, jn.stop) )
    
    return dict((key, numpy.array(val)) for key, val in jns.iteritems() )

def main():
    exon_fps, tss_exon_fps, tes_exon_fps, junction_fps, \
        single_exon_fps, threads, out_fp, log_fp = parse_arguments()
    
    if VERBOSE:
        print >> sys.stderr, 'Parsing input...'
    
    
    internal_exons = parse_exons_files( exon_fps )
    tss_exons = parse_exons_files( tss_exon_fps )
    tes_exons = parse_exons_files( tes_exon_fps )
    jns = load_jns( junction_fps )
    
    if None != single_exon_fps:
        se_genes = parse_exons_files( single_exon_fps )
        gene_id = 1
        for (chrm, strand), exons in se_genes.iteritems():
            for start, stop in exons:
                gene_name = "single_exon_gene_%i" % gene_id
                out_fp.write( build_gtf_line( 
                    gene_name, chrm, strand, gene_name, 1, start, stop ) + "\n")
                gene_id += 1
    
    build_genes(internal_exons, tss_exons, tes_exons, jns, out_fp)
            
    return

if __name__ == "__main__":
    main()
