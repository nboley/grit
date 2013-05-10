# Copyright (c) 2011-2012 Nathan Boley

MAX_NUM_TRANSCRIPTS = 100000
VERBOSE = False
MIN_VERBOSE = False

import sys
import os
import numpy
import time
import tempfile

from collections import namedtuple
CategorizedExons = namedtuple( "CategorizedExons", ["TSS", "TES", "internal", 
                                                    "full_transcript"] )

import igraph
from cluster_exons import \
    find_overlapping_exons, find_jn_connected_exons, find_paths

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../sparsify/" ) )
import transcripts as transcripts_module

sys.path.append( os.path.join(os.path.dirname(__file__), "../file_types") )
from exons_file import parse_exons_files
from junctions_file import parse_jn_gffs, Junctions

from itertools import product, izip, chain

from collections import defaultdict
import multiprocessing
import multiprocessing.sharedctypes
import Queue

num_threads = 1

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

def cluster_exons( tss_exons, internal_exons, tes_exons, se_transcripts, 
                   jns, strand ):
    assert isinstance( tss_exons, set )
    assert isinstance( internal_exons, set )
    assert isinstance( tes_exons, set )
    assert isinstance( se_transcripts, set )
    
    all_exons = sorted( chain(tss_exons, internal_exons, 
                              tes_exons, se_transcripts) )
    
    genes_graph = igraph.Graph()
    genes_graph.add_vertices( xrange(len(all_exons)) )
    genes_graph.add_edges( find_overlapping_exons(all_exons) )
    genes_graph.add_edges( find_jn_connected_exons(all_exons, jns, strand ) )
    
    genes = genes_graph.components()
    for gene in genes:
        exons = [ all_exons[exon_i] for exon_i in gene ]
        yield CategorizedExons( tss_exons.intersection( exons ),
                                tes_exons.intersection( exons ),
                                internal_exons.intersection( exons ),
                                se_transcripts.intersection( exons ) )
    
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
    num_transcripts = 0
    for tss in exons.TSS:
        for tes in exons.TES:
            for transcript in nx.all_simple_paths(graph, tss, tes):
                num_transcripts += 1
                if num_transcripts > MAX_NUM_TRANSCRIPTS: 
                    print >> sys.stderr, "TOO COMPLEX"
                    return []

                transcripts.append( transcript )
    
    return transcripts

def build_genes_in_contig( (chrm, strand), (gene_id, gene_id_lock), ofp,
                           tss_exons, internal_exons, tes_exons, 
                           se_transcripts, jns ):
    for gene_exons in cluster_exons( 
            tss_exons, internal_exons, tes_exons, se_transcripts, jns, strand):
        transcripts = build_transcripts( 
            gene_exons, jns, strand )
        gene_id_lock.acquire()
        gene_name = "GENE%i" % gene_id.value
        gene_id.value += 1
        gene_id_lock.release()
        for trans_id, exons in enumerate(transcripts):
            trans_name = gene_name + "_TRANS%i" % trans_id
            for line in build_gtf_lines( 
                gene_name, chrm, strand, trans_name, exons):
                ofp.write( line + "\n" )
    
    ofp.flush()
    return ofp

def build_genes(all_internal_exons, all_tss_exons, all_tes_exons, 
                all_se_transcripts, all_jns, ofp):
    keys = sorted( set(chain(all_tss_exons.iterkeys(), 
                             all_internal_exons.iterkeys(), 
                             all_tes_exons.iterkeys(),
                             all_se_transcripts)) )
    gene_id = multiprocessing.sharedctypes.Value( "i", 0 )
    gene_id_lock = multiprocessing.Lock()
    
    ps = [None]*num_threads
    ofps = {}
    while len( keys ) > 0:
        finished_i = None
        for i, p in enumerate(ps):
            if p == None or p.exitcode != None:
                finished_i = i
                break
        if finished_i == None:
            time.sleep(1)
            continue
        
        chrm, strand = keys.pop()
        contig_ofp = tempfile.TemporaryFile()
        ofps[(chrm, strand)] = contig_ofp
        tss_exons = set( map(tuple, all_tss_exons[(chrm, strand)].tolist()))
        internal_exons = set(map(tuple, all_internal_exons[(chrm, strand)].tolist()))
        tes_exons = set(map(tuple, all_tes_exons[(chrm, strand)].tolist()))
        se_transcripts = set(map(tuple, all_se_transcripts[(chrm, strand)].tolist()))
        jns = all_jns[(chrm, strand)]
        p = multiprocessing.Process( 
            target=build_genes_in_contig, 
            args=[ (chrm, strand), (gene_id, gene_id_lock), contig_ofp,
                   tss_exons, internal_exons, tes_exons, 
                   se_transcripts, jns ] )
        p.start()
        ps[finished_i] = p
    
    for p in ps:
        if p != None:
            p.join()
    
    for ( chrm, strand), contig_ofp in sorted( ofps.iteritems() ):
        if VERBOSE: print >> sys.stderr, "Writing out ", chrm, strand, contig_ofp
        contig_ofp.seek(0)
        ofp.write( contig_ofp.read() )
        contig_ofp.close()
    
    ofp.flush()
    return
    
def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Determine valid transcripts.')
    parser.add_argument( 'elements', type=file, \
        help='BED file containing elements ( output of find_exons ).')
    
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
    
    assert args.threads == 1 or args.out_fname, "An output file must be specified if this is run multi-threaded."
    out_fp = open( args.out_fname, "w" ) if args.out_fname else sys.stdout
    
    log_fp = multi_safe_file( args.log_fname, 'w' ) \
        if args.log_fname else sys.stderr
    
    global num_threads
    num_threads = args.threads
    
    return args.elements, out_fp, log_fp

def load_elements( fp ):
    all_elements = defaultdict( lambda: defaultdict(set) )
    for line in fp:
        if line.startswith( 'track' ): continue
        chrm, start, stop, element_type, score, strand = line.split()[:6]
        all_elements[ element_type ][ (chrm, strand) ].add( 
            (int(start), int(stop)) )
    
    # convert into array
    all_array_elements = defaultdict( 
        lambda: defaultdict(lambda: numpy.zeros(0)) )
    for element_type, elements in all_elements.iteritems():
        for key, contig_elements in elements.iteritems():
            all_array_elements[element_type][key] \
                = numpy.array( sorted( contig_elements ) )

    return all_array_elements

def main():
    elements_fp, out_fp, log_fp = parse_arguments()
    
    if VERBOSE:
        print >> sys.stderr, 'Parsing input...'
    elements = load_elements( elements_fp )
    elements_fp.close()
        
    # internal_exons, tss_exons, tes_exons, jns
    build_genes( elements['internal_exon'], elements['tss_exon'], 
                 elements['tes_exon'], elements['single_exon_gene'], 
                 elements['intron'], 
                 out_fp)
            
    return

if __name__ == "__main__":
    main()
