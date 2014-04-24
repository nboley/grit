import sys, os

import time
import traceback

import numpy

from pysam import Fastafile

from itertools import izip, chain
from collections import defaultdict, namedtuple

import multiprocessing
from multiprocessing.sharedctypes import RawArray, RawValue
from lib.multiprocessing_utils import Pool

import networkx as nx

from grit.transcript import Transcript, Gene
from grit.files.reads import fix_chrm_name_for_ucsc
from grit.proteomics.ORF import find_cds_for_gene

import config

import cPickle as pickle
import tempfile


GeneElements = namedtuple('GeneElements', 
                          ['id', 'chrm', 'strand',
                           'tss_exons', 'internal_exons', 'tes_exons',
                           'se_transcripts', 'promoter', 'polyas', 
                           'introns'])

def find_overlapping_exons(exons):
    overlapping_exons_mapping = set()
    for o_i, (o_start, o_stop) in enumerate(exons):
        for i_i, (i_start, i_stop) in enumerate(exons):
            if i_i > o_i: break
            if not (i_stop < o_start or i_start > o_stop):
                overlapping_exons_mapping.add( (min(i_i, o_i), max(i_i, o_i)) )
    
    return list(overlapping_exons_mapping)


def find_jn_connected_exons(exons, jns, strand):
    edges = set()
    
    # build mappings from exon starts to indices
    exon_starts_map = defaultdict(list)
    exon_stops_map = defaultdict(list)
    for start, stop in exons:
        exon_starts_map[start].append( (start, stop) )
        exon_stops_map[stop].append( (start, stop ) )
    
    for jn in jns:
        for start in exon_stops_map[jn[0]-1]:
            for stop in exon_starts_map[jn[1]+1]:
                if strand == '+':
                    edges.add((tuple(jn), start, stop))
                else:
                    edges.add((tuple(jn), stop, start))
    
    return edges

def iter_nonoverlapping_exons(exons):
    if len(exons) == 0: return
    G = nx.Graph()
    G.add_nodes_from(xrange(len(exons)))
    overlapping_exons = find_overlapping_exons(exons)
    G.add_edges_from( overlapping_exons )
    
    for gene in nx.connected_components(G):
        clustered_exons = [ exons[exon_i] for exon_i in gene ]
        if len( clustered_exons ) == 1:
            yield clustered_exons[0]
    
    return

def cluster_elements( tss_exons, internal_exons, tes_exons, se_transcripts, 
                   promoters, polyas, jns, strand ):
    assert isinstance( tss_exons, set )
    assert isinstance( internal_exons, set )
    assert isinstance( tes_exons, set )
    assert isinstance( se_transcripts, set )
    assert isinstance( promoters, set )
    assert isinstance( polyas, set )
    
    all_exons = sorted( chain(tss_exons, internal_exons, 
                              tes_exons, se_transcripts,
                              promoters, polyas) )
    if len(all_exons) == 0: return

    G = nx.Graph()
    G.add_nodes_from(xrange(len(all_exons)))
    overlapping_exons = find_overlapping_exons(all_exons)
    G.add_edges_from( overlapping_exons )
    jns_and_connected_exons = find_jn_connected_exons(all_exons, jns, strand)
    G.add_edges_from(
        (start, stop) for jn, start, stop in jns_and_connected_exons)
    edge_jn_map = dict( ((start, stop), jn)
                        for jn, start, stop in jns_and_connected_exons )
    for gene in nx.connected_component_subgraphs(G):
        exons = gene.nodes()
        jns = set(edge_jn_map[edge] for edge in gene.edges() 
                  if edge in edge_jn_map)
        
        yield ( tss_exons.intersection( exons ),
                tes_exons.intersection( exons ),
                internal_exons.intersection( exons ),
                se_transcripts.intersection( exons ),
                promoters.intersection( promoters ),
                polyas.intersection( polyas ),
                sorted(jns) 
              )
    
    return

def build_transcripts_from_elements( 
        tss_exons, internal_exons, tes_exons, se_transcripts, jns, strand ):
                       
    # build a directed graph, with edges leading from exon to exon via junctions
    all_exons = sorted(chain(tss_exons, internal_exons, tes_exons))
    graph = nx.DiGraph()
    graph.add_nodes_from( tss_exons )
    graph.add_nodes_from( internal_exons )
    graph.add_nodes_from( tes_exons )
    
    edges = find_jn_connected_exons(all_exons, jns, strand )
    graph.add_edges_from( (start, stop) for jn, start, stop in edges )
    transcripts = []
    for tss in tss_exons:
        for tes in tes_exons:
            cutoff = config.MAX_NUM_CANDIDATE_TRANSCRIPTS-len(transcripts)+1
            for transcript in nx.all_simple_paths(graph, tss, tes, cutoff):
                transcripts.append( sorted(transcript) )
                if len(transcripts) > config.MAX_NUM_CANDIDATE_TRANSCRIPTS:
                    raise ValueError, "Too many candidate transcripts"
    
    return transcripts + [ [x,] for x in se_transcripts ]

class MaxIterError( ValueError ):
    pass

def extract_elements_from_genes( genes ):
    all_elements = defaultdict( lambda: defaultdict(set) )
    for gene in genes:
        for key, val in gene.extract_elements().iteritems():
            all_elements[(gene.chrm, gene.strand)][key].update(val)

    
    return convert_elements_to_arrays( all_elements )

class ThreadSafeFile( file ):
    def __init__( *args ):
        file.__init__( *args )
        args[0].lock = multiprocessing.Lock()

    def write( self, string ):
        with self.lock:
            file.write( self, string )
            self.flush()

def write_gene_to_gtf( ofp, gene ):
    lines = []
    for index, transcript in enumerate(gene.transcripts):
        meta_data = {}
        
        if config.FIX_CHRM_NAMES_FOR_UCSC:
            transcript.chrm = fix_chrm_name_for_ucsc(transcript.chrm)
        lines.append( transcript.build_gtf_lines(
                gene.id, meta_data, source="grit") + "\n" )
    
    ofp.write( "".join(lines) )
    ofp.flush()
    
    return

def convert_elements_to_arrays(all_elements):
    # convert into array
    all_array_elements = defaultdict( 
        lambda: defaultdict(lambda: numpy.zeros(0)) )
    for key, elements in all_elements.iteritems():
        for element_type, contig_elements in elements.iteritems():
            all_array_elements[key][element_type] \
                = numpy.array( sorted( contig_elements ) )

    return all_array_elements

def load_elements( fp ):
    all_elements = defaultdict( lambda: defaultdict(set) )
    for line in fp:
        if line.startswith( 'track' ): continue
        chrm, start, stop, element_type, score, strand = line.split()[:6]
        # subtract 1 from stop becausee beds are closed open, and we 
        # wnat everything in 0-based closed-closed
        all_elements[(chrm, strand)][element_type].add( 
            (int(start), int(stop)-1) )
    
    return convert_elements_to_arrays(all_elements)

def find_matching_promoter_for_transcript(transcript, promoters):
    # find the promoter that starts at the same basepair
    # If it extends beyond the first exon, we truncate the
    # promoter at the end of the first exon
    tss_exon = transcript.exons[0] if transcript.strand == '+' \
        else transcript.exons[-1] 
    matching_promoter = None
    for promoter in promoters:
        if transcript.strand == '-' and promoter[1] == tss_exon[1]:
            matching_promoter = (max(promoter[0], tss_exon[0]), promoter[1])
        elif transcript.strand == '+' and promoter[0] == tss_exon[0]:
            matching_promoter = (promoter[0], min(promoter[1], tss_exon[1]))
    
    return matching_promoter

def find_matching_polya_region_for_transcript(transcript, polyas):
    # find the polya that ends at the same basepair
    # If it extends beyond the tes exon, we truncate the
    # polya region
    tes_exon = transcript.exons[-1] if transcript.strand == '+' \
        else transcript.exons[0] 
    matching_polya = None
    for polya in polyas:
        if transcript.strand == '+' and polya[1] == tes_exon[1]:
            matching_polya = (max(polya[0], tes_exon[0]), polya[1])
        elif transcript.strand == '-' and polya[0] == tes_exon[0]:
            matching_polya = (polya[0], min(polya[1], tes_exon[1]))
    
    return matching_polya

def build_gene(elements, fasta=None):
    gene_min = min( min(e) for e in chain(
            elements.tss_exons, elements.tes_exons, elements.se_transcripts))
    gene_max = max( max(e) for e in chain(
            elements.tss_exons, elements.tes_exons, elements.se_transcripts))
    
    transcripts = []
    for i, exons in enumerate( build_transcripts_from_elements( 
            elements.tss_exons, elements.internal_exons, elements.tes_exons,
            elements.se_transcripts, elements.introns, elements.strand ) ):
        transcript = Transcript(
            "%s_%i" % ( elements.id, i ), elements.chrm, elements.strand, 
            exons, cds_region=None, gene_id=elements.id)
        transcript.promoter = find_matching_promoter_for_transcript(
            transcript, elements.promoter)
        transcript.polya_region = find_matching_polya_region_for_transcript(
            transcript, elements.polyas)
        transcripts.append( transcript )

    gene = Gene(elements.id, elements.chrm, elements.strand, 
                gene_min, gene_max, transcripts)

    if fasta != None:
        gene.transcripts = find_cds_for_gene( 
            gene, fasta, only_longest_orf=True )
    return gene

def worker( elements, elements_lock, 
            output, output_lock, 
            gtf_ofp, fasta_fp ):
    # if appropriate, open the fasta file
    if fasta_fp != None: fasta = Fastafile(fasta_fp.name)
    else: fasta = None
    
    while True:
        try:
            with elements_lock:
                gene_elements = elements.pop()
        except IndexError:
            return

        config.log_statement(
            "Building transcript and ORFs for Gene %s" % gene_elements.id)
        
        # build the gene with transcripts, and optionally call orfs
        try:
            gene = build_gene(gene_elements, fasta)
            config.log_statement(
                "FINISHED Building transcript and ORFs for Gene %s" % gene.id)

            # dump a pickel of the gene to a temp file, and set that in the 
            # output manager
            ofname = gene.write_to_temp_file()
            with output_lock: 
                output.append((gene.id, len(gene.transcripts), ofname))
            
            write_gene_to_gtf(gtf_ofp, gene)
        except Exception, inst:
            config.log_statement(
                "ERROR building transcript in %s: %s"%(gene_elements.id, inst) )
            raise
    
    return

def add_elements_for_contig_and_strand((contig, strand), grpd_exons,
                                       elements, elements_lock):
    gene_id_num = 1
    config.log_statement( 
        "Clustering elements into genes for %s:%s" % ( contig, strand ) )
    for ( tss_es, tes_es, internal_es, 
          se_ts, promoters, polyas, jns ) in cluster_elements( 
            set(map(tuple, grpd_exons['tss_exon'].tolist())), 
            set(map(tuple, grpd_exons['internal_exon'].tolist())), 
            set(map(tuple, grpd_exons['tes_exon'].tolist())), 
            set(map(tuple, grpd_exons['single_exon_gene'].tolist())),
            set(map(tuple, grpd_exons['promoter'].tolist())), 
            set(map(tuple, grpd_exons['polya'].tolist())), 
            set(map(tuple, grpd_exons['intron'].tolist())), 
            strand):
        # skip genes without all of the element types
        if len(se_ts) == 0 and (
                len(tes_es) == 0 
                or len( tss_es ) == 0 ):
            continue

        gene_id_num += 1
        gene_id = "%s_%s_%i" % ( 
            contig, 'm' if strand == '-' else 'p', gene_id_num )

        gene_data = GeneElements( gene_id, contig, strand,
                                  tss_es, internal_es, tes_es,
                                  se_ts, promoters, polyas, 
                                  grpd_exons['intron'] )
        
        with elements_lock:
            elements.append(gene_data)
    
    config.log_statement("")
    return    

def build_transcripts(exons_bed_fp, ofprefix, fasta_fp=None):
    """Build transcripts
    """
    
    # make sure that we're starting from the start of the 
    # elements files
    config.log_statement( "Loading %s" % exons_bed_fp.name )
    exons_bed_fp.seek(0)
    raw_elements = load_elements( exons_bed_fp )
    config.log_statement( "Finished Loading %s" % exons_bed_fp.name )
    
    gtf_ofp = ThreadSafeFile("%s.gtf" % ofprefix, "w")
    gtf_ofp.write("track name=%s useScore=1\n" % ofprefix)

    manager = multiprocessing.Manager()
    elements = manager.list()
    elements_lock = manager.Lock()    
    
    all_args = []
    for (contig, strand), grpd_exons in raw_elements.iteritems():
        all_args.append([
                (contig, strand), grpd_exons, elements, elements_lock])
    
    if config.NTHREADS in (None, 1):
        for args in all_args:
            add_elements_for_contig_and_strand(*args)
    else:
        p = Pool(config.NTHREADS)
        p.apply(add_elements_for_contig_and_strand, all_args)

    output = manager.list()
    output_lock = manager.Lock()    

    with open("%s.gtf" % ofprefix, "w") as ofp:
        args = [elements, elements_lock, output, output_lock, ofp, fasta_fp]
        if config.NTHREADS in (None, 1):
            worker(*args)
        else:
            ps = []
            for i in xrange(config.NTHREADS):
                p = multiprocessing.Process(target=worker, args=args)
                p.daemon = True
                p.start()
                ps.append(p)
            for p in ps:
                p.join()
    
    genes = sorted(output)
    manager.shutdown()
    config.log_statement("Finished building transcripts")
    
    return genes

def main():
    assert False
    with open(sys.argv[1]) as fp:
        build_transcripts(fp, "tmp", sys.argv[2])

if __name__ == '__main__':
    main()
