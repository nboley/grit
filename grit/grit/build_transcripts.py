import sys, os
import shutil

import time
import traceback

import numpy

from pysam import Fastafile

from itertools import chain
from collections import namedtuple

import threading
import multiprocessing
from multiprocessing.sharedctypes import RawArray, RawValue
from lib.multiprocessing_utils import Pool

import networkx as nx

from lib.multiprocessing_utils import ThreadSafeFile
from transcript import Transcript, Gene
from files.reads import fix_chrm_name_for_ucsc
from proteomics.ORF import find_cds_for_gene
from elements import \
    load_elements, cluster_elements, find_jn_connected_exons

import config

import Queue

import random

GeneElements = namedtuple('GeneElements', 
                          ['id', 'chrm', 'strand',
                           'tss_exons', 'internal_exons', 'tes_exons',
                           'se_transcripts', 'promoter', 'polyas', 
                           'introns'])

SAMPLE_TYPE = None
REP_ID = None

class TooManyCandidateTranscriptsError(Exception):
    pass

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
                    raise TooManyCandidateTranscriptsError, "Too many candidate transcripts"
                
    return transcripts + [ [x,] for x in se_transcripts ]

class MaxIterError( ValueError ):
    pass

def write_gene_to_gtf( ofp, gene ):
    lines = []
    for index, transcript in enumerate(gene.transcripts):
        meta_data = {}
        
        if config.FIX_CHRM_NAMES_FOR_UCSC:
            transcript.chrm = fix_chrm_name_for_ucsc(transcript.chrm)
        assert transcript.gene_id != None
        lines.append( transcript.build_gtf_lines(
                meta_data, source="grit") + "\n" )
    
    ofp.write( "".join(lines) )
    
    return

def write_gene_to_tracking_file( ofp, gene):
    lines = []
    contig_name = gene.chrm
    if config.FIX_CHRM_NAMES_FOR_UCSC:
        contig_name = fix_chrm_name_for_ucsc(contig_name)
    
    for t in gene.transcripts:
        if t.gene_name != None:
            gene_short_name = t.gene_name
        elif t.ref_gene != None:
            gene_short_name = t.ref_gene
        else:
            gene_short_name = '-'
        
        line = [
            # tracking ID
            (t.id).ljust(20), 
            # class code
            ('-' if t.ref_match_class_code == None 
             else t.ref_match_class_code).ljust(10), 
            # nearest ref id
            ('-' if t.ref_trans == None else t.ref_trans).ljust(20),
            # gene unique id
            (t.gene_id).ljust(20), 
            # gene short name
            ('-' if t.gene_name == None else t.gene_name).ljust(20),
            # TSS ID
            ('-').ljust(10), 
            ("%s:%s:%i-%i"%(contig_name, t.strand, t.start, t.stop)).ljust(30),
             # transcript length
            str(t.calc_length()) ]
            
        lines.append("\t".join(line))
    
    ofp.write( "\n".join(lines) + "\n" )
    return

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

def rename_transcripts(gene, ref_genes):
    # find the ref genes that overlap gene
    ref_genes = list(ref_genes.iter_overlapping_genes(
            gene.chrm, gene.strand, gene.start, gene.stop))
    if len(ref_genes) == 0:
        return gene
    
    max_tes_offset = ( 
        config.MAX_DISTAL_SIZE_FOR_MATCH_OFFSET+config.TES_EXON_MERGE_DISTANCE)
    max_tss_offset = ( 
        config.MAX_DISTAL_SIZE_FOR_MATCH_OFFSET+config.TSS_EXON_MERGE_DISTANCE)
    
    for t in gene.transcripts:
        best_match = None
        best_match_score = (0, -1e9, 0)
        introns = set(t.introns)
        for ref_gene in ref_genes:
            for ref_t in ref_gene.transcripts:
                score = (
                    len(introns) - len(set(introns)-set(ref_t.introns)),
                    -(abs(ref_t.start-t.start)+abs(ref_t.stop-t.stop)),
                    0 )
                if score > best_match_score:
                    best_match_score = score
                    best_match = ref_t

        if best_match == None: continue
        t.ref_gene = best_match.gene_id
        t.ref_trans = best_match.id
        t.ref_match_class_code = None
        
        t.gene_name = ( t.ref_gene if best_match.gene_name == None 
                        else best_match.gene_name )
        if len(introns)  == len(best_match.introns) == best_match_score[0] and \
                best_match_score[1] > -400:
            t.ref_match_class_code = '='
        elif len(introns) == len(best_match.introns) == best_match_score[0]:
            t.ref_match_class_code = '='
        elif ( len(introns) == best_match_score[0] and 
               len(best_match.introns) >  best_match_score[0] ):
            t.ref_match_class_code = 'c'
        elif best_match_score[0] > 0:
            t.ref_match_class_code = 'j'
        else:
            t.ref_match_class_code = 'o'
        #print t.id, ref_t.id, best_match_score, len(introns)

    gene_names = set(t.ref_gene for t in gene.transcripts)
    gene.name = "\\".join(gene_names)
    return gene

def build_gene(elements, fasta=None, ref_genes=None):
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

    if len(transcripts) == 0:
        return None
    
    gene = Gene(elements.id, elements.id,
                elements.chrm, elements.strand, 
                gene_min, gene_max, transcripts)

    if fasta != None:
        gene.transcripts = find_cds_for_gene( 
            gene, fasta, only_longest_orf=True )
    
    if ref_genes != None:
        gene = rename_transcripts(gene, ref_genes)
    
    return gene

def build_and_write_gene(gene_elements, output,
                         gtf_ofp, tracking_ofp,
                         fasta, ref_genes ):
    # build the gene with transcripts, and optionally call orfs
    try:
        config.log_statement(
            "Building transcript and ORFs for Gene %s" % gene_elements.id)
        
        gene = build_gene(gene_elements, fasta, ref_genes)
        if gene == None: 
            return
        config.log_statement(
            "FINISHED Building transcript and ORFs for Gene %s" % gene.id)

        # dump a pickle of the gene to a temp file, and set that in the 
        # output manager
        ofname = gene.write_to_file(
            config.get_gene_tmp_fname(gene.id, SAMPLE_TYPE, REP_ID))
                
        output.put((gene.id, len(gene.transcripts), ofname))
        write_gene_to_gtf(gtf_ofp, gene)
        write_gene_to_tracking_file(tracking_ofp, gene)
    except TooManyCandidateTranscriptsError:
        start = min(x[0] for x in chain(
                gene_elements.promoter, gene_elements.polyas))
        stop = max(x[1] for x in chain(
                gene_elements.promoter, gene_elements.polyas))
        config.log_statement(
            "Too many candidate transcripts in %s(%s:%s:%i-%i)" % (
                gene_elements.id, gene_elements.chrm, gene_elements.strand, 
                start, stop), 
            log=True)
        return
    except Exception, inst:
        start = min(x[0] for x in chain(
                gene_elements.promoter, gene_elements.polyas))
        stop = max(x[1] for x in chain(
                gene_elements.promoter, gene_elements.polyas))
        config.log_statement(
            "ERROR building transcript in %s(%s:%s:%i-%i): %s" % (
                gene_elements.id, gene_elements.chrm, gene_elements.strand, 
                start, stop, inst), 
            log=True)
        if config.DEBUG_VERBOSE:
            config.log_statement( traceback.format_exc(), log=True )
    
    return

def build_transcripts_worker( elements, 
                              output,
                              gtf_ofp, tracking_ofp,
                              fasta_fp, ref_genes ):
    # if appropriate, open the fasta file
    if fasta_fp != None: fasta = Fastafile(fasta_fp.name)
    else: fasta = None
    while True:
        config.log_statement("Waiting for gene to process. (%i)" % elements.qsize())
        gene_elements = elements.get()
        if gene_elements == 'FINISHED':
            config.log_statement("")
            return
        build_and_write_gene( gene_elements, output, 
                              gtf_ofp, tracking_ofp,
                              fasta, ref_genes)
    return

def group_elements_in_gene(grpd_exons):    
    for g_start, g_stop in grpd_exons['gene']:
        args = []
        for key in ('tss_exon', 'internal_exon', 'tes_exon', 
                    'single_exon_gene', 'promoter', 'polya', 'intron'):
            if key not in grpd_exons: 
                args.append(set())
            else:
                exons = [tuple(x) for x in grpd_exons[key].tolist()
                         if x[0] >= g_start and x[1] <= g_stop]
                args.append(set(exons))
        yield args

def add_elements_for_contig_and_strand((contig, strand), 
                                       grpd_exons, elements, gene_id_cntr,
                                       output, gtf_ofp, tracking_ofp, 
                                       fasta_fp, ref_genes):
    if fasta_fp != None: fasta = Fastafile(fasta_fp.name)
    else: fasta = None
    
    config.log_statement( 
        "Clustering elements into genes for %s:%s" % ( contig, strand ) )

    """ old code that actually clustered elements
    args = []
    for key in ('tss_exon', 'internal_exon', 'tes_exon', 
                'single_exon_gene', 'promoter', 'polya', 'intron'):
        if key not in grpd_exons: 
            args.append(set())
        else:
            args.append(
                set(map(tuple, grpd_exons[key].tolist())))
    args.append(strand)
    """
    for ( tss_es, internal_es, tes_es,
          se_ts, promoters, polyas, jns ) in group_elements_in_gene(grpd_exons):
        # skip genes without all of the element types
        if len(se_ts) == 0 and (
                len(tes_es) == 0 
                or len( tss_es ) == 0 ):
            continue
        
        with gene_id_cntr.get_lock():
            gene_id = "XLOC_%i" % gene_id_cntr.value
            gene_id_cntr.value += 1
        
        gene_data = GeneElements( gene_id, contig, strand,
                                  tss_es, internal_es, tes_es,
                                  se_ts, promoters, polyas, 
                                  jns )

        try: 
            elements.put(gene_data, timeout=0.1)
        except Queue.Full: 
            build_and_write_gene( gene_data, output, 
                                  gtf_ofp, tracking_ofp,
                                  fasta, ref_genes)
            config.log_statement( 
                "Clustering elements into genes for %s:%s" % ( 
                    contig, strand ) )
    
    config.log_statement( 
        "FINISHED Clustering elements into genes for %s:%s" % (contig, strand))
    return    

def add_elements_for_contig_and_strand_worker(
        args_queue, elements, gene_id_cntr,
        output, gtf_ofp, tracking_ofp, 
        fasta_fp, ref_genes):
    while True:
        args = args_queue.get()
        if args == 'FINISHED': 
            config.log_statement("")
            return
        (contig, strand), grpd_exons = args
        add_elements_for_contig_and_strand(
            (contig, strand), grpd_exons,
            elements, gene_id_cntr,
            output, gtf_ofp, tracking_ofp, 
            fasta_fp, ref_genes)

def feed_elements(raw_elements, elements, 
                  output, gtf_ofp, tracking_ofp, 
                  fasta_fp, ref_genes ):
    all_args = multiprocessing.Queue()
    for (contig, strand), grpd_exons in raw_elements.iteritems():
        all_args.put([(contig, strand), dict(grpd_exons)])
    for i in xrange(config.NTHREADS):
        all_args.put('FINISHED')

    num_add_element_threads = min(len(raw_elements), config.NTHREADS)
    gene_id_cntr = multiprocessing.Value('i', 0)
    nthreads_remaining = multiprocessing.Value('i', num_add_element_threads)
    worker_args = [ all_args, elements, gene_id_cntr,
                    output, gtf_ofp, tracking_ofp, 
                    fasta_fp, ref_genes ]
    cluster_pids = []
    for i in xrange(num_add_element_threads):
        pid = os.fork()
        if pid == 0:
            add_elements_for_contig_and_strand_worker(*worker_args)
            with nthreads_remaining.get_lock():
                nthreads_remaining.value -= 1
                config.log_statement("Finished adding elements (%i left)" 
                                     % nthreads_remaining.value)
            build_transcripts_worker( elements, 
                                      output,
                                      gtf_ofp, tracking_ofp,
                                      fasta_fp, ref_genes )      
            os._exit(0)

        cluster_pids.append(pid)

    while True:
        with nthreads_remaining.get_lock():
            if nthreads_remaining.value == 0:
                for i in xrange(config.NTHREADS+1):
                    elements.put('FINISHED')
                break
        time.sleep(1.0)

    for pid in cluster_pids:
        os.waitpid(pid, 0) 
    
    config.log_statement("Finished adding elements")
    return

def build_transcripts(exons_bed_fp, gtf_ofname, tracking_ofname, 
                      fasta_fp=None, ref_genes=None):
    """Build transcripts
    """    
    # make sure that we're starting from the start of the 
    # elements files
    config.log_statement( "Loading %s" % exons_bed_fp.name, log=True )
    exons_bed_fp.seek(0)
    raw_elements = load_elements( exons_bed_fp )
    config.log_statement( "Finished Loading %s" % exons_bed_fp.name )
    
    gtf_ofp = ThreadSafeFile(gtf_ofname + ".unfinished", "w")
    gtf_ofp.write("track name=%s useScore=1\n" 
                  % ".".join(gtf_ofname.split(".")[:-1]))
    
    tracking_ofp = ThreadSafeFile(tracking_ofname + ".unfinished", "w")
    tracking_ofp.write("\t".join(
            ["tracking_id".ljust(20), 
             "class_code", 
             "nearest_ref_id".ljust(20), 
             "gene_id".ljust(20), 
             "gene_short_name".ljust(20), 
             "tss_id".ljust(10), 
             "locus".ljust(30), 
             "length"]) + "\n")
    
    config.log_statement( "Building Transcripts", log=True )
    manager = multiprocessing.Manager()
    elements = manager.Queue(2*config.NTHREADS)
    output = manager.Queue()

    transcript_building_children_args = [
        elements, output, 
        gtf_ofp, tracking_ofp,
        fasta_fp, ref_genes]

    
    pids = []
    for i in xrange(max(1,config.NTHREADS - len(raw_elements))):
        pid = os.fork()
        if pid == 0:
            build_transcripts_worker(elements, 
                                     output,
                                     gtf_ofp, tracking_ofp,
                                     fasta_fp, ref_genes)
            os._exit(0)
        pids.append(pid)

    elements_feeder_pid = os.fork()
    if elements_feeder_pid == 0:
        feed_elements( raw_elements, elements, 
                       output, gtf_ofp, tracking_ofp, 
                       fasta_fp, ref_genes )
        os._exit(0)

    for pid in pids:
        os.waitpid(pid, 0) 

    os.waitpid(elements_feeder_pid, 0)
    
    genes = []
    while output.qsize() > 0:
        try: 
            genes.append(output.get_nowait())
        except Queue.Empty: 
            continue
    
    assert len(genes) == len(set(genes))
    config.log_statement("Finished building transcripts")

    gtf_ofp.close()
    tracking_ofp.close()

    # we store to unfinished so we know if it errors out early
    shutil.move(gtf_ofname + ".unfinished", gtf_ofname)
    shutil.move(tracking_ofname + ".unfinished", tracking_ofname)

    manager.shutdown()
    
    return genes

def main():
    assert False
    with open(sys.argv[1]) as fp:
        build_transcripts(fp, "tmp", sys.argv[2])

if __name__ == '__main__':
    main()
