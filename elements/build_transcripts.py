MAX_NUM_TRANSCRIPTS = 50000
VERBOSE = False
MIN_VERBOSE = False

import sys
import os
import numpy

from collections import namedtuple
CategorizedExons = namedtuple( "CategorizedExons", ["TSS", "TES", "internal"] )

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "./exons/" ) )
from build_genelets import cluster_exons, cluster_overlapping_exons

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../sparsify/" ) )
import transcripts as transcripts_module

sys.path.append( os.path.join(os.path.dirname(__file__), "../file_types") )
from exons_file import parse_exons_files
from junctions_file import parse_junctions_files

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

def build_gtf_lines( gene_name, chrm, gene_strand, transcript_name, transcript, exons ):
    lines = []
    for exon_num, ( start, stop ) in \
            enumerate( exons[exon_index] for exon_index in transcript ):
        lines.append( build_gtf_line(gene_name, chrm, gene_strand, \
                                         transcript_name, exon_num, start, stop) )
        
    return lines

def build_genes( exons, tss_exons, tes_exons, jns ):
    def has_exon( exons, exon ):
        """Determine if exon is in sorted numpy array exons
        """
        index = exons[:,0].searchsorted( exon[0] )
        # search through all of the exons that start with the
        # correct coordinate and determine if the other coordinate
        # corresponds also
        while index < len(exons) and exon[0] == exons[index][0]:
            if exon[1] == exons[index][1]:
                return True
            index += 1
        return False
    
    
    cluster_num = 1
    clustered_exons = {}
    
    # create genelets for each chrm, strand combination that exists 
    # in both all exon types and junctions
    keys = set( exons ).intersection( jns ).intersection( \
        tss_exons ).intersection( tes_exons )
    for key in sorted(keys):
        chrm, strand = key
        all_chrm_exons = numpy.vstack( \
            (exons[key], tss_exons[key], tes_exons[key]) )
        
        genelets = cluster_exons( all_chrm_exons, jns[key] )
        
        # find tss and tes exons and add the grouped exons to the 
        # clustered_exons list
        for all_cluster_exons in genelets:
            cluster_tss_exons = []
            cluster_tes_exons = []
            for exon in all_cluster_exons:
                if has_exon( tss_exons[key], exon ):
                    cluster_tss_exons.append(exon)
                if has_exon( tes_exons[key], exon ):
                    cluster_tes_exons.append(exon)
            
            if len( cluster_tss_exons ) == 0 or len( cluster_tes_exons ) == 0:
                continue
            
            cluster_key = ("cluster_{0:d}".format( cluster_num ), chrm, strand)
            cluster_num += 1
            clustered_exons[ cluster_key ] = CategorizedExons( \
                cluster_tss_exons, cluster_tes_exons, all_cluster_exons)

    """
    for key, val in clustered_exons.iteritems():
        print key
        print val.TSS
        print val.TES
        print val.internal
    """
    
    return clustered_exons

def find_transcripts( chrm, strand, junctions, start_exons, stop_exons, exons ):
    """Return transcripts from exons and junctions provided
    """
    def build_transcripts( exons, conn_exons, start_exons, stop_exons ):
        transcripts = []
        for transcript_index, entry in enumerate(transcripts_module.iter_transcripts( \
                exons, conn_exons, start_exons, stop_exons ) ):
            if transcript_index > MAX_NUM_TRANSCRIPTS:
                return 'TOO MANY TRANSCRIPTS'
            else:
                transcripts.append( entry )
        
        if len( transcripts ) == 0: return 'NO TRANSCRIPTS PRODUCED'
                
        return transcripts
    
    min_start, max_stop = min( min(i) for i in exons), max(max(i) for i in exons)
    conn_exons, conn_exon_scores = junctions.iter_connected_exons( \
        chrm, strand, min_start, max_stop, exons, True )
    
    if strand == '+':
        upstream_exons = start_exons
        downstream_exons = stop_exons
    else:
        assert strand == '-'
        upstream_exons = stop_exons
        downstream_exons = start_exons
    
    return build_transcripts( exons, conn_exons, upstream_exons, downstream_exons )

def build_transcripts_gtf_lines( gene_name, chrm, strand, exons, junctions, log_fp ):
    """ Build the gtf lines corresponding with exons and junctions.
    """
    def write_log( gene_name, len_exons, len_trans, error_str="" ):
        return '{0},{1:d},{2:d},chr{3}:{4:d}-{5:d}:{6},"{7}"\n'.format( \
            gene_name, len_exons, len_trans, chrm, \
                min( min(i) for i in exons), \
                max( max(i) for i in exons), \
                strand, error_str )
    
    def get_indices( exons, some_exons ):
        indices = []
        for value in some_exons:
            index = exons.index( value )
            indices.append( index )
            
        return indices
    
    
    # unpack exons and convert start and stop exons to indices of 
    # sorted exons list
    tss_exons, tes_exons, all_exons = exons
    exons = sorted( all_exons )
    tss_exons = get_indices( exons, tss_exons )
    tes_exons = get_indices( exons, tes_exons )
    
    transcripts = find_transcripts( \
        chrm, strand, junctions, tss_exons, tes_exons, exons )
    
    # if too many or no transcripts were produced
    if isinstance( transcripts, str ):
        # log the region and note that too many transcripts were produced
        log_line = write_log(gene_name, len(exons), 0, transcripts)
        log_fp.write( log_line + "\n" )
        return []
    
    # otherwise, write them to GTF
    all_lines = []
    for transcript_index, transcript in enumerate( transcripts ):
        transcript_name = "{0}_{1:d}".format( gene_name, transcript_index )
        lines = build_gtf_lines( gene_name, chrm, strand, \
                                 transcript_name, transcript, exons )
        all_lines.extend( lines )
    
    log_fp.write( write_log(gene_name, len(exons), len( all_lines ) ) )
    
    return all_lines

def enumerate_transcripts_worker( input_queue, output_queue, jns, log_fp ):
    while not input_queue.empty():
        try:
            ( gene_name, chrm, strand ), exons = input_queue.get(block=False)
        except Queue.Empty:
            break
        
        lines = build_transcripts_gtf_lines( \
            gene_name, chrm, strand, exons, jns, log_fp)
        
        output_queue.put( lines )
    
    return

def write_transcripts( genes, jns, log_fp, out_fp, threads ):
    # create queues to store input and output data
    manager = multiprocessing.Manager()
    input_queue = manager.Queue()
    output_queue = manager.Queue()
    
    for data, exons in genes.iteritems():
        input_queue.put( (data, exons) )
    
    args = ( input_queue, output_queue, jns, log_fp )
    # spawn threads to estimate genes expression
    processes = []
    for thread_id in xrange( threads ):
        p = multiprocessing.Process( target=enumerate_transcripts_worker, args=args )
        p.start()
        processes.append( p )
    
    def process_queue():
        while not output_queue.empty():
            try:
                lines = output_queue.get()
            except Queue.Empty:
                break
            
            if len( lines ) > 0:
                out_fp.write( "\n".join( lines ) + "\n" )
    
    # process output queue
    while any( p.is_alive() for p in processes ):
        process_queue()
    
    # process any remaining lines
    process_queue()
    
    return

def build_objs( exons_fps, tss_exons_fps, tes_exons_fps, junctions_fps ):
    exons = parse_exons_files( exons_fps )
    
    tss_exons = parse_exons_files( tss_exons_fps )
    
    tes_exons = parse_exons_files( tes_exons_fps )
    
    jns = parse_junctions_files( junctions_fps )
    
    return exons, tss_exons, tes_exons, jns

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

    parser.add_argument( '--single-exon', default=False, action= 'store_true', 
        help='Only build single exon genes.')    
    parser.add_argument( '--threads', '-t', type=int , default=1, \
                             help='Number of threads spawn for multithreading.')
    parser.add_argument( '--out-fname', '-o', default='',\
                             help='Output file name. (default: stdout)')
    parser.add_argument( '--log-fname', '-l',\
                             help='Output log file name. (default: sterr)')
    parser.add_argument( '--verbose', '-v', default=False, action= 'store_true', \
                             help='Whether or not to print status information.')
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    transcripts_module.VERBOSE = VERBOSE
    
    out_fp = open( args.out_fname, "w" ) if args.out_fname else sys.stdout
    
    log_fp = multi_safe_file( args.log_fname, 'w' ) if args.log_fname else sys.stderr
    
    return args.exons, args.tss, args.tes, args.junctions, \
        args.single_exon, args.threads, out_fp, log_fp


class DistalExonBndries( set ):
    def __init__( self, strand, is_tss, distal_exons ):
        self.is_tss = is_tss
        self.strand = strand
        
        if is_tss: 
            index = 1
        else: 
            index = 0
        if strand == '+':
            index = 1 - index
        
        self.update( i[1-index] for i in distal_exons )
    
    def _is_distal_exon( self, exon ):
        """Check to see if 'exon' is a distal exon.

        """
        if (self.strand == '+' and self.is_tss) \
                or (self.strand == '-' and not self.is_tss):
            return exon[1] in self
        else:
            return exon[0] in self
    
    def is_tss_exon( self, exon ):
        assert self.is_tss
        return self._is_distal_exon( exon )

    def is_tes_exon( self, exon ):
        assert not self.is_tss
        return self._is_distal_exon( exon )

def iter_single_exon_gene_gtf_lines( exons, jns, tss_exons, tes_exons ):
    se_genes = []
    # search for single exon genes
    for ( chrm, strand ), unclustered_inner_exons in exons.iteritems():
        jn_bndries = set( jns[( chrm, strand )].ravel() )
        tss_bndries = DistalExonBndries( strand, True, tss_exons[( chrm, strand )] )
        tes_bndries = DistalExonBndries( strand, False, tes_exons[( chrm, strand )] )
        exon_clusters = cluster_overlapping_exons( unclustered_inner_exons )
        
        se_gene_exons = defaultdict( list )
        for cluster_id, inner_exons in exon_clusters.iteritems():
            if len( inner_exons ) != 1: continue
            exon = inner_exons[0]
            if exon[0]-1 not in jn_bndries and exon[1]+1 not in jn_bndries:
                #if tss_bndries.is_tss_exon( exon ): # \
                #    #    and tes_bndries.is_tes_exon( exon ):
                se_gene_exons[ cluster_id ].append( (max(1,exon[0]-400), exon[1]+400) )
            
            # test to see if this is a TSS exon
            if len( inner_exons ) == 1 and 1 == len(se_gene_exons[ cluster_id ]):
                start = se_gene_exons[ cluster_id ][0][0]
                stop = se_gene_exons[ cluster_id ][0][1]
                se_genes.append( ( chrm, strand, start, stop ) )
    
    for index, (chrm, strand, start, stop) in enumerate(se_genes):
        gene_id = "SE_gene_%i" % (index + 1)
        trans_id = gene_id
        yield build_gtf_line( gene_id, chrm, strand, trans_id, 1, start, stop )

    return

def main():
    exon_fps, tss_exon_fps, tes_exon_fps, junction_fps, \
        single_exon, threads, out_fp, log_fp = parse_arguments()
    
    if VERBOSE:
        print >> sys.stderr, 'Parsing input...'
    exons, tss_exons, tes_exons, jns = build_objs( \
        exon_fps, tss_exon_fps, tes_exon_fps, junction_fps )
    
    if VERBOSE:
        print >> sys.stderr, 'Clustering exons...'
    genes = build_genes( exons, tss_exons, tes_exons, jns )
    
    if VERBOSE:
        print >> sys.stderr, 'Outputting all gene models...'
    write_transcripts( genes, jns, log_fp, out_fp, threads )
    

    se_gene_gtf_lines_iter = iter_single_exon_gene_gtf_lines( \
        exons, tss_exons, tes_exons, jns )
    out_fp.write( "\n".join(se_gene_gtf_lines_iter) + "\n" )
        
    return

if __name__ == "__main__":
    main()
