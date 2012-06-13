import sys, os

from collections import defaultdict
from itertools import izip
from random import shuffle

import multiprocessing
import time
import Queue

sys.path.append( os.path.join(os.path.dirname(__file__), "..", "file_types") )
from gtf_file import parse_gtf_line

DO_PROFILE = False

def write_unrepresented_trans( missed_known_trans, unrepresented_trans_fp ):
    for (chrm, strand), key_known_trans in missed_known_trans.iteritems():
        for known_trans in key_known_trans:
            line_base = '\t'.join( ( 
                    'chr'+chrm, known_trans[1][2], 'exon', '{0:d}', 
                    '{1:d}', '.', strand, '.', 
                    'gene_id "{0}"; transcript_id "{1}";'.format( 
                        known_trans[1][0], known_trans[1][1] ) ) )
            
            for start, stop in izip( known_trans[0][0::2], known_trans[0][1::2] ):
                unrepresented_trans_fp.write( 
                    line_base.format(start, stop) + '\n' )
    
    return

def make_gtf_lines( orf, chrm, strand, source ):
    # unpack orf data
    (fp_utr, cds, tp_utr, gene_id, trans_id, trans), sources = orf
    
    orf_lines = []
    orf_base_line = '\t'.join( (
        'chr' + chrm, source, '{0}\t{1:d}\t{2:d}\t.', strand, '{3}', 
        'gene_id "{0}"; transcript_id "{1}";'.format( gene_id, trans_id ) ) )
    
    # treat minus strand different in order to calculate frame correctly
    frame = 0
    if strand == '+':
        for start, stop in fp_utr:
            orf_lines.append( orf_base_line.format( '5UTR', start, stop, '.' ) )
        for start, stop in cds:
            orf_lines.append( 
                orf_base_line.format( 'CDS', start, stop, str(frame) ) )
            frame = ((( -(( stop-start+1 )%3 ))%3 ) + frame )%3
        for start, stop in tp_utr:
            orf_lines.append( orf_base_line.format( '3UTR', start, stop, '.' ) ) 
    else:
        assert strand == '-', 'Invalid Strand: {0}'.format( strand )
        for start, stop in reversed( tp_utr ):
            orf_lines.append( orf_base_line.format( '3UTR', start, stop, '.' ) )
        for start, stop in reversed( cds ):
            orf_lines.append( 
                orf_base_line.format( 'CDS', start, stop, str(frame) ) )
            frame = ((( -(( stop-start+1 )%3 ))%3 ) + frame )%3
        for start, stop in reversed( fp_utr ):
            orf_lines.append( orf_base_line.format( '5UTR', start, stop, '.' ) )
        orf_lines.reverse()
    
    trans_trans_id = '_'.join( trans_id.split( '_' )[:-1] )
    gtf_base_line = '\t'.join( (
        'chr' + chrm, source, '{0}\t{1:d}\t{2:d}\t.', strand, '.', 
        'gene_id "{0}"; transcript_id "{1}";'.format( gene_id, trans_trans_id ) ) )
    
    gtf_lines = []
    for start, stop in izip( trans[::2], trans[1::2] ):
        gtf_lines.append( gtf_base_line.format( 'exon', start, stop ) )
   
    return ( '\n'.join( orf_lines ), '\n'.join( gtf_lines ) )

def write_novel_orfs( disc_novel_orfs, novel_orfs_fp, novel_trans_fp ):
    for (chrm, strand), key_novel_orfs in disc_novel_orfs.iteritems():
        for novel_orf in key_novel_orfs:
            orf_lines, trans_lines = make_gtf_lines( 
                novel_orf, chrm, strand, 'novel_orf' )
            novel_orfs_fp.write( orf_lines + '\n' )
            novel_trans_fp.write( trans_lines + '\n' )
    
    novel_orfs_fp.close()
    novel_trans_fp.close()
    
    return

def write_matched_orfs( disc_matched_orfs, matched_orfs_fp, matched_trans_fp ):
    for (chrm, strand), key_matched_orfs in disc_matched_orfs.iteritems():
        for matched_orf in key_matched_orfs:
            orf_lines, trans_lines = make_gtf_lines( 
                matched_orf, chrm, strand, 'known_trans' )
            matched_orfs_fp.write( orf_lines + '\n' )
            matched_trans_fp.write( trans_lines + '\n' )
    
    matched_orfs_fp.close()
    matched_trans_fp.close()
    
    return

def get_sample_specific_novel( disc_orfs_for_novel, matched_fp_utrs, 
                               matched_tp_utrs, matched_sources ):
    novel_orfs = defaultdict( list )
    for (chrm, strand), candidate_novel_orfs in disc_orfs_for_novel.iteritems():
        for cds, utrs in candidate_novel_orfs.iteritems():
            # create initial set of utrs and sources already used for this orf 
            # from matching to known transcripts
            used_fp_utrs = matched_fp_utrs[(chrm, strand)][cds] if cds in \
                matched_fp_utrs[(chrm, strand)] else set()
            used_tp_utrs = matched_tp_utrs[(chrm, strand)][cds] if cds in \
                matched_tp_utrs[(chrm, strand)] else set()
            used_sources = matched_sources[(chrm, strand)][cds] if cds in \
                matched_sources[(chrm, strand)] else set()
            
            # randomize the order of the utrs so that fewer are used *hopefully*
            utrs = utrs.items()
            shuffle( utrs )
            for (fp_utr, tp_utr, gene_id, trans_id, trans), trans_sources in utrs:
                # skip trans where both utrs and all sources for this trans 
                # have already been represented
                if fp_utr in used_fp_utrs and tp_utr in used_tp_utrs and \
                        len( used_sources ) == len( used_sources.union( \
                        trans_sources ) ):
                    continue
                
                used_fp_utrs.add( fp_utr )
                used_tp_utrs.add( tp_utr )
                used_sources.update( trans_sources )
                
                novel_orfs[(chrm, strand)].append(
                    ( ( fp_utr, cds, tp_utr, gene_id, trans_id, trans ), 
                      trans_sources ) )
    
    return novel_orfs

def get_closest_tes_match( external_ann_coords, matched_orfs, strand ):
    """ Get orf from matched_orfs which is closest to the known tes boundary
    """
    known_tes_ext_coord = external_ann_coords[0][1] if strand == '+' \
        else external_ann_coords[0][0]
    
    match_bps_diff = None
    for trans_struct, trans_sources in matched_orfs.iteritems():
        trans_ext_tes_coord = trans_struct[2][-1][-1] if strand == '+' \
            else trans_struct[2][0][0]
        if match_bps_diff == None or \
                abs( trans_ext_tes_coord - known_tes_ext_coord ) < match_bps_diff:
            match_bps_diff = abs( trans_ext_tes_coord - known_tes_ext_coord )
            
            match_trans_struct = trans_struct
            match_trans_sources = trans_sources
    
    return match_trans_struct, match_trans_sources

def get_ann_matched_trans( disc_orfs, known_trans ):
    # lists to store for printing
    disc_matched_orfs = defaultdict( list )
    missed_known_trans = defaultdict( list )
    
    # sets to create novel orfs later
    matched_fp_utrs = defaultdict( lambda: defaultdict( set ) )
    matched_tp_utrs = defaultdict( lambda: defaultdict( set ) )
    matched_sources = defaultdict( lambda: defaultdict( set ) )
    
    for (chrm, strand), key_known_trans in known_trans.iteritems():
        if (chrm, strand) in disc_orfs:
            key_disc_orfs = disc_orfs[ (chrm, strand) ]
            for known_int_coords in key_known_trans.iterkeys():
                if known_int_coords in key_disc_orfs:
                    trans_struct, trans_sources = get_closest_tes_match( 
                        key_known_trans[known_int_coords], 
                        key_disc_orfs[ known_int_coords ], 
                        strand )
                    disc_matched_orfs[ (chrm, strand) ].append( 
                        (trans_struct, trans_sources) )
                    
                    # store items matched indexed by cds for novels later
                    matched_fp_utrs[ (chrm, strand) ][ trans_struct[1] ].add(
                        trans_struct[0] )
                    matched_tp_utrs[ (chrm, strand) ][ trans_struct[1] ].add(
                        trans_struct[2] )
                    
                    matched_sources[ (chrm, strand) ][ 
                        trans_struct[1] ].update( trans_sources )
                else:
                    missed_known_trans[(chrm, strand)].append( 
                        ( ( key_known_trans[ known_int_coords ][0][:1] +
                            known_int_coords +
                            key_known_trans[ known_int_coords ][0][-1:] ),
                          key_known_trans[ known_int_coords ][1] ) )
        else:
            for known_int_coords, known_ext_coords in key_known_trans.iteritems():
                missed_known_trans[(chrm, strand)].append( 
                    ( ( key_known_trans[ known_int_coords ][0][:1] + 
                        known_int_coords +
                        key_known_trans[ known_int_coords ][0][-1:] ),
                      key_known_trans[ known_int_coords ][1] ) )
    
    return disc_matched_orfs, missed_known_trans, \
        matched_fp_utrs, matched_tp_utrs, matched_sources

def parse_known_trans( cdna_fp, ann_fp ):
    def update_known_trans( known_trans, gtf_fp ):
        if gtf_fp == None: return known_trans
        for line in gtf_fp:
            gtfl = parse_gtf_line( line )
            if gtfl == None or gtfl.feature != 'exon': continue
            known_trans[ (gtfl.region.chr, gtfl.region.strand) ][ \
                (gtfl.gene_id, gtfl.trans_id, gtfl.source) ].append( 
                (gtfl.region.start, gtfl.region.stop) )
        
        return known_trans
    
    raw_known_trans = defaultdict( lambda: defaultdict( list ) )
    raw_known_trans = update_known_trans( raw_known_trans, cdna_fp )
    raw_known_trans = update_known_trans( raw_known_trans, ann_fp )
    
    # create a list of internal coordinates from raw_trans
    all_known_trans = defaultdict( dict )
    for key, transcripts in raw_known_trans.iteritems():
        for ids, trans in transcripts.iteritems():
            raw_trans = [ pos for exon in sorted(trans) for pos in exon ]
            # skip single exon transcripts
            if len(raw_trans) <= 2: continue
            
            all_known_trans[ key ][ tuple(raw_trans[1:-1]) ] = \
                ( ( raw_trans[0], raw_trans[-1] ), ids )
    
    return all_known_trans

def parse_discovered_orfs( orfs_fp, output_queue ):
    """Parse gtf file with known orfs
    orfs are stored as an interval searchable tree and a dict indexed by the
    tree search results
    """
    source_name = os.path.basename( orfs_fp.name )
    
    raw_orfs = defaultdict( lambda: defaultdict( lambda: defaultdict(list) ) )
    for line in orfs_fp:
        gtfl = parse_gtf_line( line )
        if gtfl == None or gtfl.feature not in \
                ( 'CDS', '5UTR', '3UTR', 'exon' ): continue
        main_key = (gtfl.gene_id, gtfl.region.chr, gtfl.region.strand)
        # losing frame and exon number info here
        raw_orfs[ main_key ][ gtfl.trans_id ][ gtfl.feature ].append( 
            ( gtfl.region.start, gtfl.region.stop ) )
    
    orfs_data = []
    for ( gene_id, chrm, strand ), transcripts in raw_orfs.iteritems():
        for trans_id, feature_types in transcripts.iteritems():
            # skip non-coding transcripts
            if 'exon' in feature_types: continue
            # skip transcripts without all 3 required regions
            if 'CDS' not in feature_types or '5UTR' not in feature_types or \
                    '3UTR'not in feature_types: continue
            
            fp_utr = tuple(sorted( feature_types[ '5UTR' ] ))
            cds = tuple(sorted( feature_types[ 'CDS' ] ))
            tp_utr = tuple(sorted( feature_types[ '3UTR' ] ))
            
            # get the internal structure of teh transcript from feature types
            upstrm_utr = fp_utr if strand == '+' else tp_utr
            dnstrm_utr = tp_utr if strand == '+' else fp_utr
            
            trans = [ pos for exon in upstrm_utr for pos in exon ][:-1]
            if upstrm_utr[-1][1] != cds[0][0] - 1:
                trans.append( upstrm_utr[-1][1] )
                trans.append( cds[0][0] )
            trans.extend( [pos for exon in cds for pos in exon][1:-1] )
            if dnstrm_utr[0][0] != cds[-1][1] + 1:
                trans.append( cds[-1][1] )
                trans.append( dnstrm_utr[0][0] )
            trans.extend( [pos for exon in dnstrm_utr for pos in exon][1:] )
            trans = tuple( trans )
            
            orfs_data.append( ( (chrm, strand), fp_utr, cds, tp_utr, 
                                gene_id, trans_id, trans ) )
    
    output_queue.put( ( source_name, orfs_data ) )
    
    return

def parse_discovered_orfs_worker( input_queue, output_queue ):
    while not input_queue.empty():
        try:
            orf_fn = input_queue.get(block=False)
        except Queue.Empty:
            break
        
        with open( orf_fn ) as orfs_fp:
            parse_discovered_orfs( orfs_fp, output_queue )
    
    return

def parse_all_discovered_orfs( orfs_fns, threads ):
    # create queues to store input and output data
    manager = multiprocessing.Manager()
    input_queue = manager.Queue()
    output_queue = manager.Queue()
    
    for orfs_fn in orfs_fns:
        input_queue.put( orfs_fn )
    
    args = ( input_queue, output_queue )
    # spawn threads to estimate genes expression
    processes = []
    for thread_id in xrange( threads ):
        p = multiprocessing.Process( target=parse_discovered_orfs_worker, 
                                     args=args )
        p.start()
        processes.append( p )
    
    disc_orfs = defaultdict( lambda: defaultdict( lambda: defaultdict(list) ) )
    disc_orfs_for_novel = defaultdict( 
        lambda: defaultdict( lambda: defaultdict(list) ) )
    
    def process_queue():
        while not output_queue.empty():
            try:
                source_name, orfs_data = \
                    output_queue.get(block=False)
            except Queue.Empty:
                break
            
            for chrm_strand, fp_utr, cds, tp_utr, gene_id, trans_id, trans \
                    in orfs_data:
                disc_orfs[ chrm_strand ][ trans[1:-1] ][ \
                    (fp_utr, cds, tp_utr, gene_id, trans_id, trans) ].append(\
                    source_name )
                
                disc_orfs_for_novel[ chrm_strand ][ cds ][ \
                    (fp_utr, tp_utr, gene_id, trans_id, trans) ].append( \
                    source_name )
        
        return
    
    
    while any( p.is_alive() for p in processes ):
        process_queue()
        time.sleep(1)
    
    process_queue()
    
    return disc_orfs, disc_orfs_for_novel

def build_objects( orfs_fns, cdna_fp, ann_fp, out_prefix, threads ):
    matched_orfs_fp = open( out_prefix + '.known_orfs.gtf', 'w' )
    matched_trans_fp = open( out_prefix + '.known_trans.gtf', 'w' )
    novel_orfs_fp = open( out_prefix + '.novel_orfs.gtf', 'w' )
    novel_trans_fp = open( out_prefix + '.novel_trans.gtf', 'w' )
    unrepresented_trans_fp = open( out_prefix + '.unrepresented_trans.gtf', 'w' )
    
    if VERBOSE: print >> sys.stderr, 'parsing discovered orfs'
    disc_orfs, disc_orfs_for_novel = \
        parse_all_discovered_orfs( orfs_fns, threads )
    
    if VERBOSE: print >> sys.stderr, 'parsing known trans'
    known_trans = parse_known_trans( cdna_fp, ann_fp )
    if cdna_fp != None: cdna_fp.close()
    if ann_fp != None: ann_fp.close()
    
    return disc_orfs, disc_orfs_for_novel, known_trans, \
        matched_orfs_fp, matched_trans_fp, novel_orfs_fp, novel_trans_fp, \
        unrepresented_trans_fp

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(
        description = 'Find the transcripts that match a cDNA (or annotation) '
        'transcript and then create the 5\' and 3\' UTR combinations sample '
        'specific.' )
    parser.add_argument(
        'orfs_gtfs', nargs='*',
        help='GTF files with CDS and UTR regions.' )
    parser.add_argument(
        '--cDNA-gtf', type=file,
        help="GTF file of cDNAs. Used to join known 5' and 3' UTRs.")
    parser.add_argument(
        '--annotation', '-a', type=file,
        help='GTF file of annotated genes and transcripts. '
        "Used to join known 5' and 3' UTRs.")
    
    parser.add_argument(
        '--threads', type=int, default=1,
        help='Number of threads.')
    
    parser.add_argument(
        '--out_prefix', '-o', required=True,
        help='Prefix of output file.')
    parser.add_argument( 
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    if args.cDNA_gtf == None and args.annotation == None:
        print parser.print_help()
        sys.exit()
    
    # set flag args
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.orfs_gtfs, args.cDNA_gtf, args.annotation, args.out_prefix, \
        args.threads

def main():
    disc_orfs_fns, cdna_fp, ann_fp, out_prefix, threads = parse_arguments()
    disc_orfs, disc_orfs_for_novel, known_trans, matched_orfs_fp, \
        matched_trans_fp, novel_orfs_fp, novel_trans_fp, \
        unrepresented_trans_fp = \
        build_objects( disc_orfs_fns, cdna_fp, ann_fp, out_prefix, threads )
    
    if VERBOSE: print >> sys.stderr, 'finding orfs that match known trans'
    disc_matched_orfs, missed_known_trans, matched_fp_utrs, matched_tp_utrs, \
        matched_sources = get_ann_matched_trans( disc_orfs, known_trans )
    
    if VERBOSE: print >> sys.stderr, 'finding novel orfs'
    disc_novel_orfs =  get_sample_specific_novel( 
        disc_orfs_for_novel, matched_fp_utrs, matched_tp_utrs, matched_sources )
    
    if VERBOSE: print >> sys.stderr, 'writing output'
    write_matched_orfs( disc_matched_orfs, matched_orfs_fp, matched_trans_fp )
    write_novel_orfs( disc_novel_orfs, novel_orfs_fp, novel_trans_fp )
    write_unrepresented_trans( missed_known_trans, unrepresented_trans_fp )
    
    return

if __name__ == '__main__':
    if DO_PROFILE:
        import cProfile
        cProfile.run( 'main()' )
    else:
        main()
