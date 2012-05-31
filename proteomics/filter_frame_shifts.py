import sys, os

DO_PROFILE = False

LOW_CDS_JN_RATIO_THRESH = 0.1
HIGH_CDS_JN_RATIO_THRESH = 0.3

ACCEPT_ALL_UTRS = True
UTR_JN_RATIO_THRESH = 0.1

import re
import numpy

from collections import defaultdict
from itertools import izip

sys.path.append( os.path.join(os.path.dirname(__file__), "..", "file_types") )
from gtf_file import parse_gtf_line

sys.path.append( os.path.join( os.path.dirname(__file__), "..", "file_types" ) )
from junctions_file import parse_junctions_file

def write_all_filtered_orfs( out_prefix, ki_fname ):
    out_fname = out_prefix + '.all_filtered.gtf'
    as_filtered_fname = out_prefix + '.alternative_spliced_orfs.gtf'
    
    # get only exact protein matches and same start stop class
    cat_cmd = "cat {0} {1} > {2}".format( as_filtered_fname, ki_fname, out_fname )
    os.system( cat_cmd )
    
    return

def write_novel_orfs( as_orfs, as_orfs_fp, disc_orf_lines ):
    for key in as_orfs:
        as_orfs_fp.write( disc_orf_lines[key] )
    return

def write_shifted_orfs( shifted_orfs, shifted_orfs_fp, disc_orf_lines ):
    for key in shifted_orfs:
        shifted_orfs_fp.write( disc_orf_lines[key] )
    return

def is_valid_mut( disc_intron, known_intron, jn_ratio_thresh, jns ):
    # shift introns to match junction format (closed-closed, interior coords)
    disc_intron = ( disc_intron[0]+1, disc_intron[1]-1 )
    known_intron = ( known_intron[0]+1, known_intron[1]-1 )
    # check if both junctions exist
    if disc_intron in jns and known_intron in jns:
        if jns[disc_intron] / float( jns[known_intron]) >= jn_ratio_thresh:
            return True
    
    return False

def is_valid_utr_don_site_mut( disc_intron, known_don_sites, jns ):
    # if this junction mutation does not change the frame and 
    # passes the low jn ratio threshold
    for known_don_site in known_don_sites:
        known_intron = ( known_don_site, disc_intron[1] )
        if is_valid_mut( disc_intron, known_intron, UTR_JN_RATIO_THRESH, jns ):
            return True
    
    return False

def is_valid_utr_acc_site_mut( disc_intron, known_acc_sites, jns ):
    # if this junction mutation does not change the frame and 
    # passes the low jn ratio threshold
    for known_acc_site in known_acc_sites:
        known_intron = ( disc_intron[0], known_acc_site )
        if is_valid_mut( disc_intron, known_intron, UTR_JN_RATIO_THRESH, jns ):
            return True
    
    return False

def is_valid_cds_don_site_mut( disc_intron, known_don_sites, jns ):
    disc_intron_mod = (disc_intron[1] - disc_intron[0] + 1) % 3
    
    # if this junction mutation does not change the frame and 
    # passes the low jn ratio threshold
    for known_don_site in [site for site in known_don_sites if
                           (disc_intron[1] - site + 1) % 3 == disc_intron_mod]:
        known_intron = ( known_don_site, disc_intron[1] )
        if is_valid_mut(disc_intron, known_intron, LOW_CDS_JN_RATIO_THRESH, jns):
            return True
    
    # if this junction mutation changes the frame and passes the high jn ratio threshold
    for known_don_site in [site for site in known_don_sites if
                           (disc_intron[1] - site + 1) % 3 != disc_intron_mod]:
        known_intron = ( known_don_site, disc_intron[1] )
        if is_valid_mut(disc_intron, known_intron, HIGH_CDS_JN_RATIO_THRESH, jns):
            return True
    
    return False

def is_valid_cds_acc_site_mut( disc_intron, known_acc_sites, jns ):
    disc_intron_mod = (disc_intron[1] - disc_intron[0] + 1) % 3
    
    # if this junction mutation does not change the frame and 
    # passes the low jn ratio threshold
    for known_acc_site in [site for site in known_acc_sites if 
                           (site - disc_intron[0] + 1 ) % 3 == disc_intron_mod]:
        known_intron = ( disc_intron[0], known_acc_site )
        if is_valid_mut( 
            disc_intron, known_intron, LOW_CDS_JN_RATIO_THRESH, jns ):
            return True
    
    # if this junction mutation changes the frame and passes the high jn ratio threshold
    for known_acc_site in [site for site in known_acc_sites if
                           (site - disc_intron[0] + 1 ) % 3 != disc_intron_mod]:
        known_intron = ( disc_intron[0], known_acc_site )
        if is_valid_mut( 
            disc_intron, known_intron, HIGH_CDS_JN_RATIO_THRESH, jns ):
            return True
    
    return False

def is_valid_trans( trans_data, known_introns, jns ):
    disc_orf = sorted( trans_data[ 'CDS' ] )
    # if cds contains a splice
    if len( disc_orf ) > 1:
        disc_cds_introns = izip( zip( *disc_orf )[1][:-1], 
                                 zip( *disc_orf )[0][1:] )
        for disc_cds_intron in disc_cds_introns:
            # if the discovered intron is a known intron
            if disc_cds_intron in known_introns['CDS']['intron']:
                continue
            
            # if there is a possible acceptor site mutation
            if disc_cds_intron[0] in known_introns['CDS']['donor']:
                known_acc_sites = known_introns['CDS']['donor'][ disc_cds_intron[0] ]
                if is_valid_cds_acc_site_mut(disc_cds_intron, known_acc_sites, jns):
                    continue
            
            # if there is a possible donor site mutation
            if disc_cds_intron[1] in known_introns['CDS']['acceptor']:
                known_don_sites = known_introns['CDS']['acceptor'][ disc_cds_intron[1] ]
                if is_valid_cds_don_site_mut(disc_cds_intron, known_don_sites, jns):
                    continue
            
            # if no valid mutation of this intron is found then return False
            return False
    
    if ACCEPT_ALL_UTRS: return True
    
    # filter utr junctions
    disc_fp_utr = sorted( trans_data[ '5UTR' ] )
    # if 5' utr contains a splice
    if len( disc_fp_utr ) > 1:
        disc_fp_utr_introns = izip( zip( *disc_fp_utr )[1][:-1], 
                                    zip( *disc_fp_utr )[0][1:] )
        for disc_fp_utr_intron in disc_fp_utr_introns:
            # if the discovered intron is a known intron
            if disc_fp_utr_intron in known_introns['UTR']['intron']:
                continue
            
            # if there is a possible acceptor site mutation
            if disc_fp_utr_intron[0] in known_introns['UTR']['donor']:
                known_acc_sites = known_introns['UTR']['donor'][ disc_fp_utr_intron[0] ]
                if is_valid_utr_acc_site_mut(disc_fp_utr_intron, known_acc_sites, jns):
                    continue
            
            # if there is a possible donor site mutation
            if disc_fp_utr_intron[1] in known_introns['UTR']['acceptor']:
                known_don_sites = known_introns['UTR']['acceptor'][disc_fp_utr_intron[1]]
                if is_valid_utr_don_site_mut(disc_fp_utr_intron, known_don_sites, jns):
                    continue
            
            # if no valid mutation of this intron is found then return False
            #print 'fails 5": chr2L:{0[0]:d}-{0[1]:d}'.format( disc_fp_utr_intron )
            return False
    
    disc_tp_utr = sorted( trans_data[ '3UTR' ] )
    # if 3' utr contains a splice
    if len( disc_tp_utr ) > 1:
        disc_tp_utr_introns = izip( zip( *disc_tp_utr )[1][:-1], 
                                    zip( *disc_tp_utr )[0][1:] )
        for disc_tp_utr_intron in disc_tp_utr_introns:
            # if the discovered intron is a known intron
            if disc_tp_utr_intron in known_introns['UTR']['intron']:
                continue
            
            # if there is a possible acceptor site mutation
            if disc_tp_utr_intron[0] in known_introns['UTR']['donor']:
                known_acc_sites = known_introns['UTR']['donor'][ disc_tp_utr_intron[0] ]
                if is_valid_utr_acc_site_mut(disc_tp_utr_intron, known_acc_sites, jns):
                    continue
            
            # if there is a possible donor site mutation
            if disc_tp_utr_intron[1] in known_introns['UTR']['acceptor']:
                known_don_sites = known_introns['UTR']['acceptor'][disc_tp_utr_intron[1]]
                if is_valid_utr_don_site_mut(disc_tp_utr_intron, known_don_sites, jns):
                    continue
            
            # if no valid mutation of this intron is found then return False
            #print 'fails 3": chr2L:{0[0]:d}-{0[1]:d}'.format( disc_tp_utr_intron )
            return False
    
    return True

def filter_shifted_trans( disc_orfs, known_introns, all_jns ):
    as_orfs = []
    shifted_orfs = []
    
    for (chrm, strand), key_disc_orfs in disc_orfs.iteritems():
        key_known_introns = known_introns[(chrm, strand)]
        key_jns = all_jns[(chrm, strand)]
        for trans_id, trans_data in key_disc_orfs.iteritems():
            if 'CDS' not in trans_data or '5UTR' not in trans_data or \
                    '3UTR' not in trans_data: continue
            
            if is_valid_trans( trans_data, key_known_introns, key_jns ):
                as_orfs.append( ( chrm, strand, trans_id ) )
            else:
                shifted_orfs.append( ( chrm, strand, trans_id ) )
    
    return as_orfs, shifted_orfs

def parse_jns( jn_fp ):
    all_jns = defaultdict( dict )
    jns = parse_junctions_file(jn_fp)
    if VERBOSE:
        print >> sys.stderr, "Finished parsing", jn_fp.name
    for jn, cnt, grp in jns.iter_jns_and_cnts_and_grps():
        all_jns[ (jn.chr, jn.strand) ][ (jn.start, jn.stop) ] = cnt
    
    return all_jns

def parse_known_orfs( ann_fp ):
    known_orfs = defaultdict( lambda: defaultdict( list ) )
    known_utrs = defaultdict( lambda: defaultdict( list ) )
    for line in ann_fp:
        gtfl = parse_gtf_line( line )
        if gtfl != None and gtfl.feature == 'CDS':
            known_orfs[(gtfl.region.chr, gtfl.region.strand)][ 
                (gtfl.trans_id)].append( (gtfl.region.start, gtfl.region.stop) ) 
        
        elif gtfl != None and gtfl.feature in \
                ('five_prime_UTR', 'three_prime_UTR'):
            known_utrs[(gtfl.region.chr, gtfl.region.strand)][ \
                (gtfl.trans_id, gtfl.feature)].append( \
                (gtfl.region.start, gtfl.region.stop) )
    
    known_introns = {}
    # initialize known introns structure
    for chrm, strand in known_orfs.iterkeys():
        known_introns[(chrm, strand)] = { 'CDS':{'intron':set(), 
                                                 'donor':defaultdict( set ),
                                                 'acceptor':defaultdict( set ) },
                                          'UTR':{'intron':set(), 
                                                 'donor':defaultdict(set),
                                                 'acceptor':defaultdict(set) } }
    
    for (chrm, strand), key_known_orfs in known_orfs.iteritems():
        for regions in key_known_orfs.itervalues():
            regions.sort()
            introns = izip( zip( *regions )[1][:-1], zip( *regions )[0][1:] )
            for intron in introns:
                known_introns[(chrm, strand)]['CDS']['intron'].add( intron )
                known_introns[(chrm, strand)]['CDS']['donor'][ \
                    intron[0] ].add( intron[1] )
                known_introns[(chrm, strand)]['CDS']['acceptor'][ \
                    intron[1] ].add( intron[0] )
    
    for (chrm, strand), key_known_utrs in known_utrs.iteritems():
        for regions in key_known_utrs.itervalues():
            regions.sort()
            introns = izip( zip( *regions )[1][:-1], zip( *regions )[0][1:] )
            for intron in introns:
                known_introns[(chrm, strand)]['UTR']['intron'].add( intron )
                known_introns[(chrm, strand)]['UTR']['donor'][ \
                    intron[0] ].add( intron[1] )
                known_introns[(chrm, strand)]['UTR']['acceptor'][ \
                    intron[1] ].add( intron[0] )
    
    return known_introns

def parse_discovered_orfs( orfs_fp ):
    """Parse gtf file with discovered orfs
    """
    disc_orfs = defaultdict( lambda: defaultdict( lambda: defaultdict(list) ) )
    disc_orf_lines = defaultdict( str )
    for line in orfs_fp:
        gtfl = parse_gtf_line( line )
        if gtfl == None or gtfl.feature not in \
                ( 'CDS', '5UTR', '3UTR' ): continue
        
        disc_orfs[ (gtfl.region.chr, gtfl.region.strand) ][ gtfl.trans_id ][ \
            gtfl.feature ].append( (gtfl.region.start, gtfl.region.stop) )
        disc_orf_lines[ (gtfl.region.chr, gtfl.region.strand, gtfl.trans_id) ] += line
    
    return disc_orfs, disc_orf_lines

def build_objects( orfs_fp, ann_fp, jns_fp, out_prefix ):
    as_orfs_fp = open( out_prefix + '.alternative_spliced_orfs.gtf', 'w' )
    shifted_orfs_fp = open( out_prefix + '.shifted_orfs.gtf', 'w' )
    
    if VERBOSE: print >> sys.stderr, 'parsing junctions'
    all_jns = parse_jns( jns_fp )
    jns_fp.close()
    
    if VERBOSE: print >> sys.stderr, 'parsing discovered orfs'
    disc_orfs, disc_orf_lines = parse_discovered_orfs( orfs_fp )
    orfs_fp.close()
    
    if VERBOSE: print >> sys.stderr, 'parsing known orfs'
    known_introns = parse_known_orfs( ann_fp )
    ann_fp.close()
    
    return disc_orfs, disc_orf_lines, known_introns, all_jns, \
        as_orfs_fp, shifted_orfs_fp

def parse_arguments():
    global LOW_CDS_JN_RATIO_THRESH
    global HIGH_CDS_JN_RATIO_THRESH
    global ACCEPT_ALL_UTRS
    global UTR_JN_RATIO_THRESH
    
    import argparse
    
    parser = argparse.ArgumentParser(
        description = 'Find discovered alternatively splice ORFs that do not ' + \
            'shift the frame of the known ORF.' )
    parser.add_argument(
        'orfs_gtf', type=file,
        help='GTF file, Alternative Splice class, with CDS and UTR regions.' )
    parser.add_argument(
        'annotation', type=file,
        help='GTF file of annotated orfs assocaited with genes and transcripts.')
    parser.add_argument(
        'junctions', type=file,
        help='Gff file of junctions with scores.')
    parser.add_argument(
        'ki_orfs_gtf',
        help='GTF file, Known Isoform class, with CDS and UTR regions.' )    
    
    parser.add_argument( 
        '--filter-utrs', default=False, action='store_true',
        help='Whether or not to filter utr introns. default: False')
    parser.add_argument( 
        '--utr-jn-ratio', type=float, default=UTR_JN_RATIO_THRESH,
        help='Threshold ratio of novel utr intron with known intron. '
        'default: %(default)f ')
    
    parser.add_argument( 
        '--low-cds-jn-ratio', type=float, default=LOW_CDS_JN_RATIO_THRESH,
        help='Threshold ratio of novel in-frame intron with known intron. '
        'default: %(default)f ')
    parser.add_argument( 
        '--high-cds-jn-ratio', type=float, default=HIGH_CDS_JN_RATIO_THRESH,
        help='Threshold ratio of novel out-of-frame intron with known intron.'
        'default: %(default)f ')
    
    parser.add_argument(
        '--out-prefix', '-o',
        help='Prefix of output file. (default: gtf)')
    parser.add_argument( 
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    ACCEPT_ALL_UTRS = not args.filter_utrs
    UTR_JN_RATIO_THRESH = args.utr_jn_ratio
    LOW_CDS_JN_RATIO_THRESH = args.low_cds_jn_ratio
    HIGH_CDS_JN_RATIO_THRESH = args.high_cds_jn_ratio
    
    # create default if no prefix provided or if same as gtf filename
    out_prefix = args.out_prefix if args.out_prefix != None else \
        os.path.basename( args.orfs_gtf.name )
    
    # set flag args
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.orfs_gtf, args.annotation, args.junctions, args.ki_orfs_gtf, out_prefix

def main():
    disc_orfs_fp, ann_fp, jns_fp, ki_fname, out_prefix = parse_arguments()
    disc_orfs, disc_orf_lines, known_introns, all_jns, \
        as_orfs_fp, shifted_orfs_fp = \
        build_objects( disc_orfs_fp, ann_fp, jns_fp, out_prefix )
    
    if VERBOSE: print >> sys.stderr, 'filtering orfs with shifted frames'
    as_orfs, shifted_orfs = filter_shifted_trans( 
        disc_orfs, known_introns, all_jns )
    
    if VERBOSE: print >> sys.stderr, 'writing out filtered transcripts'
    write_novel_orfs( as_orfs, as_orfs_fp, disc_orf_lines )
    write_shifted_orfs( shifted_orfs, shifted_orfs_fp, disc_orf_lines )
    
    write_all_filtered_orfs( out_prefix, ki_fname )
    
    return

if __name__ == '__main__':
    if DO_PROFILE:
        import cProfile
        cProfile.run( 'main()' )
    else:
        main()
