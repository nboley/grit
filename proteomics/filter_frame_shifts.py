import sys, os 

DO_PROFILE = False

APPLY_LENGTH_THRESH = True
LENGTH_RATIO_THRESH = 0.5

import re
import numpy

from collections import defaultdict
from itertools import izip

sys.path.append( os.path.join(os.path.dirname(__file__), "..", "file_types") )
from gtf_file import parse_gtf_line


""" NATHAN: This function defines a frame shift event:
disc_splice is a tuples containing the last pos in an exon and 
the first pos in the next exon

known_orf is a 2D numpy.array (for searchsorted) with the exon coordinates

I have tried to start this function that determines if an individual splice 
event is a frameshift, but I have to go now.  I am pretty sure that
everything else should work, but it is not debugged. I also have not written the
 print functions, but those are straightforward as I am storing the lines from
the original gtf files. I can work on this in Monday if you want to wait to 
add it to grit and I can try to run it by hand if we have a good run through orf
finder on Monday.
"""
def is_frame_shift_event( disc_splice, known_orf ):
    # find the index of the index that might match the splice donor site
    known_e_disc_donor_i = known_orf[:,1].searchsorted( disc_splice[0] )
    donor_site_match = known_orf[known_e_disc_donor_i][1] == disc_splice[0]
    
    # find the index of the index that might match the splice acceptor site
    known_e_disc_acceptor_i = known_orf[:,0].searchsorted( disc_splice[1] )
    acceptor_site_match = known_orf[known_e_disc_donor_i][0] == disc_splice[1]
    
    if donor_site_match and acceptor_site_match:
        # if this is a matched junction
        if known_e_disc_donor_i == known_e_disc_acceptor_i - 1:
            return False
        
        # check if removing the intervening exon(s) changes the frame
        num_bases_shifted = sum( stop - start + 1 for start, stop in known_orf[ 
                known_e_disc_donor_i+1 : known_e_disc_acceptor_i ] )
        if num_bases_shifted % 3 != 0: 
            return True
    
    # if the donor site matches, but the acceptor site is novel
    elif donor_site_match:
        num_bases_shifted = sum( stop - start + 1 for start, stop in known_orf[
                known_e_disc_donor_i+1 : known_e_disc_acceptor_i ] )
        # add or subtract the number of bases shifted on the downstream end
        num_bases_shifted += disc_splice[1] - known_orf[known_e_disc_acceptor_i][0]
        
        if num_bases_shifted % 3 != 0:
            return True
    
    elif acceptor_site_match:
        # if this is donor site extention into an intron
        if known_e_disc_donor_i == known_e_disc_acceptor_i:
            if disc_splice[0] - known_orf[ known_e_disc_acceptor_i ][1] % 3 == 0:
                return False
            else:
                return True
        
        num_bases_shifted = -1 * sum( 
            stop - start + 1 for start, stop in known_orf[ 
                known_e_disc_donor_i : known_e_disc_acceptor_i ] )
        # add the number of bases shifted on the upstream end
        num_bases_shifted += disc_splice[0] - known_orf[known_e_disc_donoe_i-1][0]
        
        if num_bases_shifted % 3 != 0:
            return True
    
    # if both the donor and acceptor are mismatched
    else:
        pass
    
    return False

def is_frame_shifted( disc_orf, known_orf ):
    disc_orf_int_coords = [pos for exon in disc_orf for pos in exon][1:-1]
    for disc_splice in izip( disc_orf[0::2], disc_orf[1::2] ):
        if is_frame_shift_event( disc_splice, known_orf ):
            ### THIS IS WHERE WE WANT TO CHECK FOR HIGHER JUNCTION THRESHOLD ###
            # if not passes_higher_jn_thresh( disc_orf, known_match, junction counts ):
            return True
    
    return False

def filter_shifted_orfs( disc_orfs, known_orfs ):
    as_orfs = []
    shifted_orfs = []
    
    for (chrm, strand), orfs in disc_orfs.iteritems():
        key_known_orfs = known_orfs[(chrm, strand)]
        for (gene_id, trans_id, known_trans_id), trans_data in orfs.iteritems():
            disc_orf = sorted( trans_data[ 'CDS' ] )
            known_orf = key_known_orfs[ known_trans_id ]
            
            # check for cds length here and filter those that are 
            # below LENGTH_RATIO_THRESH
            if APPLY_LENGTH_THRESH:
                disc_len = sum( stop - start + 1 for start, stop in disc_orf )
                known_len = sum( stop - start + 1 for start, stop in known_orf )
                if disc_orf / float( known_len ) < LENGTH_RATIO_THRESH:
                    shifted_orfs.append( trans_data[ 'lines' ] )
            
            # check if the frame is shifted
            if not is_frame_shifted( disc_orf, known_orf ):
                as_orfs.append( trans_data[ 'lines' ] )
            else:
                shifted_orfs.append( trans_data[ 'lines' ] )
    
    return as_orfs, shifted_orfs

def parse_known_trans( ann_fp ):
    raw_known_orfs = defaultdict( lambda: defaultdict( list ) )
    for line in ann_fp:
        gtfl = parse_gtf_line( line )
        if gtfl == None or gtfl.feature != 'CDS': continue
        
        raw_known_orfs[(gtfl.region.chr, gtfl.region.strand)][ 
            (gtfl.trans_id)].append( (gtfl.region.start, gtfl.region.stop) )
    
    known_orfs = defaultdict( dict )
    for (chrm, strand), transcripts in known_orfs.iteritems():
        for trans_id, cdss in transcripts.iteritems():
            known_orfs[(chrm, strand)][ trans_id ] = \
                numpy.array( sorted( cdss ))
    
    return known_orfs

def parse_discovered_orfs( orfs_fp ):
    """Parse gtf file with discovered orfs
    """
    known_match_pat = re.compile( '"AS:(.+)"' )
    disc_orfs = defaultdict( lambda: defaultdict( lambda: defaultdict(list) ) )
    
    for line in orfs_fp:
        gtfl = parse_gtf_line( line )
        if gtfl == None or gtfl.feature not in \
                ( 'CDS', '5UTR', '3UTR' ): continue
        main_key = (gtfl.region.chr, gtfl.region.strand)
        
        known_match = known_match_pat.match( 
            gtfl.meta_data[ 'orf_class' ] ).group(1)
        # Store the cds region if this is a cds
        if gtfl.feature == 'CDS':
            disc_orfs[ main_key ][ \
                ( gtlf.gene_id, gtfl.trans_id, known_match ) ][ \
                'CDS' ].append( ( gtfl.region.start, gtfl.region.stop ) )
        disc_orfs[ main_key ][ ( gtfl.gene_id, gtfl.trans_id, known_match ) ][ 
            'lines' ].append( line )
    
    return disc_orfs

def build_objects( orfs_fp, ann_fp, out_prefix ):
    as_orfs_fp = open( out_prefix + '.alternative_spliced_orfs.gtf', 'w' )
    shifted_orfs_fp = open( out_prefix + '.shifted_orfs.gtf', 'w' )
    
    if VERBOSE: print >> sys.stderr, 'parsing discovered orfs'
    disc_orfs = parse_discovered_orfs( orfs_fp )
    orfs_fp.close()
    
    if VERBOSE: print >> sys.stderr, 'parsing known orfs'
    known_orfs = parse_known_orfs( ann_fp )
    ann_fp.close()
    
    return disc_orfs, known_orfs, as_orfs_fp, shifted_orfs_fp

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(
        description = 'Find discovered alternatively splice ORFs that do not ' + \
            'shift the frame of the known ORF.' )
    parser.add_argument(
        'orfs_gtf', type=file,
        help='GTF file with CDS and UTR regions.' )
    parser.add_argument(
        'annotation', type=file,
        help='GTF file of annotated orfs assocaited with genes and transcripts.')
    
    parser.add_argument(
        '--out_prefix', '-o',
        help='Prefix of output file. (default: gtf)')
    parser.add_argument( 
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    # create default if no prefix provided or if same as gtf filename
    out_prefix = args.out_prefix if args.out_prefix != None else \
        os.path.basename( args.orfs_gtf.name )
    
    # set flag args
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.orfs_gtf, args.annotation, out_prefix

def main():
    disc_orfs_fp, ann_fp, out_prefix = parse_arguments()
    disc_orfs, known_orfs, as_orfs_fp, shifted_orfs_fp = \
        build_objects( disc_orfs_fp, ann_fp, out_prefix )
    
    if VERBOSE: print >> sys.stderr, 'filtering orfs with shifted frames'
    as_orfs, shifted_orfs = get_ann_matched_trans( 
        disc_orfs, known_trans )
    
    write_novel_orfs( as_orfs, as_orfs_fp )
    write_frameshifted_orfs( shifted_orfs, shifted_orfs_fp )
    
    return

if __name__ == '__main__':
    if DO_PROFILE:
        import cProfile
        cProfile.run( 'main()' )
    else:
        main()
