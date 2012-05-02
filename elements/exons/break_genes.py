VERBOSE = False
MIN_VERBOSE = True

MIN_CVRG_REGION_SIZE = 200
CVRG_FRACTION = 0.95
FOLD_LEN_THRSH = 5

import sys
import os
import numpy

from build_genelets import cluster_all_exons, write_all_clusters
from gene_models import parse_junctions_file

sys.path.append( os.path.join( os.path.dirname(__file__), '..', 'file_types' ) )
from wiggle import wiggle
from clustered_exons_file import parse_clustered_exons_file
from exons_file import parse_exons_file

sys.path.append( os.path.join( os.path.dirname(__file__), 'get_elements' ) )
from get_tss_exons import convert_to_genomic_intervals
from get_exons import write_exons

def cluster_and_write_exons( int_exons, tss_exons, tes_exons, jns, out_fp ):
    if jns != None:
        all_exons = {}
        keys = set( int_exons ).union( tss_exons ).union( tes_exons )
        for key in keys:
            key_int_exons = int_exons[key] if key in int_exons else []
            key_tss_exons = [ tuple(item) for item in tss_exons[key] ] if key in tss_exons else []
            key_tes_exons = [ tuple(item) for item in tes_exons[key] ] if key in tes_exons else []
            all_exons[key] = list( set( key_int_exons ).union( \
                    key_tss_exons ).union( key_tes_exons ) )
        all_genelets = cluster_all_exons( all_exons, jns )
        write_all_clusters( all_genelets, out_fp )
    else:
        exons = convert_to_genomic_intervals( int_exons )
        write_exons( exons, out_fp )
    
    return

def remove_gene_merge_exons( clustered_exons, tss_exons, tes_exons, all_cvrg, key, log_fp ):
    """Remove exons from a cluster that may be causing a gene merge
    1) Identify possible merge exons
         - shares one boundary with a TSS exon and one with a TES exon
    2) Confirm merge exons
         a) find lowest coverage region within candidate exon
         b) find region from either side of low_cvrg_region which contains 
              CVRG_FRACTION of the RNAseq coverage
         c) Check if length of region on either side is small enough to
              confirm that this is merge exon
    3) Remove all exons overlapping any confirmed merge exon
    """
    chrm, strand = key
    
    def get_internal_coords( exons, is_tss ):
        """Return the list of all internal coordinates from start or stop exons
        """
        internal_coords = []
        int_index = 1 if (strand == '+' and is_tss) or (strand == '-' and not is_tss) else 0
        for exon in exons:
            internal_coords.append( exon[int_index] )
        
        return internal_coords
    
    tss_int_coords = get_internal_coords( tss_exons, True )
    tes_int_coords = get_internal_coords( tes_exons, False )
    
    upstrm_coords = tes_int_coords if strand == '+' else tss_int_coords
    dnstrm_coords = tss_int_coords if strand == '+' else tes_int_coords
    
    def get_candidate_merge_exons( cluster ):
        """Get candidate merge exons 
        def: exons that share one bndry with a TSS the other boundary with a TES
        """
        candidate_merge_exons = []
        for start, stop in cluster:
            if start in upstrm_coords and stop in dnstrm_coords:
                candidate_merge_exons.append( (start, stop) )
        
        return candidate_merge_exons
        
    def confirm_merge_exons( candidate_merge_exons ):
        def is_gene_merge_exon( start, stop ):
            # first get the region within the start stop with the lowest coverage
            if stop - start + 1 <= MIN_CVRG_REGION_SIZE: 
                log_fp.write( 'Region is too small.\n' )
                return False
            
            region_cumsum = numpy.append( 0, numpy.cumsum( all_cvrg[key][start:stop+1] ) )
            min_cvrg_total = None
            min_cvrg_region = None
            for offset in xrange( stop - start + 1 - MIN_CVRG_REGION_SIZE ):
                if min_cvrg_region == None or \
                        region_cumsum[offset+MIN_CVRG_REGION_SIZE] - region_cumsum[offset] < min_cvrg_total:
                    min_cvrg_total = region_cumsum[offset+MIN_CVRG_REGION_SIZE] - region_cumsum[offset]
                    min_cvrg_region = (start + offset, start + offset + MIN_CVRG_REGION_SIZE)
            log_fp.write( '\tlow cvrg region:' + str(min_cvrg_region[0]) + '-' + str(min_cvrg_region[1]) + '\n' )
            
            # get region on either side of the low cvrg region that contains CVRG_FRACTION
            # then determine if that region is *small enough* to confirm that this is a merge exon
            upstrm_region = all_cvrg.get_region_w_fraction_cvrg( \
                key, start, min_cvrg_region[0]-1, CVRG_FRACTION, True )
            upstrm_len = min_cvrg_region[0] - start
            upstrm_low_cvrg_len = min_cvrg_region[0] - upstrm_region[1]
            log_fp.write( '\tupstream region:' + str(upstrm_region[0]) + '-' + str(upstrm_region[1]) + '\n' )
            if upstrm_low_cvrg_len < upstrm_len * (1-CVRG_FRACTION) * FOLD_LEN_THRSH:
                return False
            
            dnstrm_region = all_cvrg.get_region_w_fraction_cvrg( \
                key, min_cvrg_region[1]+1, stop, CVRG_FRACTION, False )
            dnstrm_len = stop - min_cvrg_region[1]
            dnstrm_low_cvrg_len =  dnstrm_region[0] - min_cvrg_region[1]
            log_fp.write( '\tdownstream region:' + str(dnstrm_region[0]) + '-' + str(dnstrm_region[1]) + '\n' )
            if dnstrm_low_cvrg_len < dnstrm_len * (1-CVRG_FRACTION) * FOLD_LEN_THRSH:
                return False
            
            return True
        
        # attempt to confirm each candidate merge exon
        confirmed_merge_exons = []
        log_fp.write( '#'*80 + '\n' )
        for start, stop in candidate_merge_exons:
            log_fp.write( 'chr' + chrm + ':' + str(start) +'-'+ str(stop) + ' ' + strand + '\n' )
            
            if is_gene_merge_exon( start, stop ):
                log_fp.write( '^'*15 + 'CONFIRMED' + '^'*15 + '\n' )
                confirmed_merge_exons.append( (start, stop) )
            
            log_fp.write( '\n' )
        
        return confirmed_merge_exons
    
    def get_nonoverlapping_exons( confirmed_merge_exons, cluster ):
        """Get all exons that *don't* overlap a merge exon
        """
        non_merge_exons = []
        for start, stop in cluster:
            # if this exon *doesn't* overlap any merge exons
            if all( start > merge_exon[1] or stop < merge_exon[0] \
                        for merge_exon in confirmed_merge_exons ):
                non_merge_exons.append( (start, stop) )
        
        return non_merge_exons
    
    
    filtered_exons = []
    for cluster in clustered_exons:
        # change cluster exons only if at least one gene merge exon is found
        new_cluster = cluster
        candidate_merge_exons = get_candidate_merge_exons( cluster )
        if len( candidate_merge_exons ) > 0:
            confirmed_merge_exons = confirm_merge_exons( candidate_merge_exons )
            if len( confirmed_merge_exons ) > 0:
                new_cluster = get_nonoverlapping_exons( confirmed_merge_exons, cluster )
        
        filtered_exons.extend( new_cluster )
    
    return filtered_exons

def remove_all_gene_merge_exons( clustered_exons, tss_exons, tes_exons, all_cvrg, log_fp ):
    """Wrapper for remove_gene_merge_exons
    """
    all_filtered_exons = {}
    keys = set( clustered_exons ).intersection( \
        tss_exons ).intersection( tes_exons ).intersection( all_cvrg )
    for key in keys:
        chrm, strand = key
        all_filtered_exons[ key ] = remove_gene_merge_exons( \
            clustered_exons[key], tss_exons[key], tes_exons[key], all_cvrg, key, log_fp )
    
    return all_filtered_exons

def build_objs( cluster_exons_fp, tss_exons_fp, tes_exons_fp, wig_fps, chrm_sizes_fp, jns_fp ):
    clustered_exons = parse_clustered_exons_file( cluster_exons_fp )
    cluster_exons_fp.close()
    tss_exons = parse_exons_file( tss_exons_fp )
    tss_exons_fp.close()
    tes_exons = parse_exons_file( tes_exons_fp )
    tes_exons_fp.close()
    
    all_cvrg = wiggle( chrm_sizes_fp )
    for wig_fp in wig_fps:
        all_cvrg.load_data_from_fname( wig_fp )

    jns = parse_junctions_file( jns_fp ) if jns_fp != None else None
    
    return clustered_exons, tss_exons, tes_exons, all_cvrg, jns

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(description='Determine valid transcripts.')
    parser.add_argument( 'clustered_exons_gff', type=file, \
        help='GTF of exons with corresponding gene_ids.')
    parser.add_argument( 'tss_exons_gff', type=file, \
        help='A GFF file containing valid transcription stop site exons. ' )
    parser.add_argument( 'tes_exons_gff', type=file, \
        help='A GFF file containing valid transcription stop site exons.')
    parser.add_argument( 'chrm_sizes_fname', type=file, \
        help='File with chromosome names and sizes.')
    parser.add_argument( 'wig_fps', type=file, nargs='+', \
        help='WIG files of coverge.')
    
    parser.add_argument( '--junctions_gff', '-j', type=file, \
        help='A gff file containing valid introns.')
    
    parser.add_argument( '--out_fname', '-o', default='',\
                             help='Output file name. (default: stdout)')
    parser.add_argument( '--log_fname', '-l',\
                             help='Output log file name. (default: sterr)')
    parser.add_argument( '--verbose', '-v', default=False, action= 'store_true', \
                             help='Whether or not to print status information.')
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    out_fp = open( args.out_fname, "w" ) if args.out_fname else sys.stdout

    log_fp = open( args.log_fname, 'w' ) if args.log_fname else sys.stderr
    
    return args.clustered_exons_gff, args.tss_exons_gff, args.tes_exons_gff, \
        args.wig_fps, args.chrm_sizes_fname, args.junctions_gff, out_fp, log_fp

def main():
    clustered_exons_fp, tss_exons_fp, tes_exons_fp, \
        wig_fps, chrm_sizes_fp, jns_fp, out_fp, log_fp = parse_arguments()
    
    if MIN_VERBOSE:
        print >> sys.stderr, 'Parsing input...'
    clustered_exons, tss_exons, tes_exons, all_cvrg, jns = build_objs( \
        clustered_exons_fp, tss_exons_fp, tes_exons_fp, wig_fps, chrm_sizes_fp, jns_fp )
    
    if MIN_VERBOSE:
        print >> sys.stderr, 'Detecting gene merges...'
    exons = remove_all_gene_merge_exons( clustered_exons, tss_exons, tes_exons, all_cvrg, log_fp )
    
    if MIN_VERBOSE:
        print >> sys.stderr, 'Writing output...'
    cluster_and_write_exons( exons, tss_exons, tes_exons, jns, out_fp )
    
    return

if __name__ == "__main__":
    main()

