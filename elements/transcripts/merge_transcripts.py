import sys, os

import numpy
from collections import defaultdict, OrderedDict
from itertools import izip, repeat

sys.path.append( os.path.join( os.path.dirname( __file__ ),
                               "..", "..", "file_types", "fast_gtf_parser" ) )
from gtf import load_gtfs

sys.path.append( os.path.join( os.path.dirname( __file__ ), 
                               "..", "..", "file_types" ) )
from gtf_file import iter_gtf_lines
from genomic_intervals import GenomicInterval

sys.path.append( os.path.join( os.path.dirname( __file__ ), "..", "exons" ) )
from build_genelets import build_edges_hash_table, \
    build_genelets_from_exons_and_edges, cluster_all_exons

VERBOSE = True

def build_gtf_lines( chrm, strand, gene_name, trans_name, exon_bndries ):
    feature="exon"
    source="merge_transcripts"
    gene_name_iter = repeat( gene_name )
    trans_name_iter = repeat( trans_name )
    
    # create an iterator over every exon
    exons_iter = enumerate( izip( exon_bndries[0::2], exon_bndries[1::2] ) )
    regions_iter = ( GenomicInterval( chrm, strand, start, stop )
                     for exon_num, (start, stop) in exons_iter )
    
    return iter_gtf_lines( regions_iter, gene_name_iter, trans_name_iter,
                           feature=feature, source=source )

def write_all_transcripts( clustered_trans, ofp, sources_fp, transcript_maps ):
    # iterate through unique transcripts, and write them to the output file
    for (chrm, strand), genes in clustered_trans.iteritems():
        for gene_id, transcripts in genes.iteritems():
            gene_name = "merged_cluster_%i" % gene_id
            for trans_id, exons in enumerate( transcripts ):
                trans_name = gene_name + "_%i" % trans_id
                
                if sources_fp != None:
                    sources = [ os.path.basename(source) for source in
                                transcript_maps[ ( chrm, strand ) ][ exons ] ]
                    source_str = ",".join( sources )
                    sources_fp.write( '\t'.join(
                            (gene_name, trans_name, source_str ) ) + '\n' )
                
                gtf_lines = build_gtf_lines(
                    chrm, strand, gene_name, trans_name, exons )
                
                ofp.write( "\n".join( gtf_lines ) + "\n" )
    
    ofp.close()
    if sources_fp != None:
        sources_fp.close()
    
    return

def get_exons_w_shared_coords( exons, key ):
    """Group exons by those that contain the same start or stop coordinate
    This is single linkage clustering
    
    Hueristic:
    For each set of overlapping exons
    Base) Store the first exon in the first index of grp_conn_exons and track
          the index at which it is stored in separate start and stop coordinate.
    Recursion)
        1) If an exon's start and stop coords are not in the grp_starts and 
           grp_stops dicts, resp'ly, append the exon to grp_conn_exons.
        2) Add exons to the group in which they belong by searching grp_starts
           and grp_stops dictionaries and adjust appropriate dicts accordingly.
        3) If the start and stop coords of an exon point to different groups
           then those groups must be merged and the grp_starts and grp_stops 
           dicts must be updated for each exon in the moved group.
    End) Exons are sorted by start coord so when an exon does not overlap the
         current group it and any exon after it can no longer share a start or
         stop coord with the current group. grp_conn_exons are dumped and the
         grp_starts and grp_stops dicts are reset.
    """
    # sort exons by start coord
    exons = exons[ exons[:,0].argsort() ]
    
    # structures for tracking the current group of exons
    # when exons can not longer share a coord these structs are reset
    # in order to avoid large dict key lookups
    grp_conn_exons = []
    # key == coord; value == current grp_conn_exons array index
    coord_to_group = {}
    
    def add_exon_to_grp( start, stop ):
        def merge_exon_groups():
            # add exon to the stop indexed position
            grp_conn_exons[ coord_to_group[stop] ].append( (start, stop) )
            
            # add all exons from the start indexed group
            grp_conn_exons[ coord_to_group[stop] ].extend( 
                grp_conn_exons[ coord_to_group[start] ] )
            
            # remember the start group index as it will change in next loop
            start_grp_index = coord_to_group[start]
            
            # point all coords that were just moved to the stop indexed 
            # array to the appropriate grp_conn_exons array index
            for exon in grp_conn_exons[ coord_to_group[start] ]:
                coord_to_group[ exon[0] ] = coord_to_group[stop]
                coord_to_group[ exon[1] ] = coord_to_group[stop]
            
            # remove contents of the start list *cleaned it up later*
            # can't delete as list indices must be maintained
            grp_conn_exons[ start_grp_index ] = []
            
            return
        
        # if start coord has been seen previously
        if start in coord_to_group:
            # both start and stop coords have been seen previously
            if stop in coord_to_group:
                # if the start and stop coord groups are already the same
                if coord_to_group[start] == coord_to_group[stop]:
                    grp_conn_exons[ coord_to_group[start] ].append( (start, stop) )
                # start and stop groups are different and must be merged
                else:
                    merge_exon_groups()
            # only start coord has been seen previously
            else:
                grp_conn_exons[ coord_to_group[start] ].append( (start, stop) )
                coord_to_group[ stop ] = coord_to_group[ start ]
        # if the stop coord, but not the start coord, has been seen previously
        elif stop in coord_to_group:
            grp_conn_exons[ coord_to_group[stop] ].append( (start, stop) )
            coord_to_group[ start ] = coord_to_group[ stop ]
        # if neither the start or stop coord have been seen yet
        else:
            coord_to_group[start] = len(grp_conn_exons)
            coord_to_group[stop] = len(grp_conn_exons)
            grp_conn_exons.append( [ (start, stop), ] )
        
        return
    
    
    # check when we can reset the grp structs
    grp_max_stop = -1
    # store all groups of connected exons
    all_connected_exons = []
    for start, stop in exons:
        # if this exon is past the last stop of the current group
        if grp_max_stop < start:
            # dump all exon groups that are not empty
            all_connected_exons.extend(
                [ exons for exons in grp_conn_exons if len( exons ) > 0 ] )
            # reset all values when starting a new group
            grp_conn_exons = []
            coord_to_group = {}
        
        # add exon information to the current group of exons
        add_exon_to_grp( start, stop )
        # adjust grp_max_stop to maintian current grp max position
        grp_max_stop = max( grp_max_stop, stop )
    
    # add final exon groups that are not empty
    all_connected_exons.extend( 
        [ exons for exons in grp_conn_exons if len( exons ) > 0 ] )
    
    # return a dict of a group id pointing to a group of exons
    return dict( zip( range(len(all_connected_exons)), all_connected_exons ) )

def get_connected_components( exons, jns ):
    all_connected_exons = {}
    keys = set( exons ).intersection( jns )
    for key in sorted( keys ):
        
        key_jns = numpy.array( jns[key] )
        
        # filtered_exons = filter_exons( exons[key], key_jns )
        filtered_exons = numpy.array( exons[ key ] )
        
        exons_w_shared_coords = get_exons_w_shared_coords( filtered_exons, key )
        
        if len(exons_w_shared_coords) > 0:
            # find the hash table of exons ( nodes ) and junctions ( edges )
            edges = build_edges_hash_table( exons_w_shared_coords, key_jns )
            
            # build the genelets by connected components
            genelets = build_genelets_from_exons_and_edges( 
                exons_w_shared_coords, edges )
            
            all_connected_exons[key] = genelets
    
    return all_connected_exons

def get_exons_and_junctions( transcript_maps ):
    # extract unique exons and jns sets from transcripts for clustering
    all_exons = {}
    all_jns = {}
    for (chrm, strand), transcripts in transcript_maps.iteritems():
        key_exons = set()
        key_jns = set()
        for exons in transcripts:
            # skip single exon genes as they will not cluster
            if len(exons) <= 2: continue
            
            key_exons.update( izip( exons[::2], exons[1::2] ) )
            key_jns.update( [(start+1, stop-1) for start, stop in 
                                 izip( exons[1::2], exons[2::2] ) ] )
        
        # convert sets to list for clustering
        all_exons[ (chrm, strand) ] = list( key_exons )
        all_jns[ (chrm, strand) ] = list( key_jns )
    
    return all_exons, all_jns

def create_exons_map( exon_groups ):
    # reverse exon groups so that an exon points to a unique id
    exons_map = {}
    # cluster id (gene_id)
    c_id = 1
    for (chrm, strand), key_exon_groups in exon_groups.iteritems():
        for exons in key_exon_groups:
            for exon in exons:
                exons_map[ (chrm, strand, exon) ] = c_id
            c_id += 1
    
    return exons_map

def get_exon_clusters( transcript_maps, cluster_by_conn ):
    # get all unique exons and junctions from all transcripts
    all_exons, all_jns = get_exons_and_junctions( transcript_maps )
    
    if cluster_by_conn:
        # get the exons connected by junctions *but NOT by overlapping coords*
        exon_groups = get_connected_components( all_exons, all_jns )
    else:
        # get exons connected by overlap *OR* connected components
        exon_groups = cluster_all_exons( all_exons, all_jns )
    
    # convert maps dict
    exons_map = create_exons_map( exon_groups )
    
    return exons_map

def cluster_transcripts( transcript_maps, cluster_by_conn ):
    # cluster all exons from transcripts
    exons_map = get_exon_clusters( transcript_maps, cluster_by_conn )
    
    # identify cluster which each transcript belongs
    clustered_trans = defaultdict( lambda: defaultdict( list ) )
    for (chrm, strand), all_trans_exons in transcript_maps.iteritems():
        for exons in all_trans_exons.iterkeys():
            try:
                # try to get cluster id from exons map
                key = (chrm, strand, tuple(exons[:2]))
                gene_id = exons_map[ key ]
            # if the key does not exist then this must be a single exon gene
            except KeyError:
                # give single exon genes their own cluster regardless of 
                # possible overlap
                assert len( exons ) == 2
                gene_id = c_id
                c_id += 1
            
            # add the gene_id to the transcripts meta_data
            clustered_trans[(chrm, strand)][gene_id].append( exons )
    
    return clustered_trans

def build_unique_transcripts( transcriptomes, gtf_fnames ):
    # build the sets of transcripts
    tmp_transcript_maps = defaultdict( lambda: defaultdict( list ) )
    for source, transcriptome in izip( gtf_fnames, transcriptomes ):
        for cluster_id, chrm, strand, start, stop, transcripts in transcriptome:
            for trans_id, exons in transcripts:
                # skip single  exon transcripts and transcripts that 
                # have an odd number of boundaries
                if len( exons ) <= 2: continue
                if len( exons ) % 2 != 0: continue
                
                tmp_transcript_maps[ ( chrm, strand ) ][ \
                    tuple( exons[1:-1] ) ].append( \
                    (exons[0], exons[-1], source) )
    
    # get the longest external coordinate for each internal structure
    transcript_maps = defaultdict( dict )
    for key, transcripts in tmp_transcript_maps.iteritems():
        for internal_exons, ends_and_sources in transcripts.iteritems():
            start_coord = min( zip( *ends_and_sources )[0] )
            stop_coord = max( zip( *ends_and_sources )[1] )
            exons = tuple(
                [start_coord,] + list( internal_exons ) + [stop_coord,])
            # uniquify the list of sources
            transcript_maps[ key ][ exons ] = \
                list(set(zip( *ends_and_sources )[2]))
    
    return transcript_maps

def parse_arguments():
    import argparse
    desc = 'Merge transcripts.'
    parser = argparse.ArgumentParser(description=desc)    
    parser.add_argument(
        "gtfs", type=str, nargs="+",
        help='GTF file contraining transcripts. (i.e. from build_transcripts)' )
    
    parser.add_argument(
        '--sources-fname', type=argparse.FileType( 'w' ),
        help='File name to write the sources for each transcript. '
        + 'Default: Do not write out sources map.')
    
    parser.add_argument(
        '--out-fname', '-o', type=argparse.FileType('w'), default=sys.stdout,
        help='Output file. default: stdout')
    
    parser.add_argument(
        '--cluster-by-overlap', default=True, action='store_false', 
        dest='cluster_by_conn',
        help='Whether or not to print status information.')
    
    parser.add_argument(
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
 
    parser.add_argument(
        '--threads', '-t', type=int, default=1,
        help='The number of threads to use.')
   
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.gtfs, args.out_fname, args.sources_fname, \
        args.cluster_by_conn, args.threads

def main():
    gtf_fnames, ofp, sources_fp, cluster_by_conn, n_threads = parse_arguments()
    
    # load all transcripts file and store corresponding data
    transcriptomes = load_gtfs( gtf_fnames, n_threads )
    
    # build maps from unique transcripts to a list of sources (gtf file names)
    transcript_maps = build_unique_transcripts( transcriptomes, gtf_fnames )
    
    # use clustered exons to recluster transcripts
    clustered_trans = cluster_transcripts( transcript_maps, cluster_by_conn )
    
    write_all_transcripts( clustered_trans, ofp, sources_fp, transcript_maps )
    
    return

if __name__ == "__main__":
    main()
