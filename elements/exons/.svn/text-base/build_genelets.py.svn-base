VERBOSE = False

import sys
import os

import numpy
from collections import defaultdict
from operator import itemgetter

# import slide modules
sys.path.append( os.path.join(os.path.dirname(__file__), "../../file_types") )
from exons_file import parse_exons_file
from junctions_file import parse_junctions_file
sys.path.append( os.path.join(os.path.dirname(__file__), "../") )
from extract_elements_from_transcripts import get_elements_from_gene, load_gtf

def write_all_clusters( all_genelets, out_fp, track_name="clustered_exons" ):
    cluster_id = 1
    out_fp.write( "track name=" + track_name + "\n" )
    for (chrm, strand), clusters in all_genelets.iteritems():
        gff_tempate_line = '\t'.join( \
            ( 'chr' + chrm, 'build_genelets', 'exon', '{0}', '{1}', \
                  '.', strand, '.', '{2}' ) )
        
        def create_gff_cluster_lines( cluster ):
            for exon_num, (start, stop) in enumerate(cluster):
                yield gff_tempate_line.format( start, stop, cluster_id )
            return
        
        
        for cluster in clusters:
            cluster_lines = create_gff_cluster_lines( cluster )
            out_fp.write( '\n'.join(cluster_lines) + '\n' )
            cluster_id += 1
    
    return

def cluster_overlapping_exons( exons ):
    """Cluster overlapping exons by sorting them comparing neighboring exons.
    """
    exons = exons[ exons[:,0].argsort() ]
    exon_groups = []
    
    grp_stop = -1
    current_grp = []
    for start, stop in exons:
        # if start is > the current group stop add the current group and 
        # start a new group
        if start > grp_stop and current_grp:
            exon_groups.append( current_grp )
            current_grp = []
        # if this stop extends the group stop do so
        if stop > grp_stop:
            grp_stop = stop
        # append this exon to the current group
        current_grp.append( (start, stop) )
    
    # if exons is not empty and not exons were added then add the last group
    if current_grp:
        exon_groups.append( current_grp )
    
    return dict( zip( xrange(len(exon_groups)), exon_groups ) )

def build_edges_hash_table( overlapping_exons, jns ):
    """ build a hash table that gives graph edges between different clusters.
    this is just the mapping between jn starts, and ends, and allows us
    to quickly follow paths between exons
    """
    # build a mutable set of nodes. Basically, we want to map every exon edge
    # to a node id.
    nodes = {}
    for grp_num, exons in overlapping_exons.iteritems():
        for start, stop in exons:
            assert start not in nodes or nodes[ start ] == grp_num
            assert stop not in nodes or nodes[ stop ] == grp_num
            nodes[ start ] = grp_num
            nodes[ stop ] = grp_num
    
    edges = defaultdict( set )
    for start, stop in jns:
        # do the subtractions so we get exon start/ends, instead of intron bnds
        if stop+1 in nodes: edges[ start-1 ].add( nodes[ stop+1 ] )
        if start-1 in nodes: edges[ stop+1 ].add( nodes[ start-1 ] )
    
    return edges

def build_genelets_from_exons_and_edges( overlapping_exons, edges ):
    def get_connected_node_ids_from_exon_bnds( exon_bnds ):
        node_ids = set()
        for start, stop in exon_bnds:
            if start in edges: node_ids.update( edges[ start ] )
            if stop in edges: node_ids.update( edges[ stop ] )
        
        return node_ids
    
    
    curr_cluster_id = 0
    genelets = defaultdict( list )
    # nodes that we havn't looked at yet. We remove nodes from here when 
    # we find a connection.
    untouched_nodes = set( range(len(overlapping_exons)) )
    while len( untouched_nodes ) > 0:
        # increment the cluster id
        curr_cluster_id += 1
        # store the connected nodes that we still need to explore
        curr_connected_nodes = set( ( untouched_nodes.pop(), ) )
        while len( curr_connected_nodes ) > 0:
            # pop off a node id, and follow the edges
            curr_node_id = curr_connected_nodes.pop()
            
            # add the exons in this node to the new genelet
            genelets[ curr_cluster_id ].extend(
                overlapping_exons[ curr_node_id ] )
            
            # find the nodes curr_node_id connects to
            new_node_ids = get_connected_node_ids_from_exon_bnds(
                overlapping_exons[ curr_node_id ] )
            
            # filter by those that we've already observed, so we take the
            # intersection of nodes that we *havnt* yet observed
            new_node_ids.intersection_update( untouched_nodes )
            
            # add the new nodes to curr_connected_nodes, so that they can 
            # be processed as part of this genelet
            curr_connected_nodes.update( new_node_ids )
            
            # remove these nodes from the untouched list
            untouched_nodes.difference_update( new_node_ids )
    
    return genelets.values()

def filter_exons( exons, junctions ):
    """Remove potential single exon genes
    """
    if len( junctions ) == 0: return []
    
    jn_starts = set( junctions[:,0] )
    jn_stops = set( junctions[:,1] )
    
    filtered_exons = []
    for start, stop in exons:
        # the junctions are inclusive intron bounds, so we subtract one
        # from exon starts to give us the intron stop, and add 1 to the 
        # exon stop to give intron starts.
        if stop+1 in jn_starts or start-1 in jn_stops:
            filtered_exons.append( (start, stop) )
    
    filtered_exons = numpy.array( filtered_exons )
    return filtered_exons

def cluster_exons( exons, junctions ):
    """filter exons such that they must have at least one junction
    then group them by overlapping exons and connected junctions
    """
    filtered_exons = filter_exons( exons, junctions )
    
    if len(filtered_exons) > 0:
        # group overlapping exons
        overlapping_exons = cluster_overlapping_exons( filtered_exons )
        
        # find the hash table of exons ( nodes ) and junctions ( edges )
        edges = build_edges_hash_table( overlapping_exons, junctions )
        
        # build the genelets
        genelets = build_genelets_from_exons_and_edges(
            overlapping_exons, edges )
        
        return genelets
    
    return {}

def cluster_all_exons( exons, jns ):
    all_clustered_exons = {}
    keys = set( exons ).intersection( jns )
    for key in sorted( keys ):
        # Note that exons and jns are sorted numpy.arrays
        all_clustered_exons[key] = cluster_exons(
            numpy.array(exons[key]),
            numpy.array(jns[key]) )
    
    return all_clustered_exons

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(description='Determine valid transcripts.')
    parser.add_argument(
        '--exons-gff', type=file,
        help='GFF of exons - strand will be ignored.')
    parser.add_argument(
        '--junctions-gff', type=file,
        help='A gff file containing valid introns.')
    parser.add_argument(
        '--transcripts-gtf', type=str,
        help='An annotation file from which to extract exons and jns.')
    
    parser.add_argument(
        '--out_fname', '-o', default='',
        help='Output file name. (default: stdout)')
    parser.add_argument(
        '--verbose', '-v', default=False, action= 'store_true',
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    if args.transcripts_gtf != None and \
            (args.junctions_gff != None or args.exons_gff != None ):
        raise ValueError, "Must include either an annotation file, " \
            + "or both an exon and jn file."
    
    out_fp = open( args.out_fname, "w" ) if args.out_fname else sys.stdout
    
    return args.transcripts_gtf, args.exons_gff, args.junctions_gff, out_fp

def main():
    ann_fname, exons_fp, jns_fp, out_fp = parse_arguments()
    if VERBOSE:
        print >> sys.stderr, 'Parsing input...'
    
    jns = defaultdict( set )
    exons = defaultdict( set )
    transcripts_mapping = {}
    if ann_fname != None:
        # load the gtf file
        genes = load_gtf( ann_fname )
        
        for gene in genes:
            # only extract the junctions and exons.
            tss_exons, introns, tes_exons, in_exons = get_elements_from_gene(
                gene, get_tss=False, get_jns=True,
                get_tes=False, get_exons=True )
            
            exons[(gene[1], gene[2])].update( \
                (reg.start, reg.stop) for reg in in_exons )
            jns[(gene[1], gene[2])].update( \
                (reg.start, reg.stop) for reg in introns )
    else:
        exons = parse_exons_file( exons_fp )
        exons_fp.close()
        
        jns = parse_junctions_file( jns_fp )
        jns_fp.close()
    
    if VERBOSE:
        print >> sys.stderr, 'Clustering exons...'
    all_genelets = cluster_all_exons( exons, jns )
    
    if VERBOSE:
        print >> sys.stderr, 'Writing exons...'
    # if we were given an annotation, write it out again
    # with the correct clusters
    if ann_fname != None:
        # build a mapping from exon to cluster id, so that we can assocaite
        # a transcript with the proper exon
        c_id = 1
        exon_cluster_id_mapping = {}
        for (chrm, strand), grouped_exons in all_genelets.iteritems():
            for exons in grouped_exons:
                for exon in exons:
                    exon_cluster_id_mapping[ exon ] = c_id
                    c_id += 1
        
        # now, loop through each trancript and find it's cluster id. We only 
        # need to look at the first exon
        for gene in genes:
            transcripts = gene[-1]
            for trans_id, trans in transcripts:
                gene_id = exon_cluster_id_mapping[tuple(trans[:2])]
                exons = zip( trans[::2], trans[1::2] )
                for exon in exons:
                    assert gene_id == exon_cluster_id_mapping[tuple(exon)]
    
    else:
        # if we weren't given a gtf, then just write out the 
        # clustered exons
        write_all_clusters( all_genelets, out_fp )
    
    return

if __name__ == '__main__':
    main()
