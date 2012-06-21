# Copyright (c) 2011-2012 Nathan Boley

import sys, os
from collections import namedtuple, defaultdict
from itertools import izip
import numpy

sys.path.append(os.path.join(os.path.dirname(__file__), "../../", 'file_types'))
from wiggle import Wiggle
from junctions_file import parse_junctions_file_dont_freeze
from build_genelets import cluster_exons

# length of contiguous space to define as an empty region
EMPTY_REGION_SPLIT_SIZE = 80

## Single Exon gene tuning paramaters
MIN_SE_GENE_SIZE = 300
MIN_SE_GENE_CVG = 50
SE_GENE_FRAG_SIZE = 200
SE_GENE_PAIR_MIN_RD_LOG_RATIO = 2.5

EXON_SPLIT_RATIO = 4.0
EXON_SPLIT_SIZE = 40
BIN_BOUNDRY_SIZE = 20

def find_initial_label( prev_bt, bt, next_bt ):
    # if this is the start of an empty region this is always
    # an empty region
    if bt == 'empty':
        return 'empty'
    
    # if this region is surrounded by empty regions, then it is an exon
    # TODO - give this a single exon label?
    if prev_bt == next_bt == 'empty':
        return 'Single'
    
    # we know the next regions is not empty ( because we would have returned 
    # at the line above )
    if bt == 'after_empty':
        assert next_bt != 'empty'
        if next_bt == 'intron': 
            return 'Exon'
        # otherwise, if the next bndry is an exon, this is an extension
        else:
            return 'exon_extension'
    
    # If this is a proper exon i.e. this boundary is an exon_start and
    # the next boundary is an intron_start or empty region
    if bt == 'exon' and next_bt in ( 'intron', 'empty' ):
        return 'Exon'
    
    # if this is a proper intron i.e. this boundary is an intron_start 
    # and the next boundary is an exon_start or an empty region
    if bt == 'intron' and next_bt in ('exon', 'empty' ):
        return 'Intron'

    # otherwise, this is a possible exon extension, but should not be an
    # exon on its own
    # print >> sys.stderr, prev_bt, bt, next_bt
    return 'exon_extension'

def is_in_empty_region( pos, zero_intervals ):
    # otherwise, we do a searchsorted into the zero_intervals array
    index = zero_intervals[:,0].searchsorted( pos, side='right' )
    
    # if the pos is before the first zero_interval or if the 
    # position does not fall within a zero region
    if index == 0 or pos > zero_intervals[index-1][1]:
        return False
    else:
        return True

def find_initial_boundaries_and_labels( read_cov, empty_regions, jns ):
    boundary_types = {}
    # add junctions to boundary types
    for start, stop in jns:
        # move the stop so that boundaries correspond to region starts
        # ( ie, the start is an intron starty and the stop+1 is an 
        # exon start )
        stop += 1

        # skip junctions in empty regions
        if is_in_empty_region( start-1, empty_regions ) \
                or is_in_empty_region( stop, empty_regions ):
            continue

        # add the start 
        boundary_types[ start ] = "intron"
        boundary_types[ stop ] = "exon"

    for start, stop in empty_regions:
        boundary_types[ start ] = "empty"
        boundary_types[ stop+1 ] = "after_empty"

    # build the full set of boundaries and label them
    boundaries = sorted( boundary_types )
    labels = []
    
    # deal with the first label
    if boundary_types[ boundaries[0] ] == 'empty':
        labels.append( 'empty' )
    elif boundary_types[ boundaries[1] ] == 'empty':
        labels.append( 'Exon' )
    else:
        labels.append( 'Exon' )
    
    # build the initial labels set
    for prev_bndry, bndry, next_bndry in izip( \
                boundaries[:-2], boundaries[1:-1], boundaries[2:] ):
        prev_bt = boundary_types[ prev_bndry ]
        bt = boundary_types[ bndry ]
        next_bt = boundary_types[ next_bndry ]

        labels.append( find_initial_label( prev_bt, bt, next_bt ) )
    
    # deal with the last label
    if boundary_types[ boundaries[-1] ] == 'empty':
        labels.append( 'empty' )
    elif boundary_types[ boundaries[-2] ] == 'empty':
        labels.append( 'Exon' )
    else:
        labels.append( 'exon_extension' )
    
    assert len( labels )  == len( boundaries )
    
    return labels, boundaries

def find_peaks( cage_cov, window_len, min_score, max_score_frac ):
    cumsum_cvg_array = \
        numpy.append(0, numpy.cumsum( cage_cov ))
    scores = cumsum_cvg_array[window_len:] - cumsum_cvg_array[:-window_len]
    indices = numpy.argsort( scores )
    
    def overlaps_prev_peak( new_loc ):
        for start, stop in peaks:
            if not( new_loc > stop or new_loc + window_len < start ):
                return True
        return False
    
    # merge the peaks
    def grow_peak( start, stop, grow_size=max(3, window_len/4), min_grow_ratio=0.5 ):
        while True:
            curr_signal = cage_cov[start:stop+1].sum()
            downstream_sig = cage_cov[max(0, start-grow_size):start].sum()
            upstream_sig = cage_cov[stop+1:stop+1+grow_size].sum()

            exp_factor = float(grow_size)/window_len

            # if neither passes the threshold, then return the current peak
            if float(max( upstream_sig, downstream_sig ))*exp_factor \
                    < min_grow_ratio: return (start, stop)
            
            # otherwise, we know one does
            if upstream_sig > downstream_sig:
                stop += grow_size
            else:
                start = max(0, start - grow_size )
            
    
    peaks = []
    peak_scores = []
    
    for index in reversed(indices):
        if not overlaps_prev_peak( index ):
            score = scores[ index ]
            new_peak = grow_peak( index, index + window_len )
            # if we are below the minimum score, then we are done
            if score < min_score:
                return peaks

            # if we have observed peaks, 
            if len( peak_scores ) > 0:
                if float(score)/peak_scores[0] < max_score_frac:
                    return peaks
                        
            peaks.append( new_peak ) 
            peak_scores.append( score )
    
    print >> sys.stderr, "EXHAUSTED EVERY REGION?!?!?!", scores
    return peaks

def refine_se_gene_labels( labels, bndrys, read1_cov, read2_cov, cage_cov ):
    new_labels = []
    for label, (left_bndry, right_bndry) in izip( 
            labels[:-1], izip(bndrys[:-1], bndrys[1:]) ):
        if label != 'Single':
            new_labels.append( label )
            continue
        
        
        cov_1 = read1_cov[ left_bndry:right_bndry-1 ]
        cov_2 = read2_cov[ left_bndry:right_bndry-1 ]
        merged_cov = cov_1 + cov_2
        
        # make sure it is long enough
        assert len( cov_1 ) == len( cov_2 )
        if len( cov_1 ) < MIN_SE_GENE_SIZE:
            new_labels.append( 'empty' )
            continue
        
        # make sure the average coverage is high enough
        if merged_cov.sum()/len( merged_cov ) < MIN_SE_GENE_CVG:
            new_labels.append( 'empty' )
            continue
        
        # make sure the maximum is in the center of the gene
        if len( merged_cov ) - 2*SE_GENE_FRAG_SIZE < 200:
            bndry_size = 100
        else: bndry_size = SE_GENE_FRAG_SIZE
        max_pos = merged_cov.argmax()
        if max_pos < bndry_size or max_pos > len( merged_cov ) - bndry_size:
            new_labels.append( 'empty' )
            continue
        
        # check paired fragment distribution
        # if the ratio of read 1's to rd 2's in the start
        # *or* the ratio of rd 2's to rd 1's in the stop
        # is too low, then this is probably noise
        start_scores = numpy.log((cov_2+1)/(cov_1+1))
        if start_scores[:bndry_size].min() > -SE_GENE_PAIR_MIN_RD_LOG_RATIO:
            new_labels.append( 'empty' )
            continue

        if start_scores[-bndry_size:].max() < SE_GENE_PAIR_MIN_RD_LOG_RATIO:
            new_labels.append( 'empty' )
            continue
        
        new_labels.append( 'Single' )
    
    return new_labels

def split_gene_merge_exons( labels, bndrys, merged_read_cov, 
                           read1_cov, read2_cov, all_cage_cov ):
    def check_gene_merge( region_start, merged_cov, cov_1, cov_2, cage_cov ):
        if len(cov_1) < 2*BIN_BOUNDRY_SIZE + EXON_SPLIT_SIZE + 1:
            return [ region_start, ], [ 'Exon', ]

        merged_cov_cumsum = merged_cov.cumsum()
        
        window_covs = merged_cov_cumsum[EXON_SPLIT_SIZE:] \
            - merged_cov_cumsum[:-EXON_SPLIT_SIZE]
        left_cov = window_covs[0]
        left_bndry = BIN_BOUNDRY_SIZE
        right_cov = window_covs[-1]
        right_bndry = len(merged_cov)-BIN_BOUNDRY_SIZE
        global_min_cov = window_covs[BIN_BOUNDRY_SIZE:-BIN_BOUNDRY_SIZE].min()+1
        
        if right_bndry - left_bndry < 1:
            return [ region_start, ], [ 'Exon', ]

        min_index = numpy.argmin( \
            window_covs[left_bndry:right_bndry+1] ) + left_bndry
        # add one to guard against divide zero, and remove low read cov
        min_cov = window_covs[ min_index ] + 1            
        
        # if this does not look like a gene merge exon...
        if left_cov/min_cov < EXON_SPLIT_RATIO \
                or right_cov/min_cov < EXON_SPLIT_RATIO:
            return [ region_start, ], [ 'Exon', ]
        
        # find the cage peak
        peaks = find_peaks( cage_cov, 10, 100, 1.0 )
        if len( peaks ) == 0:
            return [ region_start, ], [ 'Exon', ]
        
        peak = peaks[0]
        
        """
        bndry_size = 100
        start_scores = numpy.log((cov_2+1)/(cov_1+1))
        if len( start_scores ) < 3*bndry_size: continue
        
        if start_scores[:bndry_size].min() > -SE_GENE_PAIR_MIN_RD_LOG_RATIO:
            continue

        if start_scores[-bndry_size:].max() < SE_GENE_PAIR_MIN_RD_LOG_RATIO:
            continue
        """
        
        #print (cov_1 + cov_2)
        #print numpy.log(((cov_2+1)/(cov_1+1)))

        """
        import matplotlib.pyplot as plt
        plt.plot( range(len(cage_cov)), cage_cov, label="Cage Cov" )
        plt.plot( range(len(cov_1)), cov_1, label="Read 1 Cov" )
        plt.plot( range(len(cov_2)), cov_2, label="Read 2 Cov" )
        plt.legend()
        plt.show()
        """
        
        return [region_start, region_start+peak[0]-10, region_start+peak[0] ], \
            ['exon_extension', 'empty', 'exon_extension']
            
    
    new_labels = []
    new_bndrys = []

    for label, (left_bndry, right_bndry) in izip( 
            labels[:-1], izip(bndrys[:-1], bndrys[1:]) ):
        if label != 'Exon': 
            new_labels.append( label )
            new_bndrys.append( left_bndry )
            continue
        
        cov_1 = read1_cov[ left_bndry:right_bndry-1 ]
        cov_2 = read2_cov[ left_bndry:right_bndry-1 ]
        cage_cov = all_cage_cov[left_bndry:right_bndry-1]
        merged_cov = merged_read_cov[ left_bndry:right_bndry-1 ]

        gs_bndrys, gs_labels = check_gene_merge( 
            left_bndry, merged_cov, cov_1, cov_2, cage_cov  )
        
        new_bndrys.extend( gs_bndrys )
        new_labels.extend( gs_labels )
    
    new_labels.append( labels[-1] )
    new_bndrys.append( bndrys[-1] )
    
    return new_labels, new_bndrys
    
def refine_retained_introns( labels, bndrys, merged_read_cov ):
    new_labels = []
    for label, (left_bndry, right_bndry) in izip( 
            labels[:-1], izip(bndrys[:-1], bndrys[1:]) ):
        if label != 'Intron':
            new_labels.append( label )
            continue
        
        new_labels.append( 'Intron' )
    
    return new_labels

def find_genes_in_contig( strand, 
                          merged_read_cov, zero_intervals,
                          jns_w_cnts,
                          read1_cov, read2_cov, 
                          all_cage_cov, polya_cov ):
    contig_stop = len( merged_read_cov )
    
    # segment the genome
    labels, bndrys = find_initial_boundaries_and_labels( 
        merged_read_cov, zero_intervals, jns_w_cnts.keys() )
    
    if strand == '-':
        return []
        bndrys = [ len(merged_read_cov)-bndry for bndry in reversed(bndrys)]
        labels = labels[::-1]
        merged_read_cov = merged_read_cov[::-1]
        read1_cov = read1_cov[::-1]
        read2_cov = read2_cov[::-1]
        all_cage_cov = all_cage_cov[::-1]
        
    labels = refine_se_gene_labels( 
        labels, bndrys, read1_cov, read2_cov, all_cage_cov )
    
    labels, bndrys = split_gene_merge_exons( 
        labels, bndrys, merged_read_cov, read1_cov, read2_cov, all_cage_cov )
       
    labels = refine_retained_introns( labels, bndrys, merged_read_cov )
    
    clustered_labels, clustered_bndrys = cluster_labels_and_bndrys( 
        labels, bndrys, jns_w_cnts, contig_stop )
    
    return zip( clustered_labels, clustered_bndrys )

###############################################################################
#
#
#  Post initial labelling
#
#


def get_possible_exon_bndrys( labels, boundaries, chrm_stop ):
    """Get all of the possible exon boundaries corresponding to contiguous 'Ee's
    Invalid exons (starting with intron_start or ending with intron_stop) are 
    removed later.
    """
    # get all labels that are either 'Ee'
    exon_indices = []
    for index, label in enumerate( labels ):
        if label in ('TSS', 'TES', 'T_S', 'Exon', 'exon_extension'):
            exon_indices.append( index )
    
    if len( exon_indices ) == 0:
        return []
    
    # find and group contiguous 'Ee' regions
    exon_grps = [ [exon_indices[0],], ]
    for index in exon_indices[1:]:
        if index - exon_grps[-1][-1] == 1:
            exon_grps[-1].append( index )
        else:
            exon_grps.append( [ index,] )
    
    bndrys = []
    for exon_grp in exon_grps:
        # look at all combinations of E or e's in a row
        for start_i, start_exon in enumerate(exon_grp):
            for stop_i, stop_exon in enumerate(exon_grp[start_i:]):
                stop_i += start_i
                # skip grps that contain zero 'E''s
                if not any( labels[index] in ('Exon', 'TSS', 'T_S', 'TES') \
                                for index in exon_grp[start_i:stop_i+1] ):
                    continue
                # if this is the last region goes to the end of the chrm then
                # the last exon goes up to the end of the chromosome
                if stop_exon+1 == len(boundaries):
                    bndrys.append( ( boundaries[start_exon], chrm_stop ) )
                else:
                    bndrys.append( ( boundaries[start_exon], \
                                         boundaries[stop_exon+1]-1 ) )
    
    return bndrys

def build_genelets( labels, bndrys, jns_dict, chrm_stop ):
    # first build a mapping of label indices to other label indices.
    # a label connects to another label if:
    #  1) label_bndry_start - 1 is the stop of another intron
    #  2) label_bndry_stop is an intron start
    #  3) the next or previous label is an 'exon' or 'exon_extension'
    #
    # So, first make these mappings
    #
    
    intron_starts_to_stops = defaultdict(list)
    intron_stops_to_starts = defaultdict(list)
    for start, stop in jns_dict.keys():
        intron_starts_to_stops[start].append( stop )
        intron_stops_to_starts[stop].append( start )

    label_start_to_index = dict( (x, i) for i, x in enumerate(bndrys) )
    label_stop_to_index = dict( (x-1, i) for i, x in enumerate(bndrys[1:]) )
    label_stop_to_index[ chrm_stop ] = len( bndrys )

    # next, loop through each label
    # if it is an exon or exon extension, then find the places that it connects
    # to. And follow those. 
    def find_connected_region_indices( label_index ):
        connected_indices = []
        
        # find the region indices spliced to this
        region_start = bndrys[ label_index ]        
        intron_starts = intron_stops_to_starts[region_start-1]
        #print region_start, intron_starts, \
        #    [x in bndrys for x in intron_starts]
        tmp_conn_indices = [ label_stop_to_index[x-1] 
                             for x in intron_starts 
                             if x-1 in label_stop_to_index ]
        assert all( labels[i] != 'empty' for i in tmp_conn_indices )
        connected_indices.extend( tmp_conn_indices )
        region_stop = bndrys[ label_index+1 ] \
            if label_index+1 < len(bndrys) else chrm_stop
        
        # find the region indices spliced from this
        spliced_region_starts = intron_starts_to_stops[region_stop]
        tmp_conn_indices = [ label_start_to_index[x+1] 
                             for x in spliced_region_starts ]
        assert all( labels[i] != 'empty' for i in tmp_conn_indices )
        connected_indices.extend( tmp_conn_indices )
        
        if label_index+1 < len(labels) \
                and labels[label_index+1] in ('exon_extension', 'Exon' ):
            connected_indices.append( label_index+1 )
        if label_index > 0 \
                and labels[label_index-1] in ('exon_extension', 'Exon' ):
            connected_indices.append( label_index-1 )
        
        return connected_indices
    
    def find_all_connected_region_indices( start_index ):
        observed_indices = set((start_index,))
        new_indices = [start_index,]
        while len( new_indices ) > 0:
            index = new_indices.pop()
            connected_indices = find_connected_region_indices( index )
            new_indices.extend( i for i in connected_indices 
                                if i not in observed_indices  )
            observed_indices.update( connected_indices )
    
        return observed_indices
    
    unobserved_label_indices = set()
    for i, label in enumerate(labels):
        if label in ('empty', 'Intron'): continue
        unobserved_label_indices.add( i ) 
    
    clustered_indices = []
    while len( unobserved_label_indices ) > 0:
        start_index = unobserved_label_indices.pop()
        print "START INDEX:", start_index, labels[start_index]
        
        tmp_conn_indices = find_all_connected_region_indices(start_index)
        
        # skip single exon genes
        if len(tmp_conn_indices ) == 1: 
            continue
        
        clustered_indices.append( tmp_conn_indices )
        
        assert len( clustered_indices[-1] ) == len(set(clustered_indices[-1]))
        for index in clustered_indices[-1]:
            print index, bndrys[index], labels[index]
            if index != start_index:
                assert index in unobserved_label_indices
                unobserved_label_indices.remove( index )
        
    for cluster in clustered_indices:
        print cluster
    
    assert False
"""
def build_genelets( labels, bndrys, jns_dict, chrm_stop ):
    # get exon boundaries from assigned labels and boundaries
    exon_bndrys = get_possible_exon_bndrys( \
        labels, bndrys, chrm_stop )
    
    jns_wo_cnts = jns_dict.keys()
    
    # filter out any exon that start with a junction start or 
    # stops at an junction stop
    filtered_exon_bndrys = filter_exon_bndrys( \
        exon_bndrys, jns_wo_cnts )
    filtered_exon_bndrys = exon_bndrys
    
    bndries_array = numpy.array( bndrys )
    exons = numpy.array( filtered_exon_bndrys )
    clusters = cluster_exons( exons, jns_wo_cnts )
    
    clustered_labels = []
    clustered_bndrys = []
    for cluster in clusters:
        print cluster
        raw_input()
        
        cluster_bndry_indices = set()
        for exon in cluster:
            start_index, stop_index = find_bndry_indices( bndries_array, exon )
            cluster_bndry_indices.update( xrange( start_index, stop_index+1 ) )
        
        cluster_bndry_indices = sorted( cluster_bndry_indices )
        
        cluster_bndrys, cluster_labels = build_cluster_labels_and_bndrys( 
            cluster_bndry_indices, bndries_array, labels, chrm_stop )
        
        clustered_labels.append( cluster_labels )
        clustered_bndrys.append( cluster_bndrys )
    
    return clustered_labels, clustered_bndrys
"""

def cluster_labels_and_bndrys( labels, bndrys, jns_w_cnts, chrm_stop ):
    jns_dict = dict( jns_w_cnts )
    
    grpd_clusters = [ ( labels, bndrys), ]
    final_cl_labels, final_cl_bndrys = [], []
    while len( grpd_clusters ) > 0:
        labels, bndrys = grpd_clusters.pop()
                
        clustered_labels, clustered_bndrys = build_genelets( \
            labels, bndrys, jns_dict, chrm_stop )
        
        # if we didn't re cluster these, then we are done
        assert len( clustered_labels ) == len( clustered_bndrys )
        if len( clustered_labels ) == 1:             
            final_cl_labels.append( clustered_labels[0] )
            final_cl_bndrys.append( clustered_bndrys[0] )
        else:
            grpd_clusters.extend( zip(*(clustered_labels, clustered_bndrys)) )
    
    return final_cl_labels, final_cl_bndrys


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(\
        description='Find exons from wiggle and junctions files.')

    parser.add_argument( 'junctions', type=file, \
        help='GTF format file of junctions(introns).')
    parser.add_argument( 'chrm_sizes_fname', type=file, \
        help='File with chromosome names and sizes.')

    parser.add_argument( 'wigs', type=file, nargs="+", \
        help='wig files over which to search for exons.')
    
    parser.add_argument( '--cage-wigs', type=file, nargs='+', \
        help='wig files with cage reads, to identify tss exons.')
    parser.add_argument( '--polya-reads-gffs', type=file, nargs='+', \
        help='files with polya reads, to identify tes exons.')
    
    parser.add_argument( '--out-file-prefix', '-o', default="discovered_exons",\
        help='Output file name. (default: discovered_exons)')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
        help='Whether or not to print status information.')

    args = parser.parse_args()
    
    OutFPS = namedtuple( "OutFPS", [\
            "single_exon_genes", "tss_exons", \
            "internal_exons", "tes_exons", "all_exons"])
    fps = []
    for field_name in OutFPS._fields:
        fps.append(open("%s.%s.gff" % (args.out_file_prefix, field_name), "w"))
    ofps = OutFPS( *fps )
        
    # set flag args
    global VERBOSE
    VERBOSE = args.verbose
    
    rd1_plus_wigs = [ fp for fp in args.wigs 
                      if fp.name.endswith("rd1.plus.bedGraph") ]
    rd1_minus_wigs = [ fp for fp in args.wigs 
                      if fp.name.endswith("rd1.minus.bedGraph") ]
    rd2_plus_wigs = [ fp for fp in args.wigs 
                      if fp.name.endswith("rd2.plus.bedGraph") ]
    rd2_minus_wigs = [ fp for fp in args.wigs 
                      if fp.name.endswith("rd2.minus.bedGraph") ]
    
    grpd_wigs = [ rd1_plus_wigs, rd1_minus_wigs, rd2_plus_wigs, rd2_minus_wigs ]
    
    return grpd_wigs, args.junctions, args.chrm_sizes_fname, \
        args.cage_wigs, args.polya_reads_gffs, ofps

def main():
    wigs, jns_fp, chrm_sizes_fp, cage_wigs, polya_reads_gffs, ofps \
        = parse_arguments()
    
    cage_cov = Wiggle( chrm_sizes_fp, cage_wigs )

    jns = parse_junctions_file_dont_freeze( jns_fp )
    
    rd1_cov = Wiggle( 
        chrm_sizes_fp, [wigs[0][0], wigs[1][0]], ['+','-'] )
    rd2_cov = Wiggle( 
        chrm_sizes_fp, [wigs[2][0], wigs[3][0]], ['+','-'] )
    
    merged_read_cov = Wiggle(
        chrm_sizes_fp,
        [ wigs[0][0], wigs[1][0], wigs[2][0], wigs[3][0] ],
        ['+', '-', '+', '-']
    )
    merged_read_cov.calc_zero_intervals()
    
    for key in merged_read_cov:
        introns = jns[key]
        scores = jns._scores[key]
        jns_and_scores = dict(izip( introns, scores ))


        clusters = find_genes_in_contig( key[1], 
            merged_read_cov[key], merged_read_cov.zero_intervals[key], 
            jns_and_scores,
            rd1_cov[key], rd2_cov[key],
            cage_cov[key], None )
        
        print clusters
        
        for cluster in clusters:
            print cluster
        
    return

if __name__ == "__main__":
    main()
