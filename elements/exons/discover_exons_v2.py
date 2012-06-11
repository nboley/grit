import sys, os
from collections import namedtuple
from itertools import izip
import numpy

sys.path.append(os.path.join(os.path.dirname(__file__), "../../", 'file_types'))
from wiggle import Wiggle
from junctions_file import parse_junctions_file_dont_freeze

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

def check_exon_for_gene_split( rd1_cov, rd2_cov ):
    pass

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

def find_genes_in_contig( strand, 
                          merged_read_cov, zero_intervals, jns,
                          read1_cov, read2_cov, 
                          all_cage_cov, polya_cov ):
    # segment the genome
    labels, bndrys = find_initial_boundaries_and_labels( 
        merged_read_cov, zero_intervals, jns )
    
    if strand == '-':
        bndrys = [ len(merged_read_cov)-bndry for bndry in reversed(bndrys)]
        labels = labels[::-1]
        merged_read_cov = merged_read_cov[::-1]
        read1_cov = read1_cov[::-1]
        read2_cov = read2_cov[::-1]
        all_cage_cov = all_cage_cov[::-1]
    
    #labels = refine_se_gene_labels( 
    #    labels, bndrys, read1_cov, read2_cov, all_cage_cov )
    
    # find gene split exons
    for label, (left_bndry, right_bndry) in izip( 
            labels[:-1], izip(bndrys[:-1], bndrys[1:]) ):
        #if label != 'Exon': continue
        if strand != '+': continue
        
        if right_bndry - left_bndry < 2*BIN_BOUNDRY_SIZE + EXON_SPLIT_SIZE + 10:
            continue
        
        PLT_OFF = 0
        cov_1 = read1_cov[ left_bndry-PLT_OFF:right_bndry-1+PLT_OFF ]
        cov_2 = read2_cov[ left_bndry-PLT_OFF:right_bndry-1+PLT_OFF ]
        cage_cov = all_cage_cov[left_bndry-PLT_OFF:right_bndry-1+PLT_OFF]
        merged_cov = merged_read_cov[ left_bndry-PLT_OFF:right_bndry-1+PLT_OFF ]
        merged_cov_cumsum = merged_cov.cumsum()
        
        window_covs = merged_cov_cumsum[EXON_SPLIT_SIZE:] \
            - merged_cov_cumsum[:-EXON_SPLIT_SIZE]
        left_cov = window_covs[0]
        left_bndry = BIN_BOUNDRY_SIZE
        right_cov = window_covs[-1]
        right_bndry = len(merged_cov)-BIN_BOUNDRY_SIZE
        global_min_cov = window_covs[BIN_BOUNDRY_SIZE:-BIN_BOUNDRY_SIZE].min()+1
        
        print strand, label, ( left_bndry, right_bndry ), \
            left_cov, global_min_cov, right_cov
        
        if right_bndry - left_bndry < 1:
            continue
            #return [ start, ], [ 'Exon', ]

        min_index = numpy.argmin( \
            window_covs[left_bndry:right_bndry+1] ) + left_bndry
        # add one to guard against divide zero, and remove low read cov
        min_cov = window_covs[ min_index ] + 1            
        
        # if this does not look like a gene merge exon...
        if left_cov/min_cov < EXON_SPLIT_RATIO \
                or right_cov/min_cov < EXON_SPLIT_RATIO:
            continue

        # find the cage peak
        peaks = find_peaks( cage_cov, 10, 100, 1.0 )
        if len( peaks ) == 0:
            continue
        
        print peaks
        
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

        import matplotlib.pyplot as plt
        plt.plot( range(len(cage_cov)), cage_cov, label="Cage Cov" )
        plt.plot( range(len(cov_1)), cov_1, label="Read 1 Cov" )
        plt.plot( range(len(cov_2)), cov_2, label="Read 2 Cov" )
        plt.axvline( PLT_OFF, color='b' )
        plt.axvline( len(cov_1)-PLT_OFF, color='b' )
        plt.legend()
        plt.show()

        #raw_input()
        
    pass

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
        find_genes_in_contig( key[1],
            merged_read_cov[key], merged_read_cov.zero_intervals[key], jns[key],
            rd1_cov[key], rd2_cov[key],
            cage_cov[key], None )
        
    return

if __name__ == "__main__":
    main()
