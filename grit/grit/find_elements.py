import sys
import numpy
from scipy import stats
from collections import defaultdict
from itertools import chain, izip
from bisect import bisect
from copy import copy

from igraph import Graph

from files.reads import RNAseqReads, CAGEReads, RAMPAGEReads, clean_chr_name, guess_strand_from_fname, iter_coverage_intervals_for_read
from files.junctions import extract_junctions_in_contig
from files.bed import create_bed_line

USE_CACHE = False

def flatten( regions ):
    new_regions = []
    curr_start = regions[0][0]
    curr_end = regions[0][1]
    for i,(start,end) in enumerate(regions):
        if curr_end > end:
            end = curr_end
        if i+1 == len( regions ): 
            if len(new_regions) == 0:
                new_regions = [curr_start, curr_end]
                break
            if new_regions[-1][1] == end:
                break
            else:
                new_regions.append( [ curr_start, curr_end] )
                break
        if regions[i+1][0]-end <= 1:
            curr_end = max( regions[i+1][1], end ) 
        else:
            new_regions.append( [curr_start, curr_end] )
            curr_start = regions[i+1][0]
            curr_end = regions[i+1][1]
    if type(new_regions[0]) == int:
        return [new_regions]
    else:
        return new_regions

MIN_REGION_LEN = 50
EDGE_THRESHOLD_RATIO = 0.001
EMPTY_BPK = 0
MIN_BPK = 1
MIN_NUM_CAGE_TAGS = 20
MIN_EXON_BPKM = 1.0
MAX_CAGE_FRAC = 0.05
EXON_EXT_CVG_RATIO_THRESH = 4

def get_contigs_and_lens( rnaseq_reads, cage_reads ):
    chrm_lengths = {}
    for bam in chain( rnaseq_reads, cage_reads ):
        for ref_name, ref_len in zip(bam.references, bam.lengths):
            if clean_chr_name( ref_name ) not in chrm_lengths:
                chrm_lengths[clean_chr_name( ref_name )] = ref_len
            else:
                assert chrm_lengths[clean_chr_name(ref_name)] == ref_len, \
                    "Chromosome lengths do not match between bam files"
    
    return chrm_lengths

def build_empty_array():
    return numpy.array(())

def find_empty_regions( cov, thresh=1 ):
    return []
    x = numpy.diff( numpy.asarray( cov >= thresh, dtype=int ) )
    return zip(numpy.nonzero(x==1)[0],numpy.nonzero(x==-1)[0])

def merge_empty_labels( poss ):
    locs = [i[1] for i in poss ]
    labels_to_remove = []
    label_i = 0
    while label_i < len(locs):
        if locs[label_i] != "ESTART": 
            label_i += 1
            continue
        label_n = label_i+1
        while label_n < len(locs) \
                and locs[ label_n ] in ("ESTART", "ESTOP"):
            label_n += 1
        if label_n - label_i > 1:
            labels_to_remove.append( [label_i+1, label_n-1]  )
        label_i = label_n+1
    
    for start, stop in reversed( labels_to_remove ):
        del poss[start:stop+1]

    return poss

def get_qrange_long( np, w_step, w_window ):
    L = len(np)-w_window
    if L < 2:
       nm, nx = get_qrange_short( np )
       return nm, nx
    Q = []
    for pos in xrange(0, L, w_step):
        Q.append( np[pos:pos+w_window].mean() )
    Q = numpy.asarray(Q)
    lower = stats.stats.scoreatpercentile(Q, 10)
    upper = stats.stats.scoreatpercentile(Q, 90)
    return lower, upper

def get_qrange_short( np ):
    L = len(np)
    return np.min(), np.max()

def filter_exon( exon, wig, min_avg_cvg=0.01, 
                 min_short_cvg=1.0, short_exon_length=400,
                 min_long_cvg=10 ):
    '''Find all the exons that are sufficiently homogenous and expressed.
    
    '''
    start = exon.start
    end = exon.stop
    vals = wig[start:end+1]
    mean_cvg = vals.mean()
    
    # if virtually off, forget it
    if mean_cvg < min_avg_cvg: 
        return True
    
    # if its short, just make sure it's covered to at least an average of 1X 
    if end-start < short_exon_length: 
        if mean_cvg > min_short_cvg:
            return False
        else:
            return True
    
    # get the interquartile range
    low, high = get_qrange_long( vals, 50, 300 ) 
    # if the lower quartile is 0, ditch it
    if low == 0 or high < min_long_cvg: 
        return True
    
    IQ_ratio = high / max(low, min_avg_cvg)
    # if the IQ range is < 100, chalk it up to inhomogeneity and accept it
    if IQ_ratio < 100:
        return False
    
    # if the IQ range is "boarder line", but the low is high, keep it, what the hell
    if IQ_ratio >= 100 and IQ_ratio < 500 and low > 50:
        return False
    
    return True

def filter_exons( exons, rnaseq_cov ):
    for exon in exons:
        if not filter_exon( exon, rnaseq_cov ):
            yield exon
    
    return

def find_polya_sites( polya_sites_fnames ):
    locs = defaultdict( list )
    for fname in polya_sites_fnames:
        strand = guess_strand_from_fname( fname )
        with open( fname ) as fp:
            for line in fp:
                if line.startswith( "track" ): continue
                data = line.split()
                chrm, start, stop, value = \
                    data[0], int(data[1]), int(data[2]), float(data[3])
                assert start == stop
                assert value == 1
                locs[(chrm, strand)].append( start )
    
    # convert to a dict of sorted numpy arrays
    numpy_locs = defaultdict( build_empty_array )

    for (chrm, strand), polya_sites in locs.iteritems():
        # make sure they're unique
        assert len( polya_sites ) == len( set( polya_sites ) )

        polya_sites.sort()
        if chrm.startswith( 'chr' ):
            chrm = chrm[3:]
        
        numpy_locs[(chrm, strand)] = numpy.array( polya_sites )
    
    return numpy_locs

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(\
        description='Find exons from RNAseq, CAGE, and poly(A) assays.')

    parser.add_argument( 'rnaseq_reads',type=argparse.FileType('rb'),nargs='+',\
        help='BAM files containing mapped RNAseq reads ( must be indexed ).')
    
    parser.add_argument( '--cage-reads', type=file, default=[], nargs='*', \
        help='BAM files containing mapped cage reads.')
    
    parser.add_argument( '--polya-candidate-sites', type=file, nargs='*', \
        help='files with allowed polya sites.')
    
    parser.add_argument( '--out-filename', '-o', 
                         default="discovered_elements.bed",\
        help='Output file name. (default: discovered_elements.bed)')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
        help='Whether or not to print status information.')
    parser.add_argument('--write-debug-data',default=False,action='store_true',\
        help='Whether or not to print out gff files containing intermediate exon assembly data.')
    parser.add_argument( '--threads', '-t', default=1, type=int,
        help='The number of threads to use.')
        
    args = parser.parse_args()

    global num_threads
    num_threads = args.threads
    
    global WRITE_DEBUG_DATA
    WRITE_DEBUG_DATA = args.write_debug_data
    
    global VERBOSE
    VERBOSE = args.verbose
        
    ofp = open( args.out_filename, "w" )
    
    return args.rnaseq_reads, args.cage_reads, args.polya_candidate_sites, ofp

class Bin( object ):
    def __init__( self, start, stop, left_label, right_label, 
                  bin_type=None, score=1000 ):
        self.start = start
        self.stop = stop
        assert stop - start > 0
        self.left_label = left_label
        self.right_label = right_label
        self.type = bin_type
        self.score = score
    
    def mean_cov( self, cov_array ):
        return cov_array[self.start:self.stop].mean()
    
    def reverse_strand(self, contig_len):
        return Bin(contig_len-self.stop, contig_len-self.start, 
                   self.right_label, self.left_label, self.type)

    def shift(self, shift_amnt):
        return Bin(self.start+shift_amnt, self.stop+shift_amnt, 
                   self.left_label, self.right_label, self.type)
    
    def __repr__( self ):
        if self.type == None:
            return "%i-%i\t%s\t%s" % ( self.start, self.stop, self.left_label,
                                       self.right_label )

        return "%s:%i-%i" % ( self.type, self.start, self.stop )
    
    def __hash__( self ):
        return hash( (self.start, self.stop) )
    
    def __eq__( self, other ):
        return ( self.start == other.start and self.stop == other.stop )
    
    _bndry_color_mapping = {
        'CONTIG_BNDRY': '0,0,0',
        
        'POLYA': '255,255,0',

        'CAGE_PEAK': '0,255,0',
        
        'D_JN': '173,255,47',
        'R_JN': '0,0,255',
        
        'ESTART': '0,0,0',
        'ESTOP': '0,0,0'
    }
    
    def find_bndry_color( self, bndry ):
        return self._bndry_color_mapping[ bndry ]
    
    def _find_colors( self, strand ):
        if self.type != None:
            if self.type =='GENE':
                return '0,0,0'
            if self.type =='CAGE_PEAK':
                return '0,0,0'
            if self.type =='EXON':
                return '0,0,0'
            if self.type =='EXON_EXT':
                return '0,0,255'
            if self.type =='RETAINED_INTRON':
                return '255,255,0'
            if self.type =='TES_EXON':
                return '255,0,0'
            if self.type =='TSS_EXON':
                return '0,255,0'
            if self.type =='SE_GENE':
                return '255,255,0'

        if strand == '+':
            left_label, right_label = self.left_label, self.right_label
        else:
            assert strand == '-'
            left_label, right_label = self.right_label, self.left_label
        
        
        if left_label == 'D_JN' and right_label  == 'R_JN':
            return '108,108,108'
        if left_label == 'D_JN' and right_label  == 'D_JN':
            return '135,206,250'
        if left_label == 'R_JN' and right_label  == 'R_JN':
            return '135,206,250'
        if left_label == 'R_JN' and right_label  == 'D_JN':
            return '0,0,255'
        if left_label == 'R_JN' and right_label  == 'POLYA':
            return '255,0,0'
        if left_label == 'POLYA' and right_label  == 'POLYA':
            return ' 240,128,128'
        if left_label == 'D_JN' and right_label  == 'POLYA':
            return '240,128,128'
        if left_label == 'POLYA' and right_label  == 'D_JN':
            return '147,112,219'
        if left_label == 'POLYA' and right_label  == 'R_JN':
            return '159,153,87'
        if left_label == 'ESTART' and right_label  == 'ESTOP':
            return '159,153,87'
        
        return ( self.find_bndry_color(left_label), 
                 self.find_bndry_color(right_label) )

def write_unified_bed( elements, ofp ):
    assert isinstance( elements, Bins )
    
    feature_mapping = { 
        'GENE': 'gene',
        'CAGE_PEAK': 'promoter',
        'SE_GENE': 'single_exon_gene',
        'TSS_EXON': 'tss_exon',
        'EXON': 'internal_exon',
        'TES_EXON': 'tes_exon',
        'INTRON': 'intron',
        'POLYA': 'polya'
    }

    color_mapping = { 
        'GENE': '200,200,200',
        'CAGE_PEAK': '153,255,000',
        'SE_GENE': '000,000,000',
        'TSS_EXON': '140,195,59',
        'EXON': '000,000,000',
        'TES_EXON': '255,51,255',
        'INTRON': '100,100,100',
        'POLYA': '255,0,0'
    }
        
    for bin in elements:
            region = ( elements.chrm, elements.strand, bin.start, bin.stop)
            grp_id = feature_mapping[bin.type] + "_%s_%s_%i_%i" % region
            bed_line = create_bed_line( elements.chrm, elements.strand, 
                                        bin.start, bin.stop, 
                                        feature_mapping[bin.type],
                                        score=bin.score,
                                        color=color_mapping[bin.type],
                                        use_thick_lines=(bin.type != 'INTRON'))
            ofp.write( bed_line + "\n"  )
    return

class Bins( list ):
    def __init__( self, chrm, strand, iter=[] ):
        self.chrm = chrm
        self.strand = strand
        self.extend( iter )
        self._bed_template = "\t".join( ["chr"+chrm, '{start}', '{stop}', '{name}', 
                                         '1000', strand, '{start}', '{stop}', 
                                         '{color}']  ) + "\n"
        
    def reverse_strand( self, contig_len ):
        rev_bins = Bins( self.chrm, self.strand )
        for bin in reversed(self):
            rev_bins.append( bin.reverse_strand( contig_len ) )
        return rev_bins

    def shift( self, shift_amnt ):
        shifted_bins = Bins( self.chrm, self.strand )
        for bin in self:
            shifted_bins.append( bin.shift( shift_amnt ) )
        return shifted_bins
    
    def writeBed( self, ofp, contig_len, reverse_strand=True ):
        """
            chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
        """
        if reverse_strand and self.strand == '-':
            writetable_bins = self.reverse_strand( contig_len )
        else:
            writetable_bins = self
        
        for bin in writetable_bins:
            length = max( (bin.stop - bin.start)/4, 1)
            colors = bin._find_colors( self.strand )
            if isinstance( colors, str ):
                op = self._bed_template.format(
                    start=bin.start,stop=bin.stop+1,color=colors, 
                    name="%s_%s"%(bin.left_label, bin.right_label) )
                ofp.write( op )
            else:
                op = self._bed_template.format(
                    start=bin.start,stop=(bin.start + length),color=colors[0], 
                    name=bin.left_label)
                ofp.write( op )
                
                op = self._bed_template.format(
                    start=(bin.stop-length),stop=bin.stop,color=colors[1],
                    name=bin.right_label)
                ofp.write( op )
        
        return

    def writeGff( self, ofp, contig_len, filter=None ):
        """
            chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
        """
        if self.strand == '-':
            writetable_bins = self.reverse_strand( contig_len )
        else:
            writetable_bins = self
        
        for bin in writetable_bins:
            if filter != None and bin.type != filter:
                continue
            region = GenomicInterval(self.chrm, self.strand, 
                                     bin.start, bin.stop)
            grp_id = "%s_%s_%i_%i" % region
            ofp.write( create_gff_line(region, grp_id) + "\n" )
        
        return

def load_junctions( rnaseq_reads, chrm, strand ):
    if USE_CACHE:
        import cPickle
        fname = "tmp.junctions.%s.%s.obj" % ( chrm, strand )
        try:
            with open( fname ) as fp:
                junctions = cPickle.load(fp)
        except IOError:
            junctions = extract_junctions_in_contig( 
                rnaseq_reads[0], chrm, strand )
            with open( fname, "w" ) as fp:
                cPickle.dump(junctions, fp)
    else:
        junctions = extract_junctions_in_contig( 
            rnaseq_reads[0], chrm, strand )
    
    if VERBOSE: print "Finished extracting junctions for %s %s" % (chrm, strand)
    return junctions

def find_gene_boundaries((chrm, strand, contig_len), rnaseq_reads, 
                         polya_sites, junctions=None):
    def find_splits_in_region(rnaseq_reads, chrm, strand, start, stop):
        n_splits = int( (stop-start)/MIN_REGION_LEN )
        segment_len = int((stop-start)/(n_splits+1))
        return [ start + i*segment_len for i in xrange(1,n_splits+1) ]
    
    def find_initial_segmentation( chrm, strand, rnaseq_reads, polya_sites ):
        locs = {0: 'SEGMENT', contig_len-1:'SEGMENT'}
        for polya in polya_sites:
            locs[ polya ] = "POLYA"
        
        return locs
    
    def merge_segments( segments ):
        bndries = numpy.array( sorted( segments ) )
        bndries_to_delete = set()
        for i, bndry in enumerate(bndries[:-1]):
            if bndries[i+1]-bndry > 10000:
                continue
            pre_bndry_cnt = 0
            post_bndry_cnt = 0
            
            for reads in rnaseq_reads:
                cvg = reads.build_read_coverage_array( 
                    chrm, strand, max(0,bndry-50), bndry+50 )
                pre_bndry_cnt += cvg[:50].sum()
                post_bndry_cnt += cvg[50:].sum()
        
            if post_bndry_cnt > 50 and \
                (post_bndry_cnt)/(post_bndry_cnt+pre_bndry_cnt+1e-6) > 0.20:
                bndries_to_delete.add( bndry )
        
        for bndry_to_delete in bndries_to_delete:
            del segments[bndry_to_delete]
        
        return segments
    
    def cluster_segments( segments, jns ):
        boundaries = numpy.array( sorted( segments ) )
        edges = set()
        for (start, stop), cnt in jns:
            start_bin = boundaries.searchsorted( start-1 )-1
            stop_bin = boundaries.searchsorted( stop+1 )-1
            if start_bin != stop_bin:
                edges.add((min(start_bin, stop_bin), max(start_bin, stop_bin)))
        
        genes_graph = Graph( len( boundaries )-1 )
        genes_graph.add_edges( list( edges ) )
        
        segments = []
        for g in genes_graph.clusters():
            segments.append( (boundaries[min(g)]+1, boundaries[max(g)+1]-1) )
        
        return flatten( segments )
    
    # find all of the junctions
    if None == junctions:
        junctions = load_junctions( rnaseq_reads, chrm, strand )
    
    # find segment boundaries
    initial_segmentation = find_initial_segmentation( 
        chrm, strand, rnaseq_reads, polya_sites )
    if VERBOSE: print "Finished initial segmentation for %s %s" % (chrm, strand)
    
    merged_segments = merge_segments( initial_segmentation )
    clustered_segments = cluster_segments( merged_segments, junctions )
    
    # build the gene bins, and write them out to the elements file
    genes = Bins( chrm, strand, [] )
    if len( clustered_segments ) == 2:
        return genes
    
    for start, stop in clustered_segments:
        if stop - start < 100: continue
        if stop - start > 10000000: continue
        genes.append( Bin(max(0,start-10), min(stop+10,contig_len), 
                          "ESTART", "POLYA", "GENE" ) )
    
    return genes

def find_cage_peaks_in_gene( ( chrm, strand ), gene, cage_cov, rnaseq_cov ):
     raw_peaks = find_peaks( cage_cov, window_len=20, 
                             min_score=MIN_NUM_CAGE_TAGS,
                             max_score_frac=MAX_CAGE_FRAC, 
                             max_num_peaks=20 )
     if len( raw_peaks ) == 0:
         return []
     
     cage_peaks = Bins( chrm, strand )
     for peak_st, peak_sp in raw_peaks:
         cage_peaks.append( Bin( peak_st, peak_sp+1,
                                 "CAGE_PEAK_START", "CAGE_PEAK_STOP", "CAGE_PEAK") )
     return cage_peaks

def find_peaks( cov, window_len, min_score, max_score_frac, max_num_peaks ):
    cumsum_cvg_array = \
        numpy.append(0, numpy.cumsum( cov ))
    scores = cumsum_cvg_array[window_len:] - cumsum_cvg_array[:-window_len]
    indices = numpy.argsort( scores )
    
    def overlaps_prev_peak( new_loc ):
        for start, stop in peaks:
            if not( new_loc > stop or new_loc + window_len < start ):
                return True
        return False
    
    # merge the peaks
    def grow_peak( start, stop, grow_size=max(3, window_len/4), min_grow_ratio=0.5 ):
        # grow a peak at most max_num_peaks times
        for i in xrange(max_num_peaks):
            curr_signal = cov[start:stop+1].sum()
            if curr_signal < min_score:
                return ( start, stop )
            
            downstream_sig = cov[max(0, start-grow_size):start].sum()
            upstream_sig = cov[stop+1:stop+1+grow_size].sum()
            exp_factor = float( stop - start + 1 )/grow_size
            
            # if neither passes the threshold, then return the current peak
            if float(max( upstream_sig, downstream_sig ))*exp_factor \
                    < curr_signal*min_grow_ratio: return (start, stop)
            
            # otherwise, we know one does
            if upstream_sig > downstream_sig:
                stop += grow_size
            else:
                start = max(0, start - grow_size )
        
        if VERBOSE:
            print "Warning: reached max peak iteration at %i-%i ( signal %.2f )"\
                % (start, stop, cov[start:stop+1].sum() )
            print 
        return (start, stop )
    
    peaks = []
    peak_scores = []
    
    for index in reversed(indices):
        if not overlaps_prev_peak( index ):
            score = scores[ index ]
            new_peak = grow_peak( index, index + window_len )
            # if we are below the minimum score, then we are done
            if score < min_score:
                break
            
            # if we have observed peaks, and the ratio between the highest
            # and the lowest is sufficeintly high, we are done
            if len( peak_scores ) > 0:
                if float(score)/peak_scores[0] < max_score_frac:
                    break
                        
            peaks.append( new_peak ) 
            peak_scores.append( score )
    
    if len( peaks ) == 0:
        return []
    
    # merge cage peaks together
    def merge_peaks( peaks_and_scores ):
        merged_peaks = set()
        new_peaks = []
        new_scores = []
        for pk_i, (peak, score) in enumerate(peaks_and_scores):
            if pk_i in merged_peaks: continue
            curr_pk = list( peak )
            curr_score = score
            for i_pk_i, (i_peak, i_score) in enumerate(peaks_and_scores):
                if i_pk_i in merged_peaks: continue
                if i_peak[0] < curr_pk[0]: continue
                if i_peak[0] - curr_pk[1] < max( window_len, 
                                                 curr_pk[1]-curr_pk[0] ):
                    curr_pk[1] = i_peak[1]
                    curr_score += i_score
                    merged_peaks.add( i_pk_i )
                else:
                    break

            new_peaks.append( curr_pk )
            new_scores.append( curr_score )
        return zip( new_peaks, new_scores )
    
    peaks_and_scores = sorted( zip(peaks, peak_scores) )
    old_len = len( peaks_and_scores )
    for i in xrange( 99 ):
        if i == 100: assert False
        peaks_and_scores = merge_peaks( peaks_and_scores )
        if len( peaks_and_scores ) == old_len: break
    
    max_score = max( s for p, s in peaks_and_scores )
    return [ pk for pk, score in peaks_and_scores \
                 if score/max_score > max_score_frac
                 and score > min_score ]


def find_left_exon_extensions( start_index, start_bin, gene_bins, rnaseq_cov ):
    internal_exons = []
    ee_indices = []
    start_bin_cvg = start_bin.mean_cov( rnaseq_cov )
    for i in xrange( start_index-1, 0, -1 ):
        bin = gene_bins[i]

        # break at canonical exons
        if bin.left_label == 'R_JN' and bin.right_label == 'D_JN':
            break
        
        # make sure the average coverage is high enough
        bin_cvg = bin.mean_cov(rnaseq_cov)
        
        if bin_cvg < MIN_EXON_BPKM:
            break
        
        if bin.stop - bin.start > 20 and \
                start_bin_cvg/(bin_cvg+1) >EXON_EXT_CVG_RATIO_THRESH:
            break

        # update the bin coverage. In cases where the coverage increases from
        # the canonical exon, we know that the increase is due to inhomogeneity
        # so we take the conservative choice
        start_bin_cvg = max( start_bin_cvg, bin_cvg )
        
        ee_indices.append( i )
        internal_exons.append( Bin( bin.start, bin.stop, 
                                    bin.left_label, bin.right_label, 
                                    "EXON_EXT"  ) )
    
    return internal_exons

def find_right_exon_extensions( start_index, start_bin, gene_bins, rnaseq_cov ):
    exons = []
    ee_indices = []
    start_bin_cvg = start_bin.mean_cov( rnaseq_cov )
    for i in xrange( start_index+1, len(gene_bins) ):
        bin = gene_bins[i]

        # if we've reached a canonical exon, break
        if bin.left_label == 'R_JN' and bin.right_label == 'D_JN':
            break
        
        # make sure the average coverage is high enough
        bin_cvg = rnaseq_cov[bin.start:bin.stop].mean()
        if bin_cvg < MIN_EXON_BPKM:
            break
        
        if bin.stop - bin.start > 20 and \
                start_bin_cvg/(bin_cvg+1) >EXON_EXT_CVG_RATIO_THRESH:
            break
        
        # update the bin coverage. In cases where the coverage increases from
        # the canonical exon, we know that the increase is due to inhomogeneity
        # so we take the conservative choice
        start_bin_cvg = max( start_bin_cvg, bin_cvg )

        ee_indices.append( i )
                
        exons.append( Bin( bin.start, bin.stop, 
                                    bin.left_label, bin.right_label, 
                                    "EXON_EXT"  ) )
    
    return exons

def find_internal_exons_from_pseudo_exons( pseudo_exons ):
    """Find all of the possible exons from the set of pseudo exons.
    
    """
    internal_exons = []
    # condition on the number of pseudo exons
    for exon_len in xrange(1, len(pseudo_exons)+1):
        for start in xrange(len(pseudo_exons)-exon_len+1):
            start_ps_exon = pseudo_exons[start]
            stop_ps_exon = pseudo_exons[start+exon_len-1]

            if stop_ps_exon.right_label != 'D_JN': continue

            if start_ps_exon.left_label != 'R_JN': continue

            # each potential exon must have a canonical exon in it
            if not any ( exon.type == 'EXON' 
                         for exon in pseudo_exons[start:start+exon_len] ):
                continue
            internal_exons.append( Bin(start_ps_exon.start, stop_ps_exon.stop,
                                       start_ps_exon.left_label, 
                                       stop_ps_exon.right_label,
                                       'EXON') )
    
    internal_exons = sorted( internal_exons )
    
    return internal_exons

def find_tss_exons_from_pseudo_exons( tss_exon, pseudo_exons ):
    tss_exons = []
    # find the first pseudo exon that the tss exon connects to
    for i, pse in enumerate( pseudo_exons ):
        if pse.start == tss_exon.stop: 
            break
    
    for pse in pseudo_exons[i:]:
        tss_exons.append( Bin(tss_exon.start, pse.stop,
                              tss_exon.left_label, pse.right_label, "TSS_EXON"))
    
    return tss_exons

def find_pseudo_exons_in_gene( ( chrm, strand ), gene, 
                               cage_peaks, polya_sites, jns,
                               rnaseq_cov, cage_cov ):        
    locs = {}
    for polya in polya_sites:
        locs[ polya ] = "POLYA"
    
    for start, stop, cnt in jns:
        locs[start-1] = "D_JN"
        locs[stop+1] = "R_JN"

    for start, stop in find_empty_regions( rnaseq_cov ):
        if stop - start < MIN_EMPTY_REGION_LEN: continue
        if start in locs or stop in locs:
            continue
        locs[start] = "ESTART"
        locs[stop] = "ESTOP"
        
    # build all of the bins
    poss = sorted( locs.iteritems() )
    poss = merge_empty_labels( poss )
    if len( poss ) == 0:
        return [], [], [], []

    gene_bins = []
    for index, ((start, left_label), (stop, right_label)) in \
            enumerate(izip(poss[:-1], poss[1:])):
        gene_bins.append( Bin(start, stop, left_label, right_label) )

        
    # find tss exons
        # overlaps a cage peak, then 
    # find exon starts pseudo exons ( spliced_to, ... )
    pseudo_exons = []
    canonical_exon_starts = [ i for i, bin in enumerate( gene_bins )
                              if bin.left_label == 'R_JN' 
                              and bin.right_label == 'D_JN' ]
    for ce_i in canonical_exon_starts:
        canonical_bin = gene_bins[ce_i]
        if canonical_bin.mean_cov(rnaseq_cov) < MIN_EXON_BPKM:
            continue
        
        pseudo_exons.append( Bin( canonical_bin.start, canonical_bin.stop,
                                    canonical_bin.left_label, 
                                    canonical_bin.right_label, 
                                    "EXON"  ) )
        
        pseudo_exons.extend( find_left_exon_extensions(
                ce_i, canonical_bin, gene_bins, rnaseq_cov))

        
        pseudo_exons.extend( find_right_exon_extensions(
                ce_i, canonical_bin, gene_bins, rnaseq_cov))
    
    cage_peak_bin_indices = []
    tss_exons = []
    for peak in cage_peaks:
        # find the first donor junction right of the cage peaks
        for bin_i, bin in enumerate(gene_bins):
            if bin.stop < peak.stop: continue
            if bin.right_label in ('D_JN', 'POLYA', 'ESTART') : break
            cage_peak_bin_indices.append( bin_i )
            tss_exons.append( Bin( peak.start, bin.stop,
                                   "CAGE_PEAK", 
                                   bin.right_label, 
                                   "TSS_EXON"  ) )
    
    for cage_peak_i in cage_peak_bin_indices:
        bin = gene_bins[cage_peak_i]
        pseudo_exons.extend( find_right_exon_extensions(
                cage_peak_i, bin, gene_bins, rnaseq_cov))
    
    # find tes exons
    tes_exon_indices = [ i for i, bin in enumerate( gene_bins )
                         if bin.right_label in ('POLYA', ) ]

    tes_exons = []
    for tes_exon_i in tes_exon_indices:
        bin = gene_bins[tes_exon_i]
        
        tes_exons.append( Bin( bin.start, bin.stop,
                               bin.left_label, 
                               bin.right_label, 
                               "TES_EXON"  ) )
        
        pseudo_exons.append( tes_exons[-1] )
        pseudo_exons.extend( find_left_exon_extensions(
                tes_exon_i, bin, gene_bins, rnaseq_cov ))

    # build exons from the pseudo exon set
    pseudo_exons = list( set( pseudo_exons ) )
    pseudo_exons.sort( key=lambda x: (x.start, x.stop ) )
    
    return cage_peaks, tss_exons, pseudo_exons, tes_exons

def find_exons_in_gene( ( chrm, strand, contig_len ), gene, 
                        rnaseq_cov, cage_cov, 
                        polya_sites, jns ):
    ###########################################################
    # Shift all of the input data to eb in the gene region, and 
    # reverse it when necessary
    polya_sites = [x - gene.start for x in polya_sites
                   if x > gene.start and x <= gene.stop]
        
    jns = [ (x1 - gene.start, x2 - gene.start, cnt)  
            for x1, x2, cnt in jns ]
    
    if strand == '-':
        gene_len = gene.stop - gene.start + 1
        polya_sites = [ gene_len - x for x in polya_sites ]
        jns = [ (gene_len-x2, gene_len-x1, cnt) for x1, x2, cnt in jns ]
        rnaseq_cov = rnaseq_cov[::-1]
        cage_cov = cage_cov[::-1]
    
    ### END Prepare input data #########################################
    
    cage_peaks = find_cage_peaks_in_gene( 
        ( chrm, strand ), gene, cage_cov, rnaseq_cov )
    
    cage_peaks, tss_pseudo_exons, pseudo_exons, tes_pseudo_exons = \
        find_pseudo_exons_in_gene(
            ( chrm, strand ), gene, 
            cage_peaks, polya_sites, jns,
            rnaseq_cov, cage_cov  )
    
    # find the contiguous set of adjoining exons
    slices = [[0,],]
    for i, ( start_ps_exon, stop_ps_exon ) in enumerate( 
            zip(pseudo_exons[:-1], pseudo_exons[1:]) ):
        # if there is a gap...
        if start_ps_exon.stop < stop_ps_exon.start:
            slices[-1].append(i+1)
            slices.append([i+1,])
    
    # add in the last segment
    if slices[-1][0] < len( pseudo_exons ):
        slices[-1].append( len(pseudo_exons) )
    else:
        slices.pop()
    
    gpd_pseudo_exons = []    
    for start, stop in slices:
        gpd_pseudo_exons.append( pseudo_exons[start:stop] )
        
    if len( gpd_pseudo_exons ) == 0:
        print "============================================= NO EXONS IN REGION"
        print gene
        return Bins(chrm, strand, [] )
            
    
    internal_exons = Bins(chrm, strand, [] )
    for pseudo_exons_grp in gpd_pseudo_exons:
        internal_exons.extend( 
            find_internal_exons_from_pseudo_exons( pseudo_exons_grp ) )

    tss_exons = Bins(chrm, strand, [] )
    for tss_exon in tss_pseudo_exons:
        if tss_exon.right_label in 'D_JN':
            tss_exons.append( tss_exon )
        elif tss_exon.right_label == 'POLYA':
            tss_exons.append( tss_exon )
        else:
            pass
            #try:
            #    assert tss_exon.stop >= gpd_pseudo_exons[0][0].start
            #except:
            #    print tss_exon
            #    print gpd_pseudo_exons[0][0]
            #    raise
        
        for pe_grp in gpd_pseudo_exons:
            if pe_grp[0].start > tss_exon.stop:  
                continue
            if pe_grp[-1].stop <= tss_exon.start:
                break
            
            tss_exons.extend( find_tss_exons_from_pseudo_exons( 
                    tss_exon, pe_grp ) )

    single_exon_genes = set()
    for exon in tss_exons:
        if exon.right_label == 'POLYA':
            se_gene = copy(exon)
            se_gene.type = 'SE_GENE'
            single_exon_genes.add( se_gene )
    
    tes_exons = Bins(chrm, strand, [] )
    for tes_exon in tes_pseudo_exons:
        if tes_exon.left_label == 'R_JN':
            tes_exons.append( tes_exon )
        elif tes_exon.left_label == 'POLYA':
            pass
        else:
            assert tes_exon.start <= gpd_pseudo_exons[-1][-1].stop
        
        # find the overlapping intervals
        for pe_grp in gpd_pseudo_exons:
            if pe_grp[0].start >= tes_exon.stop:  
                break
            if pe_grp[-1].stop <= tes_exon.start:
                continue
            if tes_exon == pe_grp[0]:
                continue
            
            # now we know that the tes exon overlaps the pseudo exon grp
            #print tes_exon, [ (pe.left_label, pe.right_label) for pe in pe_grp ]
            #print pe_grp
            
            # find the psuedo exon that shares a boundary
            last_pe_i = max( i for i, pe in enumerate( pe_grp ) 
                             if pe.stop == tes_exon.start )

            for pe in pe_grp[:last_pe_i+1]:
                # we had this earlier, buyt I think we want tes exons to
                # extend into real exons ( although this gives us fragments )
                #if pe.left_label != 'R_JN':
                #    continue

                tes_exons.append( 
                    Bin( pe.start, tes_exon.stop, pe.left_label, 
                         tes_exon.right_label, "TES_EXON" )
                    )

    jn_bins = Bins(chrm, strand, [])
    for start, stop, cnt in jns:
        bin = Bin(start, stop, 'R_JN', 'D_JN', 'INTRON', cnt)
        jn_bins.append( bin )

    tss_exons = filter_exons( tss_exons, rnaseq_cov )
    tes_exons = filter_exons( tes_exons, rnaseq_cov )
    internal_exons = filter_exons( internal_exons, rnaseq_cov )

    elements = Bins(chrm, strand, chain(
            jn_bins, cage_peaks, tss_exons, internal_exons, tes_exons) )
    if strand == '-':
        elements = elements.reverse_strand( gene.stop - gene.start + 1 )
    elements = elements.shift( gene.start )
    elements.append( gene )
    
    return elements
        
        

def find_exons_in_contig( (chrm, strand, contig_len), ofp,
                          rnaseq_reads, cage_reads, polya_sites):
    junctions = load_junctions( rnaseq_reads, chrm, strand )
    polya_sites = polya_sites[(chrm, strand)]
    polya_bins = Bins( chrm, strand, [] )
    for x in polya_sites:
        polya_bins.append( Bin( x-10, x, "POLYA", "POLYA", "POLYA" ) )
    write_unified_bed( polya_bins, ofp)
    
    gene_bndry_bins = find_gene_boundaries( 
       (chrm, strand, contig_len), rnaseq_reads, 
       polya_sites, junctions)
        
    sorted_jns = sorted( junctions )
    jn_starts = [ i[0][0] for i in sorted_jns ]
    jn_stops = [ i[0][1] for i in sorted_jns ]
    jn_values = [ i[1] for i in sorted_jns ]
    
    for gene in gene_bndry_bins:
        # find the junctions associated with this gene
        gj_sa = bisect( jn_stops, gene.start )
        gj_so = bisect( jn_starts, gene.stop )
        gene_jns = zip( jn_starts[gj_sa:gj_so], 
                        jn_stops[gj_sa:gj_so], 
                        jn_values[gj_sa:gj_so] )
        
        # find the polyas associated with this gene
        gene_polya_sa_i = polya_sites.searchsorted( gene.start )
        gene_polya_so_i = polya_sites.searchsorted( gene.stop, side='right' )
        gene_polyas = polya_sites[gene_polya_sa_i:gene_polya_so_i]
    
        gene_rnaseq_cov = rnaseq_reads[0].build_read_coverage_array( 
            chrm, strand, gene.start, gene.stop )

        gene_cage_cov = cage_reads[0].build_read_coverage_array( 
            'chr' + chrm, strand, gene.start, gene.stop+1 )
        
        elements = \
            find_exons_in_gene( ( chrm, strand, contig_len ), gene, 
                                gene_rnaseq_cov, gene_cage_cov, 
                                gene_polyas, gene_jns )
        
        write_unified_bed( elements, ofp)
        if VERBOSE: print "FINISHED ", gene
    
    return

def main():
    rnaseq_bams, cage_bams, polya_candidate_sites_fps, ofp \
        = parse_arguments()

    ofp.write('track name="discovered_elements" visibility=2 itemRgb="On"\n')
    
    rnaseq_reads = [ RNAseqReads(fp.name).init(reverse_read_strand=True) 
                     for fp in rnaseq_bams ]
    
    cage_reads = [ RAMPAGEReads(fp.name).init(reverse_read_strand=False) 
                   for fp in cage_bams ]

    if VERBOSE: print >> sys.stderr,  'Loading candidate polyA sites'
    polya_sites = find_polya_sites([x.name for x in polya_candidate_sites_fps])
    for fp in polya_candidate_sites_fps: fp.close()
    if VERBOSE: print >> sys.stderr, 'Finished loading candidate polyA sites'
    
    contig_lens = get_contigs_and_lens( rnaseq_reads, cage_reads )
    for contig, contig_len in contig_lens.iteritems():
        if contig != "4": continue
        for strand in '+-':
            find_exons_in_contig( (contig, strand, contig_len), ofp,
                                  rnaseq_reads, cage_reads, polya_sites)
    
if __name__ == '__main__':
    main()



















#
