import sys

import numpy

from collections import defaultdict
from itertools import chain, izip
from bisect import bisect
from copy import copy

import multiprocessing
import Queue

from igraph import Graph

from files.reads import RNAseqReads, CAGEReads, RAMPAGEReads, clean_chr_name, \
    guess_strand_from_fname, iter_coverage_intervals_for_read
from files.junctions import extract_junctions_in_region, extract_junctions_in_contig
from files.bed import create_bed_line
from files.gtf import parse_gtf_line

USE_CACHE = False
NTHREADS = 1

class ThreadSafeFile( file ):
    def __init__( self, *args ):
        args = list( args )
        args.insert( 0, self )
        file.__init__( *args )
        self.lock = multiprocessing.Lock()

    def write( self, line ):
        self.lock.acquire()
        file.write( self, line )
        self.flush()
        self.lock.release()

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
MIN_EMPTY_REGION_LEN = 100
MIN_EXON_BPKM = 1
EXON_EXT_CVG_RATIO_THRESH = 5
POLYA_MERGE_SIZE = 100

CAGE_PEAK_WIN_SIZE = 20
MIN_NUM_CAGE_TAGS = 10
MAX_CAGE_FRAC = 0.01

def get_contigs_and_lens( rnaseq_reads, cage_reads ):
    """Get contigs and their lengths from a set of bam files.
    
    We make sure that the contig lengths are consistent in all of the bam files, and
    we remove contigs that dont have at least 1 read in at least one rnaseq file 
    and one promoter reads file.
    """
    chrm_lengths = {}
    contigs = None
    for bam in chain( rnaseq_reads, cage_reads ):
        bam_contigs = set()
        for ref_name, ref_len in zip(bam.references, bam.lengths):
            # add the contig to the chrm lengths file, checking to
            # make sure that lengths match if it has already been added
            if clean_chr_name( ref_name ) not in chrm_lengths:
                chrm_lengths[clean_chr_name( ref_name )] = ref_len
            else:
                assert chrm_lengths[clean_chr_name(ref_name)] == ref_len, \
                    "Chromosome lengths do not match between bam files"
            bam_contigs.add( clean_chr_name(ref_name) )
        
        if contigs == None:
            contigs = bam_contigs
        else:
            contigs = contigs.intersection( bam_contigs )
    
    # remove contigs that dont have reads in at least one file
    def at_least_one_bam_has_reads( chrm, bams ):
        for bam in  chain( rnaseq_reads, cage_reads ):
            try:
                next( bam.fetch( chrm ) )
            except StopIteration:
                continue
            except KeyError:
                continue
            else:
                return True
    
    # produce the final list of contigs
    rv =  {}
    for key, val in chrm_lengths.iteritems():
        if at_least_one_bam_has_reads(key, rnaseq_reads) \
                and at_least_one_bam_has_reads(key, cage_reads) \
                and key in contigs:
            rv[key] = val
    
    return rv

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
    Q.sort()
    return Q[int(len(Q)*0.1)], Q[int(len(Q)*0.9)]

def get_qrange_short( np ):
    L = len(np)
    return np.min(), np.max()

def filter_exon( exon, wig, min_avg_cvg=0.01, 
                 min_short_cvg=1, short_exon_length=400,
                 min_long_cvg=10 ):
    '''Find all the exons that are sufficiently homogenous and expressed.
    
    '''
    start = exon.start
    end = exon.stop
    vals = wig[start:end+1]
    n_div = max( 1, int(len(vals)/100) )
    div_len = len(vals)/n_div
    for i in xrange(n_div):
        seg = vals[i*div_len:(i+1)*div_len]
        if seg.mean() < MIN_EXON_BPKM:
            return True

    return False
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
    if IQ_ratio < 20:
        return False
    
    # if the IQ range is "boarder line", but the low is high, keep it, what the hell
    if IQ_ratio >= 20 and IQ_ratio < 100 and low > 50:
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
    
    def length( self ):
        return self.stop - self.start + 1
    
    def mean_cov( self, cov_array ):
        return numpy.median(cov_array[self.start:self.stop])
        return cov_array[self.start:self.stop].mean()
    
    def reverse_strand(self, contig_len):
        return Bin(contig_len-self.stop, contig_len-self.start, 
                   self.right_label, self.left_label, self.type)

    def shift(self, shift_amnt):
        return Bin(self.start+shift_amnt, self.stop+shift_amnt, 
                   self.left_label, self.right_label, self.type)
    
    def __repr__( self ):
        if self.type == None:
            return "%i-%i:%s:%s" % ( self.start, self.stop, self.left_label,
                                       self.right_label )

        return "%s:%i-%i" % ( self.type, self.start, self.stop )
    
    def __hash__( self ):
        return hash( (self.start, self.stop, self.type, 
                      self.left_label, self.right_label) )
    
    def __eq__( self, other ):
        return ( self.start == other.start and self.stop == other.stop )
    
    _bndry_color_mapping = {
        'CONTIG_BNDRY': '0,0,0',
        'GENE_BNDRY': '0,0,0',
        
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
        'SE_GENE': '000,000,200',
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
                                        bin.start-1, bin.stop, 
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

    def shift(self, shift_amnt ):
        shifted_bins = Bins( self.chrm, self.strand )
        for bin in self:
            shifted_bins.append( bin.shift( shift_amnt ) )
        return shifted_bins
    
    def writeBed( self, ofp ):
        """
            chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
        """
        for bin in self:
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

    def writeGff( self, ofp ):
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

def load_junctions_worker(all_jns, all_jns_lock, args):
    jns = extract_junctions_in_region( *args )
    all_jns_lock.acquire()
    all_jns.extend( jns )
    all_jns_lock.release()
    del jns
    return

def load_junctions( rnaseq_reads, (chrm, strand, contig_len) ):
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
        if NTHREADS == 1:
            junctions = extract_junctions_in_contig( 
                rnaseq_reads[0], chrm, strand )
        else:
            seg_len = int(contig_len/NTHREADS)
            segments =  [ [i*seg_len, (i+1)*seg_len] for i in xrange(NTHREADS) ]
            segments[0][0] = 0
            segments[-1][1] = contig_len
            
            from multiprocessing import Process, Manager
            manager = Manager()
            all_jns = manager.list()
            all_jns_lock = manager.Lock()
            
            ps = []
            for start, stop in segments:
                p = Process(target=load_junctions_worker,
                            args=( all_jns, all_jns_lock, 
                                   (rnaseq_reads[0], chrm, strand, start, stop, True) 
                                   ) )
                
                p.start()
                ps.append( p )
            
            import time
            while True:
                if all( not p.is_alive() for p in ps ):
                    break
                time.sleep( 0.1 )
            
            junctions = defaultdict( int )
            for jn, cnt in all_jns:
                junctions[jn] += cnt
            
            junctions = sorted(junctions.iteritems())
    
    if VERBOSE: print "Finished extracting junctions for %s %s" % (chrm, strand)
    
    # filter junctions
    jn_starts = defaultdict( int )
    jn_stops = defaultdict( int )
    for (start, stop), cnt in junctions:
        jn_starts[start] = max( jn_starts[start], cnt )
        jn_stops[stop] = max( jn_stops[stop], cnt )
    
    filtered_junctions = []
    for (start, stop), cnt in junctions:
        if float(cnt)/jn_starts[start] < 0.01: continue
        if float(cnt)/jn_stops[stop] < 0.01: continue
        if stop - start > 10000000: continue
        filtered_junctions.append( ((start, stop), cnt) )
    
    return filtered_junctions

def find_gene_boundaries((chrm, strand, contig_len), rnaseq_reads, 
                         polya_sites, junctions=None):
    def find_splits_in_region(rnaseq_reads, chrm, strand, start, stop):
        n_splits = int( (stop-start)/MIN_REGION_LEN )
        segment_len = int((stop-start)/(n_splits+1))
        return [ start + i*segment_len for i in xrange(1,n_splits+1) ]
    
    def find_initial_segmentation( chrm, strand, rnaseq_reads, polya_sites ):
        locs = {0: 'SEGMENT', contig_len:'SEGMENT'}
        for polya in polya_sites:
            locs[ polya ] = "POLYA"
        return locs
    
    def merge_segments( segments, strand, window_len = 10 ):
        bndries = numpy.array( sorted(segments, reverse=(strand!='-')) )
        bndries_to_delete = set()
        for i, bndry in enumerate(bndries[1:-1]):
            if bndry - bndries[i+1-1] > 1000000:
                continue
            pre_bndry_cnt = 0
            post_bndry_cnt = 0
            
            for reads in rnaseq_reads:
                cvg = reads.build_read_coverage_array( 
                    chrm, strand, max(0,bndry-window_len), bndry+window_len )
                pre_bndry_cnt += cvg[:window_len].sum()
                post_bndry_cnt += cvg[window_len:].sum()
            
            if strand == '-':
                pre_bndry_cnt, post_bndry_cnt = post_bndry_cnt, pre_bndry_cnt
            if (post_bndry_cnt)/(post_bndry_cnt+pre_bndry_cnt+1e-6) > 0.20:
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
                edges.add((int(min(start_bin, stop_bin)), 
                           int(max(start_bin, stop_bin))))
        
        genes_graph = Graph( len( boundaries )-1 )
        genes_graph.add_edges( list( edges ) )
        
        segments = []
        for g in genes_graph.clusters():
            segments.append( (boundaries[min(g)]+1, boundaries[max(g)+1]-1) )
        
        return flatten( segments )
    
    # find all of the junctions
    if None == junctions:
        junctions = load_junctions( rnaseq_reads, (chrm, strand, contig_len) )
    
    # find segment boundaries
    initial_segmentation = find_initial_segmentation( 
        chrm, strand, rnaseq_reads, polya_sites )
    if VERBOSE: print "Finished initial segmentation for %s %s" % (chrm, strand)    
    merged_segments = merge_segments( initial_segmentation, strand )
    if VERBOSE: print "Finished merging segments for %s %s" % (chrm, strand)
    clustered_segments = cluster_segments( merged_segments, junctions )
    if VERBOSE: print "Finished clustering segments for %s %s" % (chrm, strand)    
    
    # build the gene bins, and write them out to the elements file
    genes = Bins( chrm, strand, [] )
    if len( clustered_segments ) == 2:
        return genes
    
    for start, stop in clustered_segments:
        if stop - start < 100: continue
        if stop - start > 1000000: continue
        genes.append( Bin(max(1,start-10), min(stop+10,contig_len), 
                          "ESTART", "POLYA", "GENE" ) )
    
    return genes

def filter_polya_sites( (chrm, strand), gene, polya_sites, rnaseq_cov ):
    if len(polya_sites) == 0:
        return polya_sites
    
    polya_sites.sort()

    new_polya_sites = []
    for site in polya_sites[:-1]:
        pre_cvg = rnaseq_cov[max(0,site-10):site].sum()       
        post_cvg = rnaseq_cov[site+10:site+20].sum()
        if pre_cvg > 10 and pre_cvg/(post_cvg+1.0) < 5:
            continue
        else:
            new_polya_sites.append( site )
    
    new_polya_sites.append( polya_sites[-1] )
    polya_sites = new_polya_sites

    # merge sites that are close
    new_polya_sites = [polya_sites[0],]
    for site in polya_sites[1:]:
        if site - new_polya_sites[-1] < 20:
            new_polya_sites[-1] = site
        else:
            new_polya_sites.append( site )
    
    return new_polya_sites


def find_cage_peaks_in_gene( ( chrm, strand ), gene, cage_cov, rnaseq_cov ):
     raw_peaks = find_peaks( cage_cov, window_len=CAGE_PEAK_WIN_SIZE, 
                             min_score=MIN_NUM_CAGE_TAGS,
                             max_score_frac=MAX_CAGE_FRAC, 
                             max_num_peaks=20 )
     if len( raw_peaks ) == 0:
         return []
     
     cage_peaks = Bins( chrm, strand )
     for peak_st, peak_sp in raw_peaks:
         # make sure there is *some* rnaseq coverage post peak
         if rnaseq_cov[peak_st:peak_sp+10].sum() < MIN_NUM_CAGE_TAGS: continue
         # make sure that there is an increase in coverage from pre to post peak
         pre_peak_cov = rnaseq_cov[peak_st-10:peak_st].sum()
         post_peak_cov = rnaseq_cov[peak_st:peak_sp+10].sum()
         if post_peak_cov/(pre_peak_cov+1e-6) < 5: continue
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
    def grow_peak( start, stop, grow_size=max(1, window_len/4), min_grow_ratio=0.5 ):
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
        
        # break at canonical introns
        if bin.left_label == 'D_JN' and bin.right_label == 'R_JN':
            break
        
        # make sure the average coverage is high enough
        bin_cvg = bin.mean_cov(rnaseq_cov)   
        if bin_cvg < MIN_EXON_BPKM:
            break

        if bin.stop - bin.start > 20 and \
                (start_bin_cvg+1e-6)/(bin_cvg+1e-6) > EXON_EXT_CVG_RATIO_THRESH:
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

def find_right_exon_extensions( start_index, start_bin, gene_bins, rnaseq_cov,
                                min_ext_ratio=EXON_EXT_CVG_RATIO_THRESH, 
                                min_bpkm=MIN_EXON_BPKM):
    exons = []
    ee_indices = []
    start_bin_cvg = start_bin.mean_cov( rnaseq_cov )
    for i in xrange( start_index+1, len(gene_bins) ):
        bin = gene_bins[i]

        # if we've reached a canonical intron, break
        #if bin.left_label == 'D_JN' and bin.right_label == 'R_JN':
        #    break
        if bin.left_label == 'D_JN' and bin.right_label == 'R_JN':
            break
        
        # make sure the average coverage is high enough
        bin_cvg = bin.mean_cov(rnaseq_cov)
        if bin_cvg < min_bpkm:
            break
        
        if bin.stop - bin.start > 20 and \
                (start_bin_cvg+1e-6)/(bin_cvg+1e-6) > min_ext_ratio:
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


def build_labeled_segments( rnaseq_cov, jns, cage_peaks=[], polya_sites=[] ):
    locs = defaultdict(set)    
    for start, stop in find_empty_regions( rnaseq_cov ):
        if stop - start < MIN_EMPTY_REGION_LEN: continue
        locs[start].add( "ESTART" )
        locs[stop].add( "ESTOP" )
    
    for start, stop, cnt in jns:
        if start < 1 or stop > len(rnaseq_cov): continue
        #assert start-1 not in locs, "%i in locs" % (start-1)
        locs[start-1].add( "D_JN" )
        locs[stop+1].add( "R_JN" )

    for polya in sorted(polya_sites, reverse=True):
        locs[ polya ].add( "POLYA" )
    
    # build all of the bins
    poss = sorted( locs.iteritems() )
    poss = merge_empty_labels( poss )
    if len( poss ) == 0: 
        return []

    if poss[0][0] > 1:
        poss.insert( 0, (1, set(["GENE_BNDRY",])) )
    if poss[-1][0] < len(rnaseq_cov)-1:
        poss.append( (len(rnaseq_cov)-1, set(["GENE_BNDRY",])) )
    
    bins = []
    for index, ((start, left_labels), (stop, right_labels)) in \
            enumerate(izip(poss[:-1], poss[1:])):
        for left_label in left_labels:
            for right_label in right_labels:
                bins.append( Bin(start, stop, left_label, right_label) )
    
    return bins

def find_canonical_and_internal_exons( (chrm, strand), rnaseq_cov, jns, cage_peaks, polya_sites  ):
    bins = build_labeled_segments( rnaseq_cov, jns )    
    
    def iter_canonical_exons_and_indices():
        for i, bin in enumerate( bins ):
            if bin.left_label == 'R_JN' and bin.right_label == 'D_JN':
                yield i, bin
    
    canonical_exons = Bins( chrm, strand )
    internal_exons = Bins( chrm, strand )
    for ce_i, ce_bin in iter_canonical_exons_and_indices():
        ce_bin.type = 'EXON'
        canonical_exons.append( ce_bin )
        internal_exons.append( ce_bin )
        
        for r_ext in find_right_exon_extensions(
                ce_i, ce_bin, bins, rnaseq_cov):
            exon = copy(ce_bin)
            exon.right_label = r_ext.right_label
            exon.stop = r_ext.stop
            internal_exons.append( exon )

        for l_ext in find_left_exon_extensions(
                ce_i, ce_bin, bins, rnaseq_cov):
            exon = copy(ce_bin)
            exon.left_label = l_ext.left_label
            exon.start = l_ext.start
            internal_exons.append( exon )
    
    return canonical_exons, internal_exons

def find_tss_exons_and_se_genes( 
        (chrm, strand), rnaseq_cov, jns, polya_sites, cage_peaks ):
    bins = build_labeled_segments( rnaseq_cov, jns, polya_sites=polya_sites )
    tss_exons = Bins( chrm, strand )
    se_genes = Bins( chrm, strand )
    if len(bins) == 0: 
        return tss_exons, se_genes

    for peak in cage_peaks:
        # should this be peak start
        try:
            tss_index, tss_bin = next( (i, bin) for i, bin in enumerate(bins)
                                       if bin.stop > peak.stop
                                       and bin.right_label in ( 'POLYA', 'D_JN') )
        except StopIteration:
            continue
        
        tss_exon = Bin( 
            peak.start, tss_bin.stop, 
            "CAGE_PEAK", tss_bin.right_label, "TSS_EXON" )
        if tss_exon.right_label == 'POLYA':
            tss_exon.type = 'SE_GENE'
            se_genes.append( tss_exon )
        elif tss_exon.right_label == 'D_JN':
            tss_exons.append( tss_exon )
        
        for r_ext in find_right_exon_extensions(
                tss_index, tss_bin, bins, rnaseq_cov):
            tss_exon = Bin( 
                peak.start, r_ext.stop, 
                "CAGE_PEAK", r_ext.right_label, "TSS_EXON" )
            if tss_exon.right_label == 'POLYA':
                tss_exon.type = 'SE_GENE'
                se_genes.append( tss_exon )
            elif tss_exon.right_label == 'D_JN':
                tss_exons.append( tss_exon )
    
    return tss_exons, se_genes

def find_tes_exons( (chrm, strand), rnaseq_cov, jns, polya_sites ):
    bins = build_labeled_segments( rnaseq_cov, jns, polya_sites=polya_sites )
    
    # find the segments that contain a 
    canonical_tes_exon_indices = [ i for i, bin in enumerate( bins )
                                   if bin.right_label == 'POLYA'
                                   and bin.left_label in ('R_JN', 'D_JN')]

    tes_exons = Bins( chrm, strand )
    tes_exon_indices = []
    for tes_exon_i in canonical_tes_exon_indices:
        tes_bin = bins[tes_exon_i]
        tes_bin.type = 'TES_EXON'
        if tes_bin.left_label == 'R_JN':
            tes_exons.append( tes_bin )
            tes_exon_indices.append( (tes_exon_i, tes_exon_i) )
        
        for l_ext_i, l_ext in enumerate(find_left_exon_extensions(
                tes_exon_i, tes_bin, bins, rnaseq_cov)):
            if l_ext.left_label == 'R_JN':
                tes_exon = Bin( 
                    l_ext.start, tes_bin.stop, 
                    l_ext.left_label, tes_bin.right_label, "TES_EXON" )
                tes_exons.append( tes_exon )
                tes_exon_indices.append( (tes_exon_i, tes_exon_i+(l_ext_i+1)) )
    
    for tes_exon, (start_i, stop_i) in izip( tes_exons, tes_exon_indices ):
        for r_ext in bins[stop_i+1:]:
            if r_ext.right_label != 'POLYA': break
            if r_ext.mean_cov(rnaseq_cov) < MIN_EXON_BPKM: break
            # filter extensions with a sufficiently low coverage ratio
            if (tes_exon.mean_cov(rnaseq_cov)+1e-6)/\
                    (r_ext.mean_cov(rnaseq_cov)+1e-6) > 10:
                continue
            tes_exon.stop = r_ext.stop
    
    return Bins( chrm, strand, sorted(set(tes_exons)))

def find_exons_in_gene( ( chrm, strand, contig_len ), gene, 
                        rnaseq_reads, cage_reads, 
                        polya_sites, jns ):
    if DEBUG_VERBOSE: 
        print "STARTING Chrm %s Strand %s Pos %i-%i" % (
            chrm, strand, gene.start, gene.stop)
    ###########################################################
    # Shift all of the input data to eb in the gene region, and 
    # reverse it when necessary
    polya_sites = [x - gene.start for x in polya_sites
                   if x > gene.start and x <= gene.stop]
        
    jns = [ (x1 - gene.start, x2 - gene.start, cnt)  
            for x1, x2, cnt in jns ]

    rnaseq_cov = rnaseq_reads[0].build_read_coverage_array( 
        chrm, strand, gene.start, gene.stop )
    
    cage_cov = cage_reads[0].build_read_coverage_array( 
        chrm, strand, gene.start, gene.stop+1 )
    
    if strand == '-':
        gene_len = gene.stop - gene.start + 1
        polya_sites = [ gene_len - x for x in polya_sites ]
        jns = [ (gene_len-x2, gene_len-x1, cnt) for x1, x2, cnt in jns ]
        rnaseq_cov = rnaseq_cov[::-1]
        cage_cov = cage_cov[::-1]
    
    filtered_junctions = []
    for (start, stop, cnt) in jns:
        if start < 0 or stop > contig_len: continue
        left_intron_cvg = rnaseq_cov[start+10:start+30].sum()/20
        right_intron_cvg = rnaseq_cov[stop-30:stop-10].sum()/20        
        if cnt*10 < left_intron_cvg or cnt*10 < right_intron_cvg:
            continue
        filtered_junctions.append( (start, stop, cnt) )

    jns = filtered_junctions
    
    ### END Prepare input data #########################################
    
    cage_peaks = find_cage_peaks_in_gene( 
        (chrm, strand), gene, cage_cov, rnaseq_cov )
    #polya_sites = filter_polya_sites( (chrm, strand), gene, polya_sites, rnaseq_cov )
    
    canonical_exons, internal_exons = find_canonical_and_internal_exons(
        (chrm, strand), rnaseq_cov, jns, cage_peaks, polya_sites)
    tss_exons, se_genes = find_tss_exons_and_se_genes( 
        (chrm, strand), rnaseq_cov, jns, polya_sites, cage_peaks)
    tes_exons = find_tes_exons( 
        (chrm, strand), rnaseq_cov, jns, polya_sites)

    gene_bins = Bins(chrm, strand, build_labeled_segments( 
            rnaseq_cov, jns, polya_sites=polya_sites ) )
    
    jn_bins = Bins(chrm, strand, [])
    for start, stop, cnt in jns:
        bin = Bin(start, stop, 'R_JN', 'D_JN', 'INTRON', cnt)
        jn_bins.append( bin )
    
    tss_exons = filter_exons( tss_exons, rnaseq_cov )
    tes_exons = filter_exons( tes_exons, rnaseq_cov )
    internal_exons = filter_exons( internal_exons, rnaseq_cov )
    se_genes = filter_exons( se_genes, rnaseq_cov )

    #print len(tss_exons), len(internal_exons), len(tes_exons)
    
    elements = Bins(chrm, strand, chain(
            jn_bins, cage_peaks, tss_exons, internal_exons, tes_exons, se_genes) )
    if strand == '-':
        elements = elements.reverse_strand( gene.stop - gene.start + 1 )
    elements = elements.shift( gene.start )
    elements.append( gene )

    if strand == '-':
        gene_bins = gene_bins.reverse_strand( gene.stop - gene.start + 1 )
    gene_bins = gene_bins.shift( gene.start )
    
    return elements, gene_bins
        

def find_exons_worker( (genes_queue, genes_queue_lock), ofp, (chrm, strand, contig_len),
                        jns, rnaseq_reads, cage_reads, polya_sites ):
    jn_starts = [ i[0][0] for i in jns ]
    jn_stops = [ i[0][1] for i in jns ]
    jn_values = [ i[1] for i in jns ]
    
    rnaseq_reads = [ x.reload() for x in rnaseq_reads ]
    cage_reads = [ x.reload() for x in cage_reads ]
    
    while True:
        if genes_queue_lock != None:
            genes_queue_lock.acquire()
            if len(genes_queue) == 0:
                genes_queue_lock.release()
                break
            gene = genes_queue.pop()
            genes_queue_lock.release()
        else:
            if len( genes_queue ) == 0: 
                break
            assert NTHREADS == 1
            gene = genes_queue.pop()
        
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
                
        elements, pseudo_exons = \
            find_exons_in_gene( ( chrm, strand, contig_len ), gene, 
                                rnaseq_reads, cage_reads,
                                gene_polyas, gene_jns )
    
        elements.extend( Bin( x-10, x, "POLYA", "POLYA", "POLYA" ) 
                         for x in sorted(gene_polyas) )
        write_unified_bed( elements, ofp)
        
        if WRITE_DEBUG_DATA:
            pseudo_exons.writeBed( ofp )
        
        if DEBUG_VERBOSE: print "FINISHED Chrm %s Strand %s Pos %i-%i" % (chrm, strand, gene.start, gene.stop)
    
    return

def find_exons_in_contig( (chrm, strand, contig_len), ofp,
                          rnaseq_reads, cage_reads, polya_sites,
                          gene_bndry_bins=None):
    junctions = load_junctions( rnaseq_reads, (chrm, strand, contig_len) )
    polya_sites = polya_sites[(chrm, strand)]
    
    if gene_bndry_bins == None:
        gene_bndry_bins = find_gene_boundaries( 
            (chrm, strand, contig_len), rnaseq_reads, 
            polya_sites, junctions)
    
    if NTHREADS > 1:
        manager = multiprocessing.Manager()
        genes_queue = manager.list()
        genes_queue_lock = manager.Lock()
    else:
        genes_queue, genes_queue_lock = [], None
        
    genes_queue.extend( gene_bndry_bins )
    sorted_jns = sorted( junctions )
    args = [ (genes_queue, genes_queue_lock), ofp, (chrm, strand, contig_len),
             sorted_jns, rnaseq_reads, cage_reads, polya_sites ]

    if NTHREADS == 1:
        find_exons_worker(*args)
    else:
        ps = []
        for i in xrange( NTHREADS ):
            p = multiprocessing.Process(target=find_exons_worker, args=args)
            p.start()
            ps.append( p )

        for p in ps:
            p.join()
    
    return
def load_gene_bndry_bins( gtf_fname, contig_lens ):
    gene_bndry_bins = {}
    with open(gtf_fname) as fp:
        for line in fp:
            data = parse_gtf_line( line )
            if data.region.chr not in contig_lens:
                continue
            if data == None or data.feature != 'gene': 
                continue
            key = (data.region.chr, data.region.strand) 
            if key not in gene_bndry_bins:
                gene_bndry_bins[key] = Bins( data.region.chr, data.region.strand )
            gene_bin = Bin(max(1,data.region.start-2500), 
                           min(data.region.stop+2500, 
                               contig_lens[data.region.chr]), 
                           'GENE', 'GENE', 'GENE' )
            gene_bndry_bins[key].append( gene_bin )
    
    merged_gene_bndry_bins = gene_bndry_bins
    """
    for key, bins in gene_bndry_bins.iteritems():
        merged_gene_bndry_bins[key] = Bins( *key )
        sorted_bins = sorted(bins)
        merged_gene_bndry_bins[key].append( sorted_bins[0] )
        for bin in sorted_bins[1:]:
            if bin.start > merged_gene_bndry_bins[key][-1].stop:
                merged_gene_bndry_bins[key].append( bin )
            else:
                merged_gene_bndry_bins[key][-1].start = min( 
                    merged_gene_bndry_bins[key][-1].start, bin.start )
                merged_gene_bndry_bins[key][-1].stop = max( 
                    merged_gene_bndry_bins[key][-1].stop, bin.stop )
    """
    return merged_gene_bndry_bins

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(\
        description='Find exons from RNAseq, CAGE, and poly(A) assays.')

    parser.add_argument( 'rnaseq_reads',type=argparse.FileType('rb'),nargs='+',\
        help='BAM files containing mapped RNAseq reads ( must be indexed ).')

    parser.add_argument( '--reverse-rnaseq-strand', default=False,
                         action='store_true',
        help='Whether to reverse the RNAseq read strand (default False).')
    
    parser.add_argument( '--cage-reads', type=file, default=[], nargs='*', \
        help='BAM files containing mapped cage reads.')

    parser.add_argument( '--rampage-reads', type=file, default=[], nargs='*', \
        help='BAM files containing mapped rampage reads.')
    
    parser.add_argument( '--polya-candidate-sites', type=file, nargs='*', \
        help='files with allowed polya sites.')

    parser.add_argument( '--reference',
        help='Reference GTF ( for gene boundry information )')
    
    parser.add_argument( '--ofname', '-o', 
                         default="discovered_elements.bed",\
        help='Output file name. (default: discovered_elements.bed)')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
        help='Whether or not to print status information.')
    parser.add_argument( '--debug-verbose', default=False, action='store_true',\
        help='Whether or not to print debugging information.')
    parser.add_argument('--write-debug-data',default=False,action='store_true',\
        help='Whether or not to print out gff files containing intermediate exon assembly data.')
    parser.add_argument( '--threads', '-t', default=1, type=int,
        help='The number of threads to use.')
        
    args = parser.parse_args()

    global NTHREADS
    NTHREADS = args.threads
    
    global WRITE_DEBUG_DATA
    WRITE_DEBUG_DATA = args.write_debug_data
    
    global VERBOSE
    VERBOSE = args.verbose
    global DEBUG_VERBOSE
    DEBUG_VERBOSE = args.debug_verbose
    
    ofp = ThreadSafeFile( args.ofname, "w" )
    
    return args.rnaseq_reads, args.reverse_rnaseq_strand, \
        args.cage_reads, args.rampage_reads, \
        args.polya_candidate_sites, ofp, args.reference

def main():
    rnaseq_bams, reverse_rnaseq_strand, cage_bams, rampage_bams, \
        polya_candidate_sites_fps, ofp, gtf_fname \
        = parse_arguments()

    ofp.write('track name="%s" visibility=2 itemRgb="On"\n' % ofp.name)
    
    rnaseq_reads = [ RNAseqReads(fp.name).init(reverse_read_strand=reverse_rnaseq_strand) 
                     for fp in rnaseq_bams ]
    
    cage_reads = [ CAGEReads(fp.name).init(reverse_read_strand=True) 
                   for fp in cage_bams ]

    rampage_reads = [ RAMPAGEReads(fp.name).init(reverse_read_strand=True) 
                      for fp in rampage_bams ]
    
    promoter_reads = [] + cage_reads + rampage_reads
    assert len(promoter_reads) > 0, "Must have either CAGE or RAMPAGE reads."
    
    if VERBOSE: print >> sys.stderr,  'Loading candidate polyA sites'
    polya_sites = find_polya_sites([x.name for x in polya_candidate_sites_fps])
    for fp in polya_candidate_sites_fps: fp.close()
    if VERBOSE: print >> sys.stderr, 'Finished loading candidate polyA sites'

    contig_lens = get_contigs_and_lens( rnaseq_reads, promoter_reads )
    
    gene_bndry_bins = defaultdict( lambda: None )
    if gtf_fname != None:
        gene_bndry_bins = load_gene_bndry_bins( gtf_fname, contig_lens )
    
    if VERBOSE: print >> sys.stderr, 'Finished loading gtf'
    
    for contig, contig_len in contig_lens.iteritems():
        for strand in '+-':
            find_exons_in_contig( (contig, strand, contig_len), ofp,
                                  rnaseq_reads, promoter_reads, polya_sites,
                                  gene_bndry_bins = gene_bndry_bins[(contig, strand)])
    
if __name__ == '__main__':
    main()



















#
