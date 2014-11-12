import os, sys
import cPickle

import math
import random

import numpy
import numpy.random
from scipy.special import gammaln, gamma, cbrt
import scipy.stats

from itertools import chain

try: import grit
except ImportError: sys.path.insert(0, "/home/nboley/grit/grit/")

import files.junctions
import grit.files.reads
from grit.lib.multiprocessing_utils import ProcessSafeOPStream

import pyximport; pyximport.install()
from call_peaks_support_fns import calc_moments

from scipy.optimize import fmin_l_bfgs_b as minimize

BACKGROUND_FRACTION = 0.01
MIN_NOISE_FRAC = 0.05
MIN_RD_CNT = 5
MIN_PEAK_SIZE = 5
MAX_PEAK_SIZE = 500
MIN_EMPTY_REGION_SIZE = 1

MAX_NUM_ITERATIONS = 25

VERBOSE = False
DEBUG_VERBOSE = False

SMOOTH_WIN_LEN = 10
PEAK_THRESHOLD = 1.0
SPLIT_TYPE = 'optimal'
#SPLIT_TYPE = 'random'
N_REPS = 1
if SPLIT_TYPE == 'random': assert N_REPS > 1

MIN_MERGE_SIZE = 10
MIN_REL_MERGE_SIZE = 0.5
TRIM_FRACTION = 0.01
MAX_EXP_FRACTION = 0.01

def write_bedgraph_from_array(array, region, ofprefix):
    """
    track name=CAGE.pan..plus type=bedGraph
    chr4    89932   89933   4.00
    chr4    89955   89956   2.00
    chr4    89958   89959   2.00
   """
    chrm = region['chrm']
    start = region['start']
    ofname = "%s.%s.bedgraph" % (
        ofprefix, {'+': 'plus', '-': 'minus'}[region['strand']])
    with open(ofname, 'w') as ofp:
        print >> ofp, "track name=%s type=bedGraph" % ofname
        for i, val in enumerate(array):
            if val < 1e-6: continue
            print >> ofp, "\t".join(
                ('chr' + chrm, str(start+i), str(start+i+1), "%.2f" % val))
    return

def write_bedgraph(chrm, peaks, ofp):
    """
    track name=CAGE.pan..plus type=bedGraph
    chr4    89932   89933   4.00
    chr4    89955   89956   2.00
    chr4    89958   89959   2.00
   """
    for start, stop, value in peaks:
        ofp.write( "\t".join(
                ('chr'+chrm, str(start), str(stop+1), "%.2f" % value)) + "\n")
    return

def build_false_signal(rnaseq_reads, signal_type):
    signal_type = '5p'
    assert signal_type in ('5p', '3p')
    # get the read start coverage
    signal_cov = numpy.zeros(region['stop']-region['start']+1, dtype=float)
    for rd1, rd2 in rnaseq_reads.iter_paired_reads(**region):
        if signal_type == '3p':
            pos = max(rd1.pos, rd1.aend, rd2.pos, rd2.aend)
        else:
            pos = min(rd1.pos, rd1.aend, rd2.pos, rd2.aend)
        if pos < region['start'] or pos > region['stop']: continue
        signal_cov[pos-region['start']] += 1
    
    n_rnaseq_reads = signal_cov.sum()
    
    # add the uniform background
    signal_cov = (1-BACKGROUND_FRACTION)*signal_cov+(
        n_rnaseq_reads*BACKGROUND_FRACTION)/len(signal_cov)

def build_control(rnaseq_reads, region, control_type, smooth_win_len=SMOOTH_WIN_LEN):
    assert control_type in ('5p', '3p')
    # get the read start coverage
    cov = numpy.zeros(region['stop']-region['start']+1, dtype=float)
    for rd1, rd2 in rnaseq_reads.iter_paired_reads(**region):
        if control_type == '3p':
            pos = max(rd1.pos, rd1.aend, rd2.pos, rd2.aend)
        else:
            pos = min(rd1.pos, rd1.aend, rd2.pos, rd2.aend)
        if pos < region['start'] or pos > region['stop']: continue
        cov[pos-region['start']] += 1
    
    n_rnaseq_reads = cov.sum()
    # add the uniform background
    cov = (1-BACKGROUND_FRACTION)*cov+(
        n_rnaseq_reads*BACKGROUND_FRACTION)/len(cov)
    
    # get the region segment boundaries
    region_tuple = (region['chrm'], region['strand'], region['start'], region['stop'])
    jns = files.junctions.load_junctions_in_bam(
        rnaseq_reads, [region_tuple,] )[(region['chrm'], region['strand'])]
    bndries = set((region['start']-region['start'], region['stop']-region['start']+1))
    for (start, stop), cnt, entropy in jns:
        bndries.add(start-region['start'])
        bndries.add(stop-region['start'])
    bndries = sorted(bndries)

    # smooth the signal in each segment
    min_signal = n_rnaseq_reads*BACKGROUND_FRACTION/len(cov)
    window = numpy.ones(smooth_win_len, dtype=float)/smooth_win_len
    for start, stop in zip(bndries[:-1], bndries[1:]):
        segment_signal = cov[start:stop]
        if stop - start <= smooth_win_len:
            cov[start:stop] = segment_signal.mean()
        else:    
            cov[start:stop] = numpy.convolve(
                window,segment_signal,mode='same')
    #cov[cov < min_signal] = min_signal
    return (cov + 1e-12)/(cov.sum() + 1e-12*len(cov))

def build_control_in_gene(gene, paired_rnaseq_reads, bndries, 
                          control_type, smooth_win_len=SMOOTH_WIN_LEN):
    assert control_type in ('5p', '3p')
    # get the read start coverage
    cov = numpy.zeros(gene.stop-gene.start+1, dtype=float)
    for rd_key, mappings in paired_rnaseq_reads:
        for mapping in mappings:
            poss = chain(chain(*mapping[4].cov_regions), 
                         chain(*mapping[4].cov_regions))
            if control_type == '3p':
                pos = max(poss)
            else:
                pos = min(poss)
            if pos < gene.start or pos > gene.stop: continue
            cov[pos-gene.start] += mapping[-1]
    
    n_rnaseq_reads = len(paired_rnaseq_reads)
    # add the uniform background
    cov = (1-BACKGROUND_FRACTION)*cov+(
        n_rnaseq_reads*BACKGROUND_FRACTION)/len(cov)
    
    # smooth the signal in each segment
    min_signal = n_rnaseq_reads*BACKGROUND_FRACTION/len(cov)
    window = numpy.ones(smooth_win_len, dtype=float)/smooth_win_len
    for start, stop in zip(bndries[:-1], bndries[1:]):
        segment_signal = cov[start-gene.start:stop-gene.start+1]
        if stop - start <= smooth_win_len:
            cov[start-gene.start:stop-gene.start+1] = segment_signal.mean()
        else:    
            cov[start-gene.start:stop-gene.start+1] = numpy.convolve(
                window,segment_signal,mode='same')
    #cov[cov < min_signal] = min_signal
    return (cov + 1e-12)/(cov.sum() + 1e-12*len(cov))


class TestSignificance(object):
    def __init__(self, signal_cov, control_cov, noise_frac):
        self.noise_n = int(noise_frac*sum(signal_cov)) + 1
        self.signal_n = sum(signal_cov)
        
        #### initialize the array that we will use to pick 
        #### the split base(s)
        self.split_statistic = signal_cov
        
        x = numpy.diff( numpy.asarray( 
                signal_cov >= 1e-6, dtype=int ) )
        stops = numpy.nonzero(x==1)[0].tolist()
        if signal_cov[-1] < 1e-6: stops.append(len(x))
        starts = (numpy.nonzero(x==-1)[0]+1).tolist()
        if signal_cov[0] < 1e-6: starts.insert(0, 0)
        self.zero_intervals = [ 
            (start, stop) for start, stop in zip(starts, stops)
            if stop - start + 1 >= MIN_EMPTY_REGION_SIZE ]
        
        #### initialize data to test for region significance
        # initialize the null data
        null_means = [0.,]
        null_vars = [0.,]
        for i, p in enumerate(control_cov):
            mean, var = calc_moments(p, self.noise_n)
            null_means.append(mean)
            null_vars.append(var)

        self.null_means_cumsum = numpy.array(null_means).cumsum()
        self.null_variances_cumsum = numpy.array(null_vars).cumsum()
        
        # initialize the signal test statistic
        lhds = ( signal_cov*numpy.log(control_cov)
                 - gammaln(1+signal_cov) )
        self.signal_lhd_cumsum = numpy.hstack((
            numpy.zeros(1), lhds.cumsum()))
        self.signal_cnts_cumsum = numpy.hstack((
            numpy.zeros(1), signal_cov.cumsum()))
    
    def __call__(self, start, stop, alpha):
        # if there are more reads in this region than noise reads, 
        # then this region must include some signal
        sig_cnt = ( 
            self.signal_cnts_cumsum[stop] 
            - self.signal_cnts_cumsum[start] )
        if sig_cnt > self.noise_n: return True
        
        mean = -(self.null_means_cumsum[stop] 
                 - self.null_means_cumsum[start])
        variance = ( self.null_variances_cumsum[stop] 
                     - self.null_variances_cumsum[start] )
        
        scale = variance/mean
        shape = mean/scale
        
        dist = scipy.stats.gamma(shape, scale=scale)
        critical_value = -scipy.stats.gamma(
            shape, scale=scale).isf(alpha)
        
        # calculate the value of the observed likelihood
        obs_lhd = ( self.signal_lhd_cumsum[stop] 
                    - self.signal_lhd_cumsum[start] )
        
        return obs_lhd < critical_value
    
    def find_split_bases(self, r_start, r_stop):
        """Returns a closed,open interval of bases to split. 

        """
        r_start += MIN_PEAK_SIZE
        r_stop -= MIN_PEAK_SIZE
        assert r_stop >= r_start
        if SPLIT_TYPE == 'random':
            rv = random.randint(r_start, r_stop)
            return rv, rv
        assert SPLIT_TYPE == 'optimal'

        # find the largest zero interval
        split_interval = None
        for start, stop in self.zero_intervals:
            if stop < r_start: continue
            if start > r_stop: break
            start = max(start, r_start)
            stop = min(stop, r_stop)
            if ( split_interval == None or
                 stop-start+1 > split_interval[1] - split_interval[0] ):
                split_interval = (start, stop)
        
        # if we found one, then use it. Otherwise, find the location with
        # the minimum signal
        if split_interval != None:
            #diff = split_interval[1] - split_interval[0]
            #return split_interval[0]+diff/2, split_interval[0]+diff/2
            return split_interval[0], split_interval[1]+1
        
        # find the bases that are the most below the mean
        min_val = self.split_statistic[r_start:r_stop+1].min()
        
        # find the indices of the minimum value
        min_indices = (
            self.split_statistic[r_start:r_stop+1] == min_val).nonzero()
            
        rv = random.choice(min_indices[0]) + r_start
        return rv, rv

def find_noise_regions(signal_cov, control_cov, noise_frac, alpha):
    alpha = alpha/(2*len(signal_cov))
    is_significant = TestSignificance(signal_cov, control_cov, noise_frac)
    noise_regions = []
    if signal_cov.sum() == 0:
        return [(0, len(signal_cov)),]
    # initialize the first region to split
    # trim 0 count bases from the edges of the signal track
    start, stop = 0, len(signal_cov)
    for i, cnt in enumerate(signal_cov): 
        if cnt > 0: break
        start = i
    if start > 0: noise_regions.append((0, start))
    for i in reversed(xrange(len(signal_cov))):
        if signal_cov[i] > 0: break
        stop = i
    if stop < len(signal_cov): noise_regions.append((stop,len(signal_cov)))
    regions_to_split = [((start, stop), 1)]
    
    # if the full region isn't significant, then we are done
    if not is_significant(*regions_to_split[0][0], alpha=alpha):
        return noise_regions + [regions_to_split[0][0],]
    while len(regions_to_split) > 0:
        # get the region to split - we know that this is significant
        # XXX use a better data structure
        (start, stop), level = regions_to_split.pop(0)
        # if this region is too small, then it's already significant
        # and so there is nothing to do 
        if stop - start < 2*MIN_PEAK_SIZE: continue
        
        # build the sub regions, and test them for significance
        left_bnd, right_bnd = is_significant.find_split_bases(start, stop)
        
        # add the split bases to the noise set
        if right_bnd > left_bnd:
            noise_regions.append((left_bnd, right_bnd))
        
        r1, r2 = [(start, left_bnd), (right_bnd, stop)]
        r1_sig, r2_sig = [
            is_significant(*r1, alpha=alpha), 
            is_significant(*r2, alpha=alpha) ]
        
        # if neither sub region is significant, (and we know the parent region 
        # was significant) then we are done
        if not r1_sig and not r2_sig:
            continue
                
        # add the subregions to the appropriate locations
        if r1_sig:
            regions_to_split.append((r1, level+1))
        else: noise_regions.append(r1)
            
        if r2_sig:
            regions_to_split.append((r2, level+1))
        else: noise_regions.append(r2)
    
    return sorted(noise_regions)

def estimate_noise_frac(noise_regions, signal_cov, control_cov):    
    noise_cnt = sum(signal_cov[start:stop].sum() 
                    for start, stop in noise_regions )
    control_cnt = sum(control_cov[start:stop].sum() 
                    for start, stop in noise_regions )
    assert control_cnt <= 1.0+1e-6
    expected_noise_cnt = (1./control_cnt)*noise_cnt
    signal_cnt = signal_cov.sum()
    # because this is a MOM estimate, it can lay out of the domain.
    # however, this should only occur in insignificant genes
    rv = min(1., expected_noise_cnt/(signal_cnt+1e-6))
    return max(MIN_NOISE_FRAC, rv)

def update_control_cov_for_five_prime_bias(
        noise_regions, noise_frac, 
        signal_cov, control_cov, reads_type):
    # disable the correction
    return (0.,1.), control_cov
    positions = []
    Ys = []
    ps = []
    max_pos = float(len(control_cov))
    for start, stop in sorted(noise_regions):
        for pos in xrange(start, stop):
            positions.append(pos)
            Ys.append(signal_cov[pos])
            ps.append(control_cov[pos])
    positions = numpy.array(positions, dtype=float)
    Ys = numpy.array(Ys, dtype=float)
    ps = numpy.array(ps, dtype=float)
    
    def calc_new_ps(args, positions, ps):
        alpha, power = args
        if reads_type == '5p':
            weights = (1 - positions/(max_pos+1))**power
        elif reads_type == '3p':
            weights = (positions/(max_pos+1))**power
        else:
            assert False
        
        new_ps = (weights/weights.mean())*alpha*ps + (1-alpha)*ps
        return new_ps/new_ps.sum()
    
    def calc_lhd_for_reg_coef(args):
        new_ps = calc_new_ps(args, positions, ps)
        res = -(Ys*numpy.log(new_ps)).sum()
        return res
    
    res = minimize(
        calc_lhd_for_reg_coef, x0=(0.1,1), 
        approx_grad=True, bounds=[(1e-6, 1-1e-6),(1,2)])
    reg_coef = res[0].tolist()
    return reg_coef, calc_new_ps(reg_coef, numpy.arange(max_pos), control_cov)

def merge_adjacent_intervals(
        intervals, max_abs_merge_distance, max_merge_fraction):
    if len(intervals) == 0: return []
    intervals.sort()
    merged_intervals = [list(intervals[0]),]
    prev_stop = merged_intervals[-1][1]
    for start, stop in intervals[1:]:
        max_merge_distance = max(
            max_abs_merge_distance, 
            max_merge_fraction*(stop-start),
            max_merge_fraction*(merged_intervals[-1][1]-merged_intervals[-1][0]))
        if ( start - max_merge_distance - 1 <= prev_stop
             and stop - start + 1 < MAX_PEAK_SIZE ):
            merged_intervals[-1][1] = stop
        else:
            merged_intervals.append([start, stop])
        prev_stop = stop
    return merged_intervals

def estimate_read_cov_and_call_peaks(
        region, signal_reads, reads_type, 
        rnaseq_reads, alpha=0.01):
    assert reads_type in ('promoter', 'polya')
    reads_type = '5p' if reads_type == 'promoter' else '3p'
    if region['strand'] == '-': 
        reads_type = {'3p':'5p', '5p':'3p'}[reads_type]
    
    signal_cov = signal_reads.build_read_coverage_array(**region)    
    if DEBUG_VERBOSE:
        config.log_statement("Finished building signal coverage array")
    #signal_cov = build_false_signal(rnaseq_reads, '5p')
    
    control_cov = build_control(
        rnaseq_reads, region, reads_type, SMOOTH_WIN_LEN)
    
    if DEBUG_VERBOSE:
        config.log_statement("Finished building control coverage array")

    return call_peaks(signal_cov, control_cov, reads_type, alpha)

def call_peaks( signal_cov, original_control_cov,
                reads_type, alpha=0.01):
    signal = numpy.ones(len(signal_cov))
    for k in xrange(N_REPS):
        noise_frac = 1.0
        noise_regions = [(0, len(signal)),]
        reg_coef, control_cov = \
            update_control_cov_for_five_prime_bias(
                noise_regions, noise_frac, 
                signal_cov, original_control_cov, reads_type)
        for i in xrange(MAX_NUM_ITERATIONS):
            if DEBUG_VERBOSE: 
                write_bedgraph_from_array(
                    1000*control_cov, region, "control.%i"%i)
                config.log_statement(
                    "Iter %i: Noise Frac %.2f%%\tReg Coef: %s" % (
                        i+1, noise_frac*100, reg_coef))
            noise_regions = find_noise_regions(
                signal_cov, control_cov, 
                noise_frac, alpha=alpha )
            new_noise_frac = estimate_noise_frac(
                noise_regions, signal_cov, control_cov)
            new_reg_coef, control_cov = \
                update_control_cov_for_five_prime_bias(
                    noise_regions, noise_frac, 
                    signal_cov, original_control_cov, reads_type)
            if noise_frac - new_noise_frac <= 1e-3 \
                    and abs(reg_coef[0] - new_reg_coef[0]) < 1e-3 \
                    and abs(reg_coef[1] - new_reg_coef[1]) < 1e-3: 
                break
            else: 
                noise_frac = new_noise_frac
                reg_coef = new_reg_coef
        
        for start, stop in noise_regions: 
            signal[start:stop] -= 1./N_REPS
    
    # build a list of inclusive peak starts and stops
    peaks = []
    nonzero_bases = (signal>1e-6).nonzero()[0].tolist()
    if len(nonzero_bases) == 0: return peaks
    curr_start = nonzero_bases.pop(0)
    curr_stop = curr_start
    for base in nonzero_bases:
        if base == curr_stop+1: 
            curr_stop += 1
        else:
            peaks.append((curr_start, curr_stop))
            curr_start, curr_stop = base, base
    
    peaks.append((curr_start, curr_stop))
    while True:
        new_peaks = merge_adjacent_intervals(
            peaks, MIN_MERGE_SIZE, MIN_REL_MERGE_SIZE)
        if len(new_peaks) == len(peaks):
            peaks = new_peaks
            break
        else:
            peaks = new_peaks
    
    new_peaks = []
    for start, stop in peaks:
        cov_region = signal_cov[start:stop+1]
        total_cov = cov_region.sum()
        cov_cumsum = cov_region.cumsum()-cov_region[0]
        trim_start = numpy.flatnonzero(
            cov_cumsum <= int(TRIM_FRACTION*total_cov)).max()
        try: trim_stop = numpy.flatnonzero(
                cov_cumsum >= int((1.0-TRIM_FRACTION)*total_cov)).min()
        except: trim_stop=stop-start+1
        new_peaks.append((trim_start+start, 
                          trim_stop+start,
                          cov_region[trim_start:trim_stop+1].sum()))

    exp_filtered_peaks = []
    max_peak_cnt = float(max(cnt for start, stop, cnt in new_peaks))
    for start, stop, cnt in new_peaks:
        length = stop - start + 1
        if (cnt >= MIN_RD_CNT
            and length >= MIN_PEAK_SIZE
            and length <= MAX_PEAK_SIZE
            and cnt/max_peak_cnt > MAX_EXP_FRACTION ): 
            exp_filtered_peaks.append((start, stop, cnt))
    return exp_filtered_peaks
