import numpy
numpy.seterr(all='ignore')

import sys, os
import math
import traceback
import time
def make_time_str(et):
    hours = et//3600
    mins = et//60 - 60*hours
    secs = et - 3660*hours - 60*mins
    return "%i:%i:%.4f" % ( hours, mins, secs )

import copy

import hashlib
def hash_array( array ):
    if None == array:
        return str(hash(None))
    b = array.view(numpy.uint8)
    return str(hashlib.sha1(b).hexdigest())

from scipy.linalg import svd, inv
from scipy.stats import chi2
from scipy.optimize import fminbound, brentq, bisect, line_search
from scipy.io import savemat

from StringIO import StringIO
from itertools import izip
from collections import defaultdict

import subprocess
import multiprocessing
from multiprocessing import Process

from math import sqrt

sys.path.append( os.path.join( os.path.dirname( __file__ ),
                               "..", "file_types", "fast_gtf_parser" ) )
from gtf import load_gtf, Transcript
sys.path.append( os.path.join( os.path.dirname( __file__ ),
                               "..", "file_types", "fast_wiggle_parser" ) )
from bedgraph import load_bedgraph
sys.path.append(os.path.join(os.path.dirname(__file__), "..", 'file_types'))
from wiggle import guess_strand_from_fname, fix_chr_name


sys.path.append( os.path.join(os.path.dirname( __file__ ), "..", "file_types" ))
                               
from reads import Reads, bin_reads

from f_matrix import calc_expected_cnts, find_nonoverlapping_boundaries, \
    build_nonoverlapping_indices, find_nonoverlapping_contig_indices

from frag_len import load_fl_dists, FlDist, build_normal_density

import cvxopt
from cvxopt import solvers, matrix, spdiag, log, div, sqrt

MIN_TRANSCRIPT_FREQ = 1e-12
# finite differences step size
FD_SS = 1e-8
COMPARE_TO_DCP = False
num_threads = 1
NUM_ITER_FOR_CONV = 5
DEBUG_OPTIMIZATION = False
PROMOTER_SIZE = 50
ABS_TOL = 1e-5
DEBUG = False

def log_warning(text):
    print >> sys.stderr, text
    return

class ThreadSafeFile( file ):
    def __init__( *args ):
        file.__init__( *args )
        args[0].lock = multiprocessing.Lock()

    def write( self, string ):
        self.lock.acquire()
        file.write( self, string )
        self.flush()
        self.lock.release()

def build_observed_cnts( binned_reads, fl_dists ):
    rv = {}
    for ( read_len, read_group, bin ), value in binned_reads.iteritems():
        rv[ ( read_len, fl_dists[read_group], bin ) ] = value
    
    return rv

def build_expected_and_observed_arrays( 
        expected_cnts, observed_cnts, normalize=True ):
    expected_mat = []
    observed_mat = []
    unobservable_transcripts = set()
    
    for key, val in expected_cnts.iteritems():
        # skip bins with 0 expected reads
        if sum( val) == 0:
            continue
        
        expected_mat.append( val )
        try:
            observed_mat.append( observed_cnts[key] )
        except KeyError:
            observed_mat.append( 0 )

    if len( expected_mat ) == 0:
        raise ValueError, "No expected reads."
    
    expected_mat = numpy.array( expected_mat, dtype=numpy.double )
    
    if normalize:
        nonzero_entries = expected_mat.sum(0).nonzero()[0]
        unobservable_transcripts = set(range(expected_mat.shape[1])) \
            - set(nonzero_entries.tolist())
        observed_mat = numpy.array( observed_mat, dtype=numpy.int )
        expected_mat = expected_mat[:,nonzero_entries]
        expected_mat = expected_mat/expected_mat.sum(0)
    
    return expected_mat, observed_mat, unobservable_transcripts

def nnls( X, Y, fixed_indices_and_values={} ):    
    X = matrix(X)
    Y = matrix(Y)
    
    m, n = X.size
    num_constraint = len( fixed_indices_and_values )
    
    G = matrix(0.0, (n,n))
    G[::n+1] = -1.0
    h = matrix(-MIN_TRANSCRIPT_FREQ, (n,1))

    # Add the equality constraints
    A=matrix(0., (1+num_constraint,n))
    b=matrix(0., (1+num_constraint,1))

    # Add the sum to one constraint
    A[0,:] = 1.
    b[0,0] = 1.
    
    # Add the fixed value constraints
    for const_i, (i, val) in enumerate(fixed_indices_and_values.iteritems()):
        A[const_i+1,i] = 1.
        b[const_i+1,0] = val
    
    solvers.options['show_progress'] = DEBUG_OPTIMIZATION
    res = solvers.qp(P=X.T*X, q=-X.T*Y, G=G, h=h, A=A, b=b)
    x = numpy.array(res['x']).T[0,]
    rss = ((numpy.array(X*res['x'] - Y)[0,])**2).sum()
    
    if DEBUG_OPTIMIZATION:
        for key, val in res.iteritems():
            if key in 'syxz': continue
            print "%s:\t%s" % ( key.ljust(22), val )
        
        print "RSS: ".ljust(22), rss
    
    return x

def solve_trust_region_sub_problem( gradient, hessian ):
    """
    minimize gradient*p + (1/2)p'*hessian*p
    
    st sum(p) == 1 and p > 0
    """
    from cvxpy import maximize, minimize, geq, eq, variable, matrix, program, \
        log, sum, quad_form, square, quad_over_lin, geo_mean, leq, norm2, \
        cvxpy_second_order_cone, belongs
    import cvxpy
    # belongs(x, cvxpy_second_order_cone())
    n = len( gradient )
    x = variable( n )
    # eq( sum(x), 1), geq( x, 1e-6 ), leq(norm2(x), 0.7)
    solvers.options['show_progress'] = False
    p = program( minimize( -matrix(gradient)*x + 0.5*quad_form(x, matrix(hessian))),
                 [ eq( sum(x), 1), geq( x, MIN_TRANSCRIPT_FREQ )] )
    # +1*quad_over_lin(1,thetas[1, 0]) ), 
    # p = program( minimize( -Xs*log(ps*thetas) ),
    #             [eq( sum(thetas), 1), geq( thetas, 0 )] )

    p.options['maxiters']  = 500
    p.options['show_progress']  = False
    #res = p.solve(quiet=not VERBOSE)
    res = p.solve(quiet=True)
    
    """
    G = matrix(0.0, (n,n))
    G[::n+1] = -1.0
    h = matrix(MIN_TRANSCRIPT_FREQ, (n,1))

    # Add the equality constraints
    A=matrix(0., (1,n))
    b=matrix(0., (1,1))

    # Add the sum to one constraint
    A[0,:] = 1.
    b[0,0] = 1.
    
    solvers.options['show_progress'] = DEBUG_OPTIMIZATION
    res = solvers.qp(P=-matrix(hessian), q=-matrix(gradient), G=-G, h=-h, A=A, b=b )
    #res = solvers.qp(P=-matrix(hessian), q=-matrix(gradient), G=G, h=h, A=A, b=b)
    print res
    return res
    x = numpy.array(res['x']).T[0,]
    rss = ((numpy.array(X*res['x'] - Y)[0,])**2).sum()
    """
    return numpy.array( x.value.T )[0,:]


def build_expected_and_observed_rnaseq_counts( gene, bam_fname, fl_dists ):
    # load the bam file
    reads = Reads(bam_fname)
    
    # find the set of non-overlapping exons, and convert the transcripts to 
    # lists of these non-overlapping indices. All of the f_matrix code uses
    # this representation.     
    exon_boundaries = find_nonoverlapping_boundaries(gene.transcripts)
    transcripts_non_overlapping_exon_indices = \
        list(build_nonoverlapping_indices( 
                gene.transcripts, exon_boundaries ))
    
    binned_reads = bin_reads( 
        reads, gene.chrm, gene.strand, exon_boundaries, False, True)
    
    observed_cnts = build_observed_cnts( binned_reads, fl_dists )    
    read_groups_and_read_lens =  { (RG, read_len) for RG, read_len, bin 
                                   in binned_reads.iterkeys() }
    
    fl_dists_and_read_lens = [ (fl_dists[RG], read_len) for read_len, RG  
                               in read_groups_and_read_lens ]
    
    expected_cnts = calc_expected_cnts( 
        exon_boundaries, transcripts_non_overlapping_exon_indices, 
        fl_dists_and_read_lens)
    
    reads.close()
    
    return expected_cnts, observed_cnts

def build_expected_and_observed_promoter_counts( gene, cage_array ):
    # find the promoters
    promoters = list()
    for transcript in gene.transcripts:
        if transcript.strand == '+':
            promoter = [ transcript.exons[0][0], 
                         min( transcript.exons[0][0] + PROMOTER_SIZE, 
                              transcript.exons[0][1]) ]
        else:
            assert transcript.strand == '-'
            promoter = [ max( transcript.exons[-1][1] - PROMOTER_SIZE, 
                              transcript.exons[-1][0] ),
                         transcript.exons[-1][1] ]
        promoters.append( tuple( promoter ) )
    
    promoter_boundaries = set()
    for start, stop in promoters:
        promoter_boundaries.add( start )
        promoter_boundaries.add( stop + 1 )
    promoter_boundaries = numpy.array( sorted( promoter_boundaries ) )
    pseudo_promoters = zip(promoter_boundaries[:-1], promoter_boundaries[1:])
    
    # build the design matrix. XXX FIXME
    expected_cnts = defaultdict( lambda: [0.]*len(gene.transcripts) )
    for transcript_i, promoter in enumerate(promoters):
        nonoverlapping_indices = \
            find_nonoverlapping_contig_indices( 
                [promoter,], promoter_boundaries )
        # calculate the count probabilities, adding a fudge to deal with 0
        # frequency bins
        tag_cnt = cage_array[
            promoter[0]-gene.start:promoter[1]-gene.start].sum() \
            + 1e-6*len(nonoverlapping_indices)
        for i in nonoverlapping_indices:
            ps_promoter = pseudo_promoters[i]
            ps_tag_cnt = cage_array[
                ps_promoter[0]-gene.start:ps_promoter[1]-gene.start].sum()
            expected_cnts[ ps_promoter ][transcript_i] \
                = (ps_tag_cnt+1e-6)/tag_cnt
        
    # count the reads in each non-overlaping promoter
    observed_cnts = {}
    for (start, stop) in expected_cnts.keys():
        observed_cnts[ (start, stop) ] \
            = int(round(cage_array[start-gene.start:stop-gene.start].sum()))
    
    return expected_cnts, observed_cnts

def build_design_matrices( gene, bam_fname, fl_dists, cage_array ):
    import cPickle
    # only do this if we are in debugging mode
    if num_threads == 1:
        obj_name = "." + gene.id + "_" + hash_array(cage_array) \
            + "_" + os.path.basename(bam_fname) + ".obj"
        try:
            fp = open( obj_name )
            expected_array, observed_array, unobservable_transcripts = \
                cPickle.load( fp )
            fp.close()
            return expected_array, observed_array, unobservable_transcripts
        except IOError:
            if DEBUG_VERBOSE: 
                print "Couldnt load cached"
    
    # bin the rnaseq reads
    expected_rnaseq_cnts, observed_rnaseq_cnts = \
        build_expected_and_observed_rnaseq_counts( gene, bam_fname, fl_dists )
    expected_rnaseq_array, observed_rnaseq_array, unobservable_rnaseq_trans = \
        build_expected_and_observed_arrays( 
            expected_rnaseq_cnts, observed_rnaseq_cnts, True )

    if cage_array == None:
        return expected_rnaseq_array, \
            observed_rnaseq_array, \
            unobservable_rnaseq_trans
    
    # bin the CAGE data
    expected_promoter_cnts, observed_promoter_cnts = \
        build_expected_and_observed_promoter_counts( gene, cage_array )
    expected_prom_array, observed_prom_array, unobservable_prom_trans = \
        build_expected_and_observed_arrays( 
        expected_promoter_cnts, observed_promoter_cnts, False )
    
    # combine the arrays
    observed_rnaseq_array = numpy.delete( observed_rnaseq_array, 
                  numpy.array(list(unobservable_prom_trans)) )
    observed_prom_array = numpy.delete( observed_prom_array, 
                  numpy.array(list(unobservable_rnaseq_trans)) )
    observed_array = numpy.hstack((observed_prom_array, observed_rnaseq_array))
    if observed_array.sum() == 0:
        raise ValueError, "TOO FEW READS"

    expected_rnaseq_array = numpy.delete( expected_rnaseq_array, 
                  numpy.array(list(unobservable_prom_trans)), axis=1 )
    expected_prom_array = numpy.delete( expected_prom_array, 
                  numpy.array(list(unobservable_rnaseq_trans)), axis=1 )   
    
    expected_array = numpy.vstack((expected_prom_array, expected_rnaseq_array))
    unobservable_transcripts \
        = unobservable_rnaseq_trans.union(unobservable_prom_trans)

    if num_threads == 1:
        fp = open( obj_name, "w" )
        cPickle.dump( (expected_array,observed_array,unobservable_transcripts),
                      fp )
        fp.close()

    return expected_array, observed_array, unobservable_transcripts

try:
    from sparsify_support_fns import calc_lhd, calc_gradient, calc_hessian
except ImportError:
    raise
    def calc_lhd( freqs, observed_array, expected_array ):
        return float(observed_array*numpy.log( 
                numpy.matrix( expected_array )*numpy.matrix(freqs).T ))

    def calc_lhd_deriv( freqs, observed_array, expected_array ):
        denom = numpy.matrix( expected_array )*numpy.matrix(freqs).T
        rv = (((expected_array.T)*observed_array))*(1.0/denom)
        return -numpy.array(rv)[:,0]

def calc_taylor_lhd( freqs, observed_array, expected_array ):
    observed_prbs = (observed_array + MIN_TRANSCRIPT_FREQ)/observed_array.sum()
    bin_prbs = expected_array.dot(freqs)
    fo_taylor_penalty = (observed_prbs - bin_prbs)/observed_prbs
    taylor_penalty = numpy.log(observed_prbs) + fo_taylor_penalty \
        - numpy.power(fo_taylor_penalty, 2)/2
    return float( observed_array.dot( taylor_penalty ).sum() )

def project_onto_simplex( x, debug=False ):
    if ( x >= MIN_TRANSCRIPT_FREQ ).all() and abs( 1-x.sum()  ) < 1e-6: return x
    sorted_x = numpy.sort(x)[::-1]
    if debug: print "sorted x:", sorted_x
    n = len(sorted_x)
    if debug: print "cumsum:", sorted_x.cumsum()
    if debug: print "arange:", numpy.arange(1,n+1)
    rhos = sorted_x - (1./numpy.arange(1,n+1))*( sorted_x.cumsum() - 1 )
    if debug: print "rhos:", rhos
    rho = (rhos > 0).nonzero()[0].max() + 1
    if debug: print "rho:", rho
    theta = (1./rho)*( sorted_x[:rho].sum()-1)
    if debug: print "theta:", theta
    x_minus_theta = x - theta
    if debug: print "x - theta:", x_minus_theta
    x_minus_theta[ x_minus_theta < 0 ] = MIN_TRANSCRIPT_FREQ
    return x_minus_theta

def estimate_transcript_frequencies_line_search(  
        observed_array, full_expected_array, x0, dont_zero, abs_tol ):
    expected_array = full_expected_array.copy()
    def f_lhd(x):
        log_lhd = calc_lhd(x, observed_array, expected_array)
        return log_lhd
    
    def f_gradient(x):
        return calc_gradient( x, observed_array, expected_array )
        
    def calc_max_feasible_step_size_and_limiting_index( x0, gradient ):
        """Calculate the maximum step size to stay in the feasible region.
        
        solve y - x*gradient = MIN_TRANSCRIPT_FREQ for x
        x = (y - MIN_TRANSCRIPT_FREQ)/gradient
        """
        # we use minus because we return a positive step
        steps = (x0-MIN_TRANSCRIPT_FREQ)/gradient
        step_size = -steps[ steps < 0 ].max()
        step_size_i = ( steps == -step_size ).nonzero()[0]
        return step_size, step_size_i
    
    def calc_projected_gradient( x ):
        gradient = f_gradient( x )
        gradient = gradient/gradient.sum()
        x_next = project_onto_simplex( x + 1.*gradient )
        gradient = (x_next - x)
        return gradient

    def maximum_step_is_optimal( x, gradient, max_feasible_step_size ):
        """Check the derivative at the maximum step to determine whether or 
           not the maximum step is a maximum along the gradient line.

        """
        max_feasible_step_size, max_index = \
            calc_max_feasible_step_size_and_limiting_index(x, gradient)
        if max_feasible_step_size > FD_SS and \
                f_lhd( x + (max_feasible_step_size-FD_SS)*gradient ) \
                > f_lhd( x + max_feasible_step_size*gradient ):
            return False
        else:
            return True

    def line_search( x, gradient, max_feasible_step_size ):
        def brentq_fmin(alpha):
            return f_lhd(x + (alpha+FD_SS)*gradient) \
                - f_lhd(x + (alpha-FD_SS)*gradient)
        
        def downhill_search(step_size):
            step_size = min_step_size
            curr_lhd = f_lhd( x )
            while step_size > FD_SS and curr_lhd > f_lhd( x+step_size*gradient ):
                step_size /= 1.5
            return int(step_size> FD_SS)*step_size
        
        min_step_size = FD_SS
        max_step_size = max_feasible_step_size-FD_SS
        if brentq_fmin(max_step_size) >= 0:
            return max_step_size, False
        elif brentq_fmin(min_step_size) <= 0:
            step_size = downhill_search(min_step_size)
            return step_size, True

        # do a line search with brent
        step_size = brentq(brentq_fmin, min_step_size, max_step_size )
        if f_lhd(x) > f_lhd(x+step_size*gradient):
            step_size = downhill_search( step_size )
            return step_size, (step_size==0)

        #print "Found step size:", step_size, "vs", max_feasible_step_size, \
        #    "Graniness is %i" % graniness
        return step_size, True
    
    n = full_expected_array.shape[1]
    x = x0.copy()
    prev_lhd = 1e-10
    lhd = f_lhd(x)
    lhds = []
    zeros = set()
    zeros_counter = 0
    for i in xrange( 500 ):
        # calculate the gradient and the maximum feasible step size
        gradient = calc_projected_gradient( x )
        gradient /= numpy.absolute(gradient).sum()
        max_feasible_step_size, max_index = \
            calc_max_feasible_step_size_and_limiting_index(x, gradient)
        
        # perform the line search
        alpha, is_full_step = line_search(
            x, gradient, max_feasible_step_size)
        x += alpha*gradient
        """
        if alpha > 0 and not is_full_step:
            print "%i\t%.2f\t%.6e\t%i" % ( 
                i, f_lhd(x), f_lhd(x)-f_lhd(x-alpha*gradient), len(x) )
            continue
        """
        
        if abs( 1-x.sum() ) > 1e-6:
            x = project_onto_simplex(x)
            continue
     
        if i > 30 and (alpha == 0 or f_lhd(x) - prev_lhd < abs_tol):
            zeros_counter += 1
            if zeros_counter > 3:
                break            
        else:
            zeros_counter = 0
            if not dont_zero:
                current_nonzero_entries = (x > 1e-12).nonzero()[0]
                if len( current_nonzero_entries ) < len(x):
                    n = full_expected_array.shape[1]
                    full_x = numpy.ones(n)*MIN_TRANSCRIPT_FREQ
                    full_x[ numpy.array(sorted(set(range(n))-zeros)) ] = x
                    
                    zeros = set( (full_x <= 1e-12).nonzero()[0] )
                    # build the x set
                    x = x[ current_nonzero_entries ]
                    expected_array = expected_array[:, current_nonzero_entries]
            
        
        prev_lhd = lhd
        lhd = f_lhd(x)
        lhds.append( lhd )
        if DEBUG_OPTIMIZATION:
            print "%i\t%.2f\t%.6e\t%i" % ( 
                i, lhd, lhd - prev_lhd, len(x) )
    
    final_x = numpy.ones(n)*MIN_TRANSCRIPT_FREQ
    final_x[ numpy.array(sorted(set(range(n))-zeros)) ] = x
    final_lhd = calc_lhd(final_x, observed_array, full_expected_array)
    assert final_lhd >= f_lhd(x) - abs_tol
    return final_x, lhds

def estimate_transcript_frequencies(  
        observed_array, full_expected_array, abs_tol ):
    fp = open( "lhd_change.txt", "w" )
    if observed_array.sum() == 0:
        raise ValueError, "Too few reads."
    n = full_expected_array.shape[1]
    if n == 1:
        return numpy.ones( 1, dtype=float )
    
    x = numpy.array([1./n]*n)
    #x = project_onto_simplex( nnls( full_expected_array, observed_array ) )
    eps = 10.
    start_time = time.time()
    if DEBUG_VERBOSE:
        print "Iteration\tlog lhd\t\tchange lhd\tn iter\ttolerance\ttime (hr:min:sec)"
    for i in xrange( 500 ):
        prev_x = x.copy()
        
        x, lhds = estimate_transcript_frequencies_line_search(  
            observed_array, full_expected_array, x, 
            dont_zero=False, abs_tol=eps )
        for lhd in lhds: fp.write( "%e\n" % lhd )
        fp.flush()
        
        lhd = calc_lhd( x, observed_array, full_expected_array )
        prev_lhd = calc_lhd( prev_x, observed_array, full_expected_array )
        if DEBUG_VERBOSE:
            print "Zeroing %i\t%.2f\t%.2e\t%i\t%e\t%s" % ( 
                i, lhd, (lhd - prev_lhd)/len(lhds), len(lhds ), eps, 
                make_time_str((time.time()-start_time)/len(lhds)) )
            
        start_time = time.time()
        
        if float(lhd - prev_lhd)/len(lhds) < eps:
            if eps == abs_tol: break
            eps /= 5
            eps = max( eps, abs_tol )
        
    
    for i in xrange( 5 ):
        prev_x = x.copy()
        x, lhds = estimate_transcript_frequencies_line_search(  
            observed_array, full_expected_array, x, 
            dont_zero=True, abs_tol=abs_tol)
        for lhd in lhds: fp.write( "%e\n" % (lhd - prev_lhd) )
        fp.flush()
        lhd = calc_lhd( x, observed_array, full_expected_array )
        prev_lhd = calc_lhd( prev_x, observed_array, full_expected_array )
        if DEBUG_VERBOSE:
            print "Non-Zeroing %i\t%.2f\t%.2e\t%i\t%e\t%s" % ( 
                i, lhd, (lhd - prev_lhd)/len(lhds), len(lhds), eps,
                make_time_str((time.time()-start_time)/len(lhds)))
        
        start_time = time.time()
        if len( lhds ) < 500: break

    fp.close()    
    return x

def estimate_confidence_bound_by_bisection( 
        observed_array, expected_array,
        fixed_i, optimal_est, 
        bound_type, alpha ):    
    n_transcripts = expected_array.shape[1]
    max_test_stat = chi2.ppf( 1 - alpha, 1 )/2.
            
    def calc_test_statistic(x):
        return calc_lhd( x, observed_array, expected_array )
        
    def etf_wrapped( x ):
        constrained_est = estimate_transcript_frequencies( 
            observed_array, expected_array, 
            [fixed_i,], fixed_values=[x,] )
        return constrained_est
    
    def estimate_bound( ):
        optimum_bnd = optimal_est[fixed_i]
        if bound_type=='UPPER': 
            other_bnd = 1.0 - n_transcripts*MIN_TRANSCRIPT_FREQ
        else: 
            other_bnd = MIN_TRANSCRIPT_FREQ
        
        # check to see if the far bound is sufficiently bad ( this 
        # is really an identifiability check 
        test_stat = calc_test_statistic(etf_wrapped(other_bnd))
        if max_test_stat - test_stat > 0:
            print "TERMINATE", optimal_RSS, calc_RSS(etf_wrapped(other_bnd), observed_array, expected_array), \
                calc_RSS(etf_wrapped(other_bnd), observed_array, expected_array) - optimal_RSS
            return other_bnd
        
        lower_bnd, upper_bnd = \
            min( optimum_bnd, other_bnd ), max( optimum_bnd, other_bnd )
        def obj( x ):
            rv = max_test_stat - calc_test_statistic(etf_wrapped(x))
            return rv
        
        return brentq( obj, lower_bnd, upper_bnd)
                       
    
    bnd = estimate_bound()
    test_stat = calc_test_statistic(etf_wrapped(bnd))
    
    return test_stat, bnd

class MaxIterError( ValueError ):
    pass

def estimate_confidence_bound( observed_array, 
                               expected_array, 
                               fixed_index,
                               mle_estimate,
                               bound_type,
                               alpha = 0.025):
    try:
        raise ValueError, "TEST"
        return estimate_confidence_bounds_directly( 
            observed_array,  expected_array, fixed_index,
            mle_estimate, bound_type, alpha )
    except ValueError:
        #print "WARNING: Couldn't find bound for transcript %i %s directly: trying bisection." % (fixed_index+1, bound_type)
        return estimate_confidence_bound_by_bisection( 
            observed_array,  expected_array, fixed_index,
            mle_estimate, bound_type, alpha )
                    
    pass

def build_mle_estimate( observed_array, expected_array ):
    log_lhd, mle_estimate = estimate_transcript_frequencies( 
        observed_array, expected_array )

def build_bound(observed_array, expected_array, mle_estimate, index, bnd_type):
    p_value, bnd = estimate_confidence_bound( 
        observed_array, expected_array, index, mle_estimate, bnd_type )

def write_gene_to_gtf( ofp, gene, mles, lbs=None, ubs=None, 
                       abs_filter_value=0, rel_filter_value=0 ):
    scores = [ int(1000*mle/max( mles )) for mle in mles ]
    
    for index, (mle, transcript) in enumerate(izip(mles, gene.transcripts)):
        if mle <= abs_filter_value: continue
        if mle/max( mles ) <= rel_filter_value: continue
        meta_data = { "frac": "%.2e" % mle }
        transcript.score = max( 1, int(1000*mle/max(mles)) )
        if lbs != None:
            meta_data["conf_lo"] = "%.2e" % lbs[index]
        if ubs != None:
            meta_data["conf_hi"] = "%.2e" % ubs[index]
        
        ofp.write( transcript.build_gtf_lines(
                gene.id, meta_data, source="grit") + "\n" )
    
    return

def estimate_gene_expression_worker( ip_lock, input_queue, 
                                     op_lock, output, 
                                     estimate_confidence_bounds, 
                                     abs_filter_value,
                                     rel_filter_value):
    while True:
        op_lock.acquire()
        if 0 == len( output ):
            op_lock.release()
            return
        op_lock.release()
        
        ip_lock.acquire()
        try:
            work_type, ( gene_id, bam_fn, trans_index ) = input_queue.pop()
        except IndexError:
            ip_lock.release()
            if num_threads == 1:
                return
            time.sleep(2.0)
            continue
        ip_lock.release()
        
        if work_type == 'design_matrices':
            op_lock.acquire()
            gene = output[(gene_id, bam_fn)]['gene']
            fl_dists = output[(gene_id, bam_fn)]['fl_dists']
            cage = output[(gene_id, bam_fn)]['cage']
            op_lock.release()
                        
            try:
                expected_array, observed_array, unobservable_transcripts \
                    = build_design_matrices( gene, bam_fn, fl_dists, cage )
            except ValueError, inst:
                print "Skipping %s: %s" % ( gene.id, inst )
                if DEBUG: raise
                continue
            
            if VERBOSE: print "FINISHED DESIGN MATRICES %s\t%s" % ( 
                gene_id, bam_fn )
            
            op_lock.acquire()
            gene_data = output[(gene_id, bam_fn)]
            gene_data['design_matrices'] = \
                (observed_array, expected_array, unobservable_transcripts)
            output[(gene_id, bam_fn)] = gene_data
            op_lock.release()

            ip_lock.acquire()
            input_queue.append( ('mle', (gene.id, bam_fn, None)) )
            ip_lock.release()
        elif work_type == 'mle':
            op_lock.acquire()
            observed_array, expected_array, unobservable_transcripts = \
                output[(gene_id, bam_fn)]['design_matrices']
            op_lock.release()

            try:
                mle_estimate = estimate_transcript_frequencies( 
                    observed_array, expected_array, ABS_TOL)
            except ValueError, inst:
                print "Skipping %s: %s" % ( gene_id, inst )
                if DEBUG: raise
                continue
            
            log_lhd = calc_lhd( mle_estimate, observed_array, expected_array)
            if VERBOSE: print "FINISHED MLE %s\t%s\t%.2f\n%s" % ( 
                gene_id, bam_fn, log_lhd, 
                ["%.2e" % x for x in mle_estimate] )
            
            op_lock.acquire()
            gene_data = output[(gene_id, bam_fn)]
            gene_data['mle'] = mle_estimate
            output[(gene_id, bam_fn)] = gene_data
            op_lock.release()
            
            if estimate_confidence_bounds:
                ip_lock.acquire()
                for i in xrange(expected_array.shape[1]):
                    input_queue.append( ('lb', (gene_id, bam_fn, i)) )
                    input_queue.append( ('ub', (gene_id, bam_fn, i)) )
                ip_lock.release()
        
        elif work_type in ('lb', 'ub'):
            op_lock.acquire()
            observed_array, expected_array, unobservable_transcripts = \
                output[(gene_id, bam_fn)]['design_matrices']
            mle_estimate = output[(gene_id, bam_fn)]['mle']
            op_lock.release()
            
            bnd_type = 'LOWER' if work_type == 'lb' else 'UPPER'
            
            p_value, bnd = estimate_confidence_bound( 
                observed_array, expected_array, 
                trans_index, mle_estimate, bnd_type )
            if VERBOSE: print "FINISHED %s BOUND %s\t%s\t%i\t%.2e\t%.2e" % ( 
                bnd_type, gene_id, bam_fn, trans_index, bnd, p_value )
            
            op_lock.acquire()
            gene_data = output[(gene_id, bam_fn)]
            bnds = gene_data[work_type + 's']
            bnds[trans_index] = bnd
            gene_data[work_type + 's'] = bnds
            output[(gene_id, bam_fn)] = gene_data
            op_lock.release()

class WorkQueue():
    def __init__(self):
        pass

def gene_is_finished( gene_data, compute_confidence_bounds ):
    if gene_data['mle'] == None:
        return False
    
    if not compute_confidence_bounds:
        return True
    
    observed_array, expected_array, unobservable_transcripts \
        = gene_data['design_matrices']
    
    n_transcripts = gene_data['design_matrices'][1].shape[1]
    assert n_transcripts == len( gene_data['gene'].transcripts ) \
        - len( gene_data['design_matrices'][1] )
    if n_transcripts != len(gene_data['ubs']):
        return False
    if n_transcripts != len(gene_data['lbs']):
        return False
    
    return True

def write_finished_data_to_disk( output_dict, output_dict_lock, ofps,
                                 compute_confidence_bounds=True, 
                                 write_design_matrices=False,
                                 abs_filter_value=0.0,
                                 rel_filter_value=0.0 ):
    written_design_matrices = set()
    while True:
        output_dict_lock.acquire()
        if len( output_dict ) == 0: 
            output_dict_lock.release()
            return
        output_dict_lock.release()
        
        for key in output_dict.keys():
            # write out the design matrix
            if write_design_matrices \
                    and key not in written_design_matrices \
                    and None != output_dict[key]['design_matrices']:
                observed, expected, missed = output_dict[key]['design_matrices']
                ofname = "./%s_%s.mat" % ( key[0], os.path.basename(key[1]) )
                if DEBUG_VERBOSE: print "Writing mat to '%s'" % ofname
                savemat( ofname, {'observed': observed, 'expected': expected} )
                if DEBUG_VERBOSE: print "Finished writing mat to '%s'" % ofname
                          
                written_design_matrices.add( key )
                
            if gene_is_finished( output_dict[key], compute_confidence_bounds ):
                output_dict_lock.acquire()
                assert gene_is_finished( 
                    output_dict[key], compute_confidence_bounds )
                if VERBOSE: print "Finished processing", key

                gene = output_dict[key]['gene']
                mles = output_dict[key]['mle']
                lbs = output_dict[key]['lbs'] if compute_confidence_bounds else None
                ubs = output_dict[key]['ubs'] if compute_confidence_bounds else None

                write_gene_to_gtf( ofps[key[1]], gene, mles, 
                                   lbs, ubs, abs_filter_value, rel_filter_value)
                
                del output_dict[key]
                output_dict_lock.release()
        
        time.sleep(2.0)
    
    return

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='Determine valid transcripts and estimate frequencies.')
    parser.add_argument( 'ofname', help='Output filename.')
    parser.add_argument( 'gtf', type=file,
        help='GTF file processed for expression')

    parser.add_argument( 'bam_fn', metavar='bam', type=file,
        help='list of bam files to for which to produce expression')    
    
    parser.add_argument( '--fl-dists', type=file, nargs='+', 
       help='a pickled fl_dist object(default:generate fl_dist from input bam)')
    parser.add_argument( '--fl-dist-norm', \
        help='mean and standard deviation (format "mn:sd") from which to ' \
            +'produce a fl_dist_norm (default:generate fl_dist from input bam)')

    parser.add_argument( '--cage-fns', nargs='+', type=file,
        help='list of bedgraph files with CAGE tag counts.')
    
    parser.add_argument( '--threads', '-t', type=int , default=1,
        help='Number of threads spawn for multithreading (default=1)')

    parser.add_argument( '--estimate-confidence-bounds', '-c', default=False,
        action="store_true",
        help='Whether or not to calculate confidence bounds ( this is slow )')

    parser.add_argument( '--write-design-matrices', default=False,
        action="store_true",
        help='Write the design matrices out to a matlab-style matrix file.')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',
                             help='Whether or not to print status information.')
    parser.add_argument( '--debug-verbose', default=False, action='store_true',
                             help='Prints the optimization path updates.')

    args = parser.parse_args()

    if not args.fl_dists and not args.fl_dist_norm:
        raise ValueError, "Must specific either --fl-dists or --fl-dist-norm."

    if args.fl_dist_norm != None:
        try:
            mean, sd = args.fl_dist_norm.split(':')
            mean = int(mean)
            sd = int(sd)
            fl_dist_norm = (mean, sd)
        except ValueError:
            raise ValueError, "Mean and SD for normal fl_dist are not properly formatted. Expected '--fl-dist-norm MEAN:SD'."
        
        mean, sd = fl_dist_norm
        fl_min = max( 0, mean - (4 * sd) )
        fl_max = mean + (4 * sd)
        fl_dists = { 'mean': build_normal_density( fl_min, fl_max, mean, sd ) }
        read_group_mappings = []
    else:
        fl_dists, read_group_mappings = load_fl_dists( 
            fp.name for fp in args.fl_dists )
    
    global DEBUG_VERBOSE
    DEBUG_VERBOSE = args.debug_verbose
    
    global VERBOSE
    VERBOSE = ( args.verbose or DEBUG_VERBOSE )

    # we change to the output directory later, and these files need to opened in
    # each sub-process for thread safety, so we get the absokute path while we 
    # can. We allow for multiple bams at this stage, but insist upon one in the
    # arguments
    bam_fns = [ os.path.abspath( args.bam_fn.name ),  ]
    
    global PROCESS_SEQUENTIALLY
    if args.threads == 1:
        PROCESS_SEQUENTIALLY = True
    
    global num_threads
    num_threads = args.threads
    
    ofps = {}
    for bam_fn in bam_fns:
        ofps[bam_fn] = ThreadSafeFile( args.ofname, "w" )
        ofps[bam_fn].write(
            "track name=transcripts.%s useScore=1\n" % os.path.basename(bam_fn))
    
    cage_fns = [] if args.cage_fns == None else [ 
        fp.name for fp in args.cage_fns ]    
    
    return args.gtf, bam_fns, cage_fns, ofps, fl_dists, read_group_mappings, \
        args.estimate_confidence_bounds, args.write_design_matrices
    
def main():
    # Get file objects from command line
    gtf_fp, bam_fns, cage_fns, ofps, fl_dists, rg_mappings, \
        estimate_confidence_bounds, write_design_matrices = parse_arguments()
    abs_filter_value = 1e-12 + 1e-16
    rel_filter_value = 0
    
    manager = multiprocessing.Manager()

    input_queue_lock = manager.Lock()
    input_queue = manager.list()

    output_dict_lock = manager.Lock()    
    output_dict = manager.dict()

    cage = {}
    for cage_fn in cage_fns:
        strand = guess_strand_from_fname( cage_fn )
        data = load_bedgraph( cage_fn )
        for chrm, array in data.iteritems():
            if cage.has_key((strand, fix_chr_name(chrm))):
                raise ValueError, "Duplicated keys for promoter data."
            cage[(strand, fix_chr_name(chrm))] = numpy.array(array)
    
    if VERBOSE:
        print "Finished Loading CAGE"
    
    # add all the genes, in order of longest first. 
    genes = load_gtf( gtf_fp.name )
    if VERBOSE:
        print "Finished Loading %s" % gtf_fp.name

    # add the genes in reverse sorted order so that the longer genes are dealt
    # with first
    for gene in sorted( genes, key=lambda x: len(x.transcripts) ):
        if (gene.strand, gene.chrm) not in cage:
            print "WARNING: Could not find CAGE signal for strand %s chr %s"% (
                gene.strand, gene.chrm )
            gene_cage = None
        else:
            gene_cage = cage[(gene.strand, gene.chrm)][
                gene.start:gene.stop+1].copy()
        
        for bam_fn in bam_fns:
            input_queue.append(('design_matrices', (gene.id, bam_fn, None)))
            gene_data = manager.dict()
            gene_data[ 'gene' ] = gene
            gene_data[ 'cage' ] = gene_cage
            gene_data[ 'fl_dists' ] = fl_dists
            gene_data[ 'lbs' ] = {}
            gene_data[ 'ubs' ] = {}
            gene_data[ 'mle' ] = None
            gene_data[ 'design_matrices' ] = None
            output_dict[ (gene.id, bam_fn) ] = gene_data

    del cage
    
    if 1 == num_threads:
        estimate_gene_expression_worker( 
            input_queue_lock, input_queue, 
            output_dict_lock, output_dict, 
            estimate_confidence_bounds,
            abs_filter_value, rel_filter_value )
        write_finished_data_to_disk( 
            output_dict, output_dict_lock, ofps,
            estimate_confidence_bounds, 
            write_design_matrices,
            abs_filter_value, rel_filter_value )
    else:
        ps = []
        for i in xrange( num_threads ):
            p = Process(target=estimate_gene_expression_worker, 
                        args=( input_queue_lock, input_queue, 
                               output_dict_lock, output_dict, 
                               estimate_confidence_bounds,
                               abs_filter_value, rel_filter_value ) )
            p.start()
            ps.append( p )    

        write_p = Process(target=write_finished_data_to_disk, args=(
                output_dict, output_dict_lock, ofps,
                estimate_confidence_bounds, write_design_matrices,
                abs_filter_value, rel_filter_value)  )
        
        write_p.start()    
        write_p.join()

        for p in ps:
            p.join()

    for ofp in ofps.values():
        ofp.close()
    
    return

if __name__ == "__main__":
    #import cProfile
    #cProfile.run( 'main()' )
    main()
