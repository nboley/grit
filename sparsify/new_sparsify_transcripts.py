import numpy
import sys, os
import math
import traceback
import time
import copy

from scipy.linalg import svd
from scipy.stats import chi2
from scipy.optimize import fminbound, brentq, bisect
from scipy.io import savemat

from StringIO import StringIO
from itertools import izip

import subprocess
import multiprocessing
from multiprocessing import Process

from math import sqrt

sys.path.append( os.path.join( os.path.dirname( __file__ ),
                               "..", "file_types", "fast_gtf_parser" ) )
from gtf import load_gtf, Transcript

sys.path.append( os.path.join(os.path.dirname( __file__ ), "..", "file_types" ))
                               
from reads import Reads, bin_reads

from f_matrix import calc_expected_cnts, find_nonoverlapping_boundaries, \
    build_nonoverlapping_indices

from frag_len import load_fl_dists, FlDist, build_normal_density

import cvxpy
import cvxopt
from cvxopt import solvers, matrix, spdiag, log, div, sqrt

MIN_TRANSCRIPT_FREQ = 1e-16
# finite differences step size
FD_SS = 1e-10
COMPARE_TO_DCP = False
num_threads = 1
CONV_EPS = 1e-2
NUM_ITER_FOR_CONV = 50



DEBUG_OPTIMIZATION = True

"""
            # otherwise the max step size is *not* the maximum. So we need to 
            # search for the max
            optimal_step_size = calc_optimal_step_size( x, gradient )
            optimal_step_lhd = f_lhd( x + optimal_step_size*gradient )
            if optimal_step_size <= max_feasible_step_size - 1e-10 \
                    and optimal_step_lhd > step_lhd:
                print "Optimal step max greater than max", \
                    optimal_step_lhd - step_lhd, \
                    optimal_step_lhd - lhd, optimal_step_size
                max_feasible_step_size = optimal_step_size
                break
            
            if lhd > step_lhd:
                print "LHD greater than max step size and optimal_step_size"
                break
                        
            #print "MADE IT!!!", raw_input()
            max_feasible_step_size, max_index = \
                calc_max_feasible_step_size_and_limiting_index( x, gradient )



        step_size = max_feasible_step_size
        step_lhd = f_lhd( x + step_size*gradient )
        while lhd > step_lhd:
            print "Step changed to:", step_size/2.
            step_size /= 2
            step_lhd = f_lhd( x + step_size*gradient )
            if step_size < 1e-12:
                print "Found step size:", step_size, "vs", max_feasible_step_size
                break




        # make sure that the small step size increases the likelihood
        start_lhd = f_lhd( project_onto_simplex(x) )
        while True:
            if f_lhd( project_onto_simplex(x+step_size*gradient) ) < start_lhd:
                if step_size < 1e-16:
                    step_size = 1e-6
                    N = 50
                    for loop in xrange( N ):
                        print loop, (loop/float(N))*step_size,  f_lhd(project_onto_simplex(x+(loop/float(N))*step_size*gradient)) - start_lhd
                    raw_input()

                    return 0
                print "DECREASING STEP SIZE", step_size/2, f_lhd( project_onto_simplex(x+step_size*gradient) ) - start_lhd
                step_size /= 2
            else:
                break
        
        big_step_lhd = f_lhd( project_onto_simplex(x+2*step_size*gradient) )
        if  big_step_lhd > start_lhd and \
                big_step_lhd <f_lhd(project_onto_simplex(x+step_size*gradient)):
            print "INCREASING STEP SIZE", step_size*2
            step_size *= 2
        
        x_new = project_onto_simplex(x + step_size*gradient)
        lhd = f_lhd(x_new)
        
        print "Step Size: %e - %e %e" % ( step_size, 
                                          f_lhd(x) - lhd,
                                          f_lhd(x_new) - lhd)

        if lhd < f_lhd( x ):
            print "WARNING: F MIN BOUND DIDNT WORK", lhd, f_lhd(x)
            N = 50
            for loop in xrange( N ):
                print loop,  (loop/float(N))*step_size,  f_lhd( project_onto_simplex(x + (loop/float(N))*step_size*gradient) )
            raw_input()


"""

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

def build_expected_and_observed_arrays( expected_cnts, observed_cnts ):
    expected_mat = []
    observed_mat = []
    
    for key, val in expected_cnts.iteritems():
        # skip bins with 0 expected reads
        if sum( val) == 0:
            continue
        
        expected_mat.append( val )
        try:
            observed_mat.append( observed_cnts[key] )
        except KeyError:
            observed_mat.append( 0 )

    expected_mat = numpy.array( expected_mat )
    nonzero_entries = expected_mat.sum(0).nonzero()[0]
    unobservable_transcripts = set(range(expected_mat.shape[1])) \
        - set(nonzero_entries.tolist())
    observed_mat = numpy.array( observed_mat )
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

def build_design_matrices( gene, bam_fname, fl_dists ):
    # load the bam file
    reads = Reads(bam_fname)
    
    # find the set of non-overlapping exons, and convert the transcripts to 
    # lists of these non-overlapping indices. All of the f_matrix code uses
    # this representation.     
    exon_boundaries = find_nonoverlapping_boundaries(gene.transcripts)
    transcripts_non_overlapping_exon_indices = \
        list(build_nonoverlapping_indices( gene.transcripts, exon_boundaries ))
    
    binned_reads = bin_reads( 
        reads, gene.chrm, gene.strand, exon_boundaries, False, True)
    if len( binned_reads ) == 0:
        raise ValueError, "TOO FEW READS"
    observed_cnts = build_observed_cnts( binned_reads, fl_dists )
    
    read_groups_and_read_lens =  { (RG, read_len) for RG, read_len, bin 
                                   in binned_reads.iterkeys() }
    
    fl_dists_and_read_lens = [ (fl_dists[RG], read_len) for read_len, RG  
                               in read_groups_and_read_lens ]
    
    expected_cnts = calc_expected_cnts( 
        exon_boundaries, transcripts_non_overlapping_exon_indices, 
        fl_dists_and_read_lens)
        
    expected_array, observed_array, unobservable_transcripts = \
        build_expected_and_observed_arrays( expected_cnts, observed_cnts )
    
    return expected_array, observed_array, unobservable_transcripts


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

def estimate_transcript_frequencies(  
        observed_array, expected_array,
        fixed_indices=[], fixed_values=[] ):
    def adjust_estimate(x):
        """Ensure the estimate is within the bounds

        """
        x_full = [ max( MIN_TRANSCRIPT_FREQ, x ) for x in x ]
        x_full_sum = sum(x_full)
        x_full = [ x/x_full_sum for x in x_full ]
        return x_full
    
    def f_lhd(x):
        log_lhd = calc_lhd(x, observed_array, expected_array)
        return log_lhd
    
    def f_gradient(x):
        deriv = calc_lhd_deriv( x, observed_array, expected_array )
        return deriv/deriv.sum()
    
    def project_onto_simplex_fast( x, debug=False ):
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


    def project_onto_simplex_generator(n, debug=False ):
        I = cvxopt.matrix(0.0, (n,n))
        I[::n+1] = 1.0
    
        G = cvxopt.matrix(0.0, (n,n))
        G[::n+1] = -1.0
        h = cvxopt.matrix(-MIN_TRANSCRIPT_FREQ, (n,1))
    
        A=matrix(1., (1,n))
        b=matrix(1., (1,1))
        cvxopt.solvers.options['show_progress'] = False
        
        def f(X):
            x = cvxopt.solvers.qp(P=I, q=-cvxopt.matrix(X), 
                                  G=G, h=h, A=A, b=b)['x']
            return numpy.array(x)[:,0]
        
        return f

    project_onto_simplex = project_onto_simplex_fast # project_onto_simplex_generator(expected_array.shape[1])
    
    def calc_max_feasible_step_size_and_limiting_index( x0, gradient ):
        """Calculate the maximum step size to stay in the feasible region.
        
        solve y - x*gradient = MIN_TRANSCRIPT_FREQ for x
        x = (y - MIN_TRANSCRIPT_FREQ)/gradient
        """
        # we use minus because we return a positive step
        steps = (x0-MIN_TRANSCRIPT_FREQ)/gradient
        step_size = -steps[ steps < 0 ].max()
        step_size_i = ( steps == -step_size ).nonzero()
        #print "DEBUG", x0[ step_size_i ] + step_size*gradient[ step_size_i ], \
        #    step_size_i, x0[step_size_i], gradient[step_size_i]
        return step_size, step_size_i
    
    def calc_search_direction( x ):
        gradient = f_gradient( x )
        x_next = project_onto_simplex( x + gradient )
        gradient = (x_next - x)
        return gradient/numpy.absolute(gradient).sum()

    def calc_optimal_step_size( x, gradient, ss=1e-4 ):
        # find the gradient
        d_l_1 = ( f_lhd( x ) - f_lhd( x + ss*gradient ) )/ss
        d_l_2 = ( f_lhd( x ) - f_lhd( x + 2*ss*gradient ) )/(2*ss)
        dd_l = (d_l_1 - d_l_2)/ss
        return d_l_1/dd_l

 
    def find_next_step( x, gradient ):
        # zero out indices until the max step isnt optimal anymore.
        while True:
            # find the maximum step size
            max_feasible_step_size, max_index = \
                calc_max_feasible_step_size_and_limiting_index( x, gradient )
            # check to see if the maximum is sub-optimal
            if max_feasible_step_size > FD_SS and \
                    f_lhd( x + (max_feasible_step_size-FD_SS)*gradient ) \
                    > f_lhd( x + max_feasible_step_size*gradient ):
                break
            
            # if it is optimal, zero out that gradient entry and continue
            x += max_feasible_step_size*gradient
            print "Zeroing ", max_index[0][0], gradient.sum()
            gradient[max_index] = 0
            # if the gradient norm has decreased too much, re-project it
            # onto the simplex
            if abs( gradient.sum() ) > 1e-3:
                #print "Gradient changed direction."
                x_next = project_onto_simplex( x + gradient )
                gradient = (x_next - x)
                gradient = gradient/numpy.absolute(gradient).sum()
            
            if f_lhd( x + FD_SS*gradient ) - f_lhd( x ) <= 0:
                gradient = calc_search_direction( x_prev )
            
            if ( gradient >= 0 ).all():
                gradient = calc_search_direction( x_prev )

        # now, do a line search
        def brentq_fmin(alpha):
            return f_lhd(x + (alpha+FD_SS)*gradient) \
                - f_lhd(x + (alpha-FD_SS)*gradient)
        
        try:
            step_size = brentq( brentq_fmin, FD_SS, 
                                max_feasible_step_size-FD_SS )
            #step_size = bisect(brentq_fmin, FD_SS, max_feasible_step_size-FD_SS)
        except ValueError:
            step_size = FD_SS
            lhd = f_lhd( x )
            while step_size > 1e-16 and lhd > f_lhd( x+step_size*gradient ):
                print "Setting step size to", step_size
                step_size /= 2
        
        print "Found step size:", step_size, "vs", max_feasible_step_size

        # code to over-project
        #if f_lhd( project_onto_simplex(x + 1e-2*gradient) ) \
        #        > f_lhd( x + step_size*gradient ):
        #    return project_onto_simplex(x + 1e-2*gradient)
        
        return x + step_size*gradient
    
    num_transcripts = expected_array.shape[1] - len( fixed_indices )
    bounds = [(MIN_TRANSCRIPT_FREQ,1)]*num_transcripts
    eq_value = 1. - sum( fixed_values )
    
    #x = nnls( expected_array, 1.*observed_array/observed_array.sum()  )
    #if DEBUG_OPTIMIZATION:
    #    print "FINISHED NNLS Estimate"
    x = numpy.ones(num_transcripts)/float(num_transcripts)
    lhds = [ -1e10 ]*NUM_ITER_FOR_CONV
    for i in xrange( 5000 ):
        x_prev = x
        # find the projected gradient
        gradient = calc_search_direction( x_prev )
        # calculate the next step
        x = find_next_step( x_prev, gradient )
        # sometimes the current estimate will wander from the simplex, so
        # we need to project it back into the simplex
        if abs( 1-x.sum() ) > 1e-6:
            x = project_onto_simplex(x)
            print "RE-PROJECTING"
            continue
        
        # stop conditions
        lhd = f_lhd( x )
        if lhd - min(lhds[-NUM_ITER_FOR_CONV:]) < CONV_EPS:
            break
        lhds.append( lhd )

        if DEBUG_OPTIMIZATION:
            print "%i\t%.2f\t%.2e\t%.2e" % ( 
                i, lhd, lhd - lhds[-2], lhd - min(lhds[-NUM_ITER_FOR_CONV:]) )
        
        #raw_input()
        
    log_lhd = calc_lhd( x, observed_array, expected_array )
    
    return x

def estimate_transcript_frequencies_new(  
        observed_array, expected_array,
        fixed_indices=[], fixed_values=[] ):
    """
    Until Convergence:

    Find the direction of maximum ascent
    Find the maximum step
    Check to see if the maximum step is the maximum
       - calculate f(max) - f(max-1e-6)
       if it is:
          zero out that coordinate, and re-project the gradient if necessary.
          continue
    # we know that the maximum isnt the max so, find the actual maximum via 
    # bisection. Continue...
        
    """
    pass

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

def estimate_gene_expression( expected_array, observed_array, unobservable_transcripts ):
    mle_estimate = estimate_transcript_frequencies( 
        observed_array, expected_array )
        
    if False and estimate_confidence_bounds:
        bnds = []
        for index, mle_value in enumerate( mle_estimate ):
            if index in unobservable_transcripts: continue
            
            transcript = gene.transcripts[index]
            bnds.append( [] )
            for bnd_type in ( "LOWER", "UPPER" ):
                p_value, bnd = estimate_confidence_bound( 
                    observed_array, expected_array, index, mle_estimate, bnd_type )
                bnds[-1].append( bnd )
            if VERBOSE: 
                print "Gene %s\tTranscript %s (%i/%i)\tBam %s\tEst %.2e\tBnds [%.2e %.2e]" % (
                    gene.id, transcript.id, index+1, len(gene.transcripts), 
                    os.path.basename(bam_fname), 
                    mle_value, bnds[-1][0], bnds[-1][1])

        lower_bnds, upper_bnds = zip( *bnds )
    else:
        lower_bnds, upper_bnds = [None]*len(mle_estimate), [None]*len(mle_estimate)
    
    return mle_estimate, lower_bnds, upper_bnds

def build_mle_estimate( observed_array, expected_array ):
    log_lhd, mle_estimate = estimate_transcript_frequencies( 
        observed_array, expected_array )

def build_bound(observed_array, expected_array, mle_estimate, index, bnd_type):
    p_value, bnd = estimate_confidence_bound( 
        observed_array, expected_array, index, mle_estimate, bnd_type )

def write_gene_to_gtf( ofp, gene, mles, lbs=None, ubs=None, filter_value=0 ):
    scores = [ int(1000*mle/max( mles )) for mle in mles ]
    
    for index, (mle, transcript) in enumerate(izip(mles, gene.transcripts)):
        if mle/max( mles ) < filter_value: continue
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
                                     estimate_confidence_bounds=False, 
                                     filter_value=1e-3):
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
            op_lock.release()
                        
            expected_array, observed_array, unobservable_transcripts \
                = build_design_matrices( gene, bam_fn, fl_dists )
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

            mle_estimate = \
                estimate_transcript_frequencies( observed_array, expected_array)
            log_lhd = calc_lhd( mle_estimate, observed_array, expected_array)
            taylor_log_lhd = calc_taylor_lhd( mle_estimate, observed_array, expected_array)
            if VERBOSE: print "FINISHED MLE %s\t%s\t%.2e\t%.2e\n%s" % ( 
                gene_id, bam_fn, log_lhd, taylor_log_lhd, 
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
                                 filter_value=0.0 ):
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
                print "Writing mat to '%s'" % ofname
                savemat( ofname, {'observed': observed, 'expected': expected} )
                print "Finished writing mat to '%s'" % ofname
                          
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
                                   lbs, ubs, filter_value)
                
                del output_dict[key]
                output_dict_lock.release()
        
        time.sleep(2.0)
    
    return

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='Determine valid transcripts and estimate frequencies.')
    parser.add_argument( 'ofprefix', \
        help='Output file name prefix. Output files will be ofprefix.bam_fn.gtf')
    parser.add_argument( 'gtf', type=file, \
        help='GTF file processed for expression')
    parser.add_argument( 'bam_fns', nargs='+', metavar='bam',\
        help='list of bam files to for which to produce expression')
    
    parser.add_argument( '--fl-dists', nargs='+', \
       help='a pickled fl_dist object(default:generate fl_dist from input bam)')
    parser.add_argument( '--fl-dist-norm', \
        help='mean and standard deviation (format "mn:sd") from which to ' \
            +'produce a fl_dist_norm (default:generate fl_dist from input bam)')

    parser.add_argument( '--threads', '-t', type=int , default=1, \
        help='Number of threads spawn for multithreading (default=1)')

    parser.add_argument( '--estimate-confidence-bounds', '-c', default=False,
        action="store_true",
        help='Whether or not to calculate confidence bounds ( this is slow )')

    parser.add_argument( '--write-design-matrices', default=False,
        action="store_true",
        help='Write the design matrices out to a matlab-style matrix file.')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
                             help='Whether or not to print status information.')
    parser.add_argument( '--debug-verbose', default=False, action='store_true',\
                             help='Prints the optimization path updates.')

    args = parser.parse_args()

    if not args.fl_dists and not args.fl_dist_norm:
        raise ValueError, "Must specific either --fl-dist or --fl-dist-norm."

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
        fl_dists, read_group_mappings = load_fl_dists( args.fl_dists )
    
    global DEBUG_VERBOSE
    DEBUG_VERBOSE = args.debug_verbose
    
    global VERBOSE
    VERBOSE = ( args.verbose or DEBUG_VERBOSE )

    # we change to the output directory later, and these files need to opened in
    # each sub-process for thread safety, so we get the absokute path while we 
    # can.
    bam_fns = [ os.path.abspath( bam_fn ) for bam_fn in args.bam_fns ]

    global PROCESS_SEQUENTIALLY
    if args.threads == 1:
        PROCESS_SEQUENTIALLY = True
    
    global num_threads
    num_threads = args.threads
    
    ofps = {}
    for bam_fn in bam_fns:
        ofps[bam_fn] = ThreadSafeFile(
            "%s.%s.gtf" % (args.ofprefix, os.path.basename(bam_fn)), "w")
        ofps[bam_fn].write(
            "track name=transcripts.%s useScore=1\n" % os.path.basename(bam_fn))
        
    return args.gtf, bam_fns, ofps, fl_dists, read_group_mappings, \
        args.estimate_confidence_bounds, args.write_design_matrices
    
def main():
    # Get file objects from command line
    gtf_fp, bam_fns, ofps, fl_dists, rg_mappings, \
        estimate_confidence_bounds, write_design_matrices = parse_arguments()
    filter_value = 0.01
    
    manager = multiprocessing.Manager()

    input_queue_lock = manager.Lock()
    input_queue = manager.list()

    output_dict_lock = manager.Lock()    
    output_dict = manager.dict()
    
    # add all the genes, in order of longest first. 
    genes = load_gtf( gtf_fp.name )
    for gene in sorted( genes, key=lambda x: -len(x.transcripts) ):
        for bam_fn in bam_fns:
            input_queue.append(('design_matrices', (gene.id, bam_fn, None)))
            gene_data = manager.dict()
            gene_data[ 'gene' ] = gene
            gene_data[ 'fl_dists' ] = fl_dists
            gene_data[ 'lbs' ] = {}
            gene_data[ 'ubs' ] = {}
            gene_data[ 'mle' ] = None
            gene_data[ 'design_matrices' ] = None
            output_dict[ (gene.id, bam_fn) ] = gene_data
    
    if 1 == num_threads:
        estimate_gene_expression_worker( 
            input_queue_lock, input_queue, 
            output_dict_lock, output_dict, 
            estimate_confidence_bounds,
            filter_value )
        write_finished_data_to_disk( 
            output_dict, output_dict_lock, ofps,
            estimate_confidence_bounds, 
            write_design_matrices,
            filter_value )
    else:
        ps = []
        for i in xrange( num_threads ):
            p = Process(target=estimate_gene_expression_worker, 
                        args=( input_queue_lock, input_queue, 
                               output_dict_lock, output_dict, 
                               estimate_confidence_bounds,
                               filter_value ) )
            p.start()
            ps.append( p )    

        write_p = Process(target=write_finished_data_to_disk, args=(
                output_dict, output_dict_lock, ofps,
                estimate_confidence_bounds,
                filter_value)  )

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
