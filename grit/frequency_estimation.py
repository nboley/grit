import os, sys
import time
import math

from itertools import izip

import numpy
numpy.seterr(all='ignore')

from scipy.linalg import svd, inv
from scipy.stats import chi2
from scipy.optimize import brentq, fminbound, nnls
from scipy.io import savemat

import config

import time
def make_time_str(et):
    hours = et//3600
    mins = et//60 - 60*hours
    secs = et - 3660*hours - 60*mins
    return "%i:%i:%.4f" % ( hours, mins, secs )

MIN_TRANSCRIPT_FREQ = 1e-12
# finite differences step size
FD_SS = 1e-10
NUM_ITER_FOR_CONV = 5
DEBUG_OPTIMIZATION = False
PROMOTER_SIZE = 50
LHD_ABS_TOL = 1e-5
PARAM_ABS_TOL = 1e-8

MAX_NUM_ITERATIONS = 1000

class TooFewReadsError( ValueError ):
    pass

"""
try:
    import sparsify_support_fns
    from sparsify_support_fns import (
        calc_lhd, calc_gradient, calc_hessian, 
        calc_penalized_lhd, calc_penalized_gradient )
except ImportError:
    raise
    def calc_lhd( freqs, observed_array, expected_array ):
        return float(observed_array*numpy.log( 
                numpy.matrix( expected_array )*numpy.matrix(freqs).T ))

    def calc_lhd_deriv( freqs, observed_array, expected_array ):
        denom = numpy.matrix( expected_array )*numpy.matrix(freqs).T
        rv = (((expected_array.T)*observed_array))*(1.0/denom)
        return -numpy.array(rv)[:,0]
"""

import sparsify_support_fns

def calc_lhd( freqs, observed_array, expected_array, 
              sparse_penalty=0, sparse_index=None ):
    rv = sparsify_support_fns.calc_lhd(freqs, observed_array, expected_array)

    if sparse_penalty > 0:
        if sparse_index != None:
            penalty = math.log(sparse_penalty) - math.log(freqs[sparse_index])
            rv -= math.exp(penalty)
        else:
            #rv -= sparse_penalty*(freqs**2).sum()
            rv += sparse_penalty*(freqs*numpy.log(freqs)).sum()

    return rv

def calc_gradient( freqs, observed_array, expected_array,
                   sparse_penalty=0, sparse_index=None):
    rv = sparsify_support_fns.calc_gradient(
        freqs, observed_array, expected_array)
    if sparse_penalty > 0:
        if sparse_index != None:
            penalty = math.log(sparse_penalty) - 2*math.log(freqs[sparse_index])
            rv[sparse_index] -= math.exp(penalty)
        else:
            #rv += sparse_penalty*2*freqs
            rv -= sparse_penalty*(1 + numpy.log(freqs))
    
    return rv

def is_row_identifiable(X, i_to_check):
    import scipy.optimize
    
    indices = numpy.array([i for i in xrange(X.shape[1]) if i != i_to_check]) 
    A = X[:,indices]
    b = X[:,i_to_check]
    res = scipy.optimize.nnls(A, b)
    return res[1] > 1e-6

def find_identifiable_transcripts( expected_array ):
    identifiable_transcripts = []
    for i in xrange(expected_array.shape[1]):
        if is_row_identifiable( expected_array, i ):
            identifiable_transcripts.append(i)
    return identifiable_transcripts


def nnls_cvxopt( X, Y, fixed_indices_and_values={} ):    
    from cvxopt import matrix, solvers
    X = matrix(X)
    Y = matrix(Y)
    
    m, n = X.size
    num_constraint = len( fixed_indices_and_values )
    
    G = matrix(0.0, (n,n))
    G[::n+1] = -1.0
    h = matrix(-MIN_TRANSCRIPT_FREQ, (n,1))

    # Add the equality constraints
    A= matrix(0., (1+num_constraint,n))
    b= matrix(0., (1+num_constraint,1))

    # Add the sum to one constraint
    A[0,:] = 1.
    b[0,0] = 1.
    
    # Add the fixed value constraints
    for const_i, (i, val) in enumerate(fixed_indices_and_values.iteritems()):
        A[const_i+1,i] = 1.
        b[const_i+1,0] = val
    
    solvers.options['show_progress'] = DEBUG_OPTIMIZATION
    res = solvers.qp(P=X.T*X, q=-X.T*Y, G=G, h=h, A=A, b=b)
    x =  numpy.array(res['x']).T[0,]
    rss = ((numpy.array(X*res['x'] - Y)[0,])**2).sum()
    
    if DEBUG_OPTIMIZATION:
        for key, val in res.iteritems():
            if key in 'syxz': continue
            config.log_statement( "%s:\t%s" % ( key.ljust(22), val ) )
        
        config.log_statement( "RSS: ".ljust(22) + str(rss) )

    x[x < MIN_TRANSCRIPT_FREQ] = MIN_TRANSCRIPT_FREQ
    x = project_onto_simplex(x)
    return x

def line_search( x, f, gradient, max_feasible_step_size ):
    """Calculate the optimal step to maximize f in the direction of gradient.

    """
    alpha = fminbound(lambda a: -f(project_onto_simplex(x + a*gradient)), 
                      0, max_feasible_step_size, 
                      xtol=1e-6)
    lhd0 = f(x)
    while alpha > 0 and f(project_onto_simplex(x+alpha*gradient)) < lhd0:
        alpha = alpha/2 - FD_SS
    return max(0, alpha)

def project_onto_simplex( x, debug=False ):
    if ( x >= 0 ).all() and abs( 1-x.sum()  ) < 1e-6: return x
    sorted_x = numpy.sort(x)[::-1]
    if debug: config.log_statement( "sorted x: %s" % sorted_x )
    n = len(sorted_x)
    if debug: config.log_statement( "cumsum: %s" % sorted_x.cumsum() )
    if debug: config.log_statement( "arange: %s" % numpy.arange(1,n+1) )
    rhos = sorted_x - (1./numpy.arange(1,n+1))*( sorted_x.cumsum() - 1 )
    if debug: config.log_statement( "rhos: %s" % rhos )
    try: 
        rho = (rhos > 0).nonzero()[0].max() + 1
    except: 
        raise
        x[x<MIN_TRANSCRIPT_FREQ] = MIN_TRANSCRIPT_FREQ
        x = x/x.sum()
        return x
    if debug: config.log_statement( "rho: %s" % rho )
    theta = (1./rho)*( sorted_x[:rho].sum()-1)
    if debug: config.log_statement( "theta: %s" % theta )
    x_minus_theta = x - theta
    if debug: config.log_statement( "x - theta: %s" % x_minus_theta )
    x_minus_theta[ x_minus_theta < 0 ] = MIN_TRANSCRIPT_FREQ
    return x_minus_theta

def calc_projected_gradient( x, expected_array, observed_array, step_size,
                             sparse_penalty=None, sparse_index=None  ):
    gradient = calc_gradient( 
        x, observed_array, expected_array, sparse_penalty, sparse_index )
    gradient = gradient/gradient.sum()
    x_next = project_onto_simplex( x + 10.*step_size*gradient )
    gradient = (x_next - x)
    return gradient


def calc_max_feasible_step_size_and_limiting_index_BAD( x0, gradient ):
    #Calculate the maximum step size to stay in the feasible region.
    #
    #solve y - x*gradient = MIN_TRANSCRIPT_FREQ for x
    #x = (y - MIN_TRANSCRIPT_FREQ)/gradient
    #
    # we use minus because we return a positive step
    try:
        steps = (x0-MIN_TRANSCRIPT_FREQ)/(gradient+1e-12)
        step_size = -steps[ steps < 0 ].max()
        step_size_i = ( steps == -step_size ).nonzero()[0]
    except:
        print steps
        raise
    return step_size, step_size_i

def calc_max_feasible_step_size_and_limiting_index( x0, gradient ):
    """Calculate the maximum step size to stay in the feasible region.

    solve y - x*gradient = MIN_TRANSCRIPT_FREQ for x
    x = (y - MIN_TRANSCRIPT_FREQ)/gradient
    """
    # we use minus because we return a positive step
    max_ss, max_i = -1e50, None
    for i, (x, dx) in enumerate(izip(x0, gradient)):
        if dx == 0: continue
        ss = (x - MIN_TRANSCRIPT_FREQ)/dx
        if ss >= 0: continue
        if ss > max_ss: 
            max_ss = ss
            max_i = i
    
    if max_i == None:
        return 0, 0
    return -max_ss, max_i

def build_zero_eliminated_matrices(x, full_expected_array, 
                                   curr_zeros, sparse_index):
    """Return x and an expected array without boundry points.
    
    """
    #current_nonzero_entries = (x > 1e-12).nonzero()[0]
    #if len( current_nonzero_entries ) > 1 \
    #        and len( current_nonzero_entries ) < len(x):
    n = full_expected_array.shape[1]
    full_x = numpy.ones(n)*MIN_TRANSCRIPT_FREQ
    full_x[ numpy.array(sorted(set(range(n))-curr_zeros)) ] = x

    zeros = set( (full_x <= 1e-12).nonzero()[0] ) - set((sparse_index,))
    nonzeros = sorted(set(range(n))-zeros)
    
    new_x = full_x[ nonzeros ]
    new_expected_array = full_expected_array[:, nonzeros]

    if sparse_index != None:
        sparse_index = sparse_index - sum(1 for i in zeros if i < sparse_index)
    return new_x, new_expected_array, zeros, sparse_index
   
def estimate_transcript_frequencies_line_search_OLD(  
        observed_array, full_expected_array, x0, 
        sparse_penalty, full_sparse_index,
        dont_zero, abs_tol ):    
    x = x0.copy()
    expected_array = full_expected_array.copy()
    n = full_expected_array.shape[1]
    
    prev_lhd = calc_lhd(x, observed_array, full_expected_array, 
                        sparse_penalty, full_sparse_index)
    lhds = [prev_lhd,]
    
    zeros = set()
    zeros_counter = 0

    sparse_index = full_sparse_index
    for i in xrange( MAX_NUM_ITERATIONS ):
        # calculate the gradient and the maximum feasible step size
        gradient = calc_projected_gradient( 
            x, expected_array, observed_array, 
            abs_tol,
            sparse_penalty, sparse_index )
        gradient_size = numpy.absolute(gradient).sum()
        # if the gradient is zero, then we have reached a maximum
        if gradient_size == 0: break
        gradient /= gradient_size
        max_feasible_step_size, max_index = \
            calc_max_feasible_step_size_and_limiting_index(x, gradient)
        
        # perform the line search
        alpha = line_search(
            x, lambda x: calc_lhd(x, observed_array, expected_array, 
                                  sparse_penalty, sparse_index), 
            gradient, max_feasible_step_size )
        x += alpha*gradient
        if abs( 1-x.sum() ) > 1e-6:
            x = project_onto_simplex(x)
            continue
     
        curr_lhd = calc_lhd( x, observed_array, expected_array, 
                             sparse_penalty, sparse_index)
        if i > 3 and (alpha == 0 or curr_lhd - prev_lhd < abs_tol):
            zeros_counter += 1
            if zeros_counter > 3:
                break            
        else:
            zeros_counter = 0
            if not dont_zero:
                x, expected_array, zeros, sparse_index = build_zero_eliminated_matrices(
                    x, full_expected_array, zeros, full_sparse_index)
                if len(x) <= 2:
                    break
        
        prev_lhd = curr_lhd
        lhds.append( curr_lhd )
        if DEBUG_OPTIMIZATION:
            config.log_statement( "%i\t%.2f\t%.6e\t%i" % ( 
                    i, curr_lhd, curr_lhd - prev_lhd, len(x) ) )
    
    final_x = numpy.ones(n)*MIN_TRANSCRIPT_FREQ
    final_x[ numpy.array(sorted(set(range(n))-zeros)) ] = x
    final_lhd = calc_lhd(final_x, observed_array, full_expected_array, 
                         sparse_penalty, full_sparse_index)

    if final_lhd < prev_lhd:
        return x0, [prev_lhd,]
    #assert numpy.isnan(prev_lhd) or final_lhd - prev_lhd > -abs_tol, \
    #    "Final: %e\tPrev: %e\tDiff: %e\tTol: %e\t %s" % (
    #    final_lhd, prev_lhd, final_lhd-prev_lhd, abs_tol, lhds)
        
    return final_x, lhds

def estimate_transcript_frequencies_line_search(  
        observed_array, full_expected_array, x0, 
        sparse_penalty, full_sparse_index,
        dont_zero, abs_tol ):    
    x = x0.copy()
    expected_array = full_expected_array.copy()
    n = full_expected_array.shape[1]
    
    prev_lhd = calc_lhd(x, observed_array, full_expected_array, 
                        sparse_penalty, full_sparse_index)
    lhds = [prev_lhd,]
    
    zeros = set()
    zeros_counter = 0
    
    sparse_index = full_sparse_index
    for i in xrange( MAX_NUM_ITERATIONS ):
        gradient = calc_gradient( 
            x, observed_array, expected_array, sparse_penalty, sparse_index )
        gradient = gradient/(gradient.sum() + 1e-12)

        def line_search_f(alpha):
            new_x = project_onto_simplex(x + alpha*gradient)
            return calc_lhd(new_x, observed_array, expected_array, 
                             sparse_penalty, sparse_index)
        
        # bisection on alpha
        min_alpha, min_val = 0, prev_lhd
        max_alpha, max_val = 10, line_search_f(10.)
        while True:
            if max_alpha - min_alpha < abs_tol and min_alpha > 1e-12:
                break
            alpha = (max_alpha + min_alpha)/2.
            if alpha < 1e-12: break
            
            val = line_search_f(alpha)
            if val > min_val: 
                min_alpha = alpha
                min_val = val
            else: 
                max_alpha = alpha
                max_val = val
        
        new_x = project_onto_simplex(x + alpha*gradient)
        curr_lhd = calc_lhd( new_x, observed_array, expected_array, 
                             sparse_penalty, sparse_index)
        if curr_lhd > prev_lhd: x = new_x
        else: curr_lhd = prev_lhd
        
        assert curr_lhd >= prev_lhd, "%e %e %e" % (curr_lhd, prev_lhd, curr_lhd - prev_lhd)
        if i > 3 and (alpha == 0 or curr_lhd - prev_lhd < abs_tol):
            zeros_counter += 1
            if zeros_counter > 3:
                break
        else:
            zeros_counter = 0
            if not dont_zero:
                x, expected_array, zeros, sparse_index = build_zero_eliminated_matrices(
                    x, full_expected_array, zeros, full_sparse_index)
                if len(x) <= 2:
                    break
        
        prev_lhd = curr_lhd
        lhds.append( curr_lhd )
        if DEBUG_OPTIMIZATION:
            config.log_statement( "%i\t%.2f\t%.6e\t%i" % ( 
                    i, curr_lhd, curr_lhd - prev_lhd, len(x) ) )
    
    final_x = numpy.ones(n)*MIN_TRANSCRIPT_FREQ
    final_x[ numpy.array(sorted(set(range(n))-zeros)) ] = x
    final_x = project_onto_simplex(final_x)
    final_lhd = calc_lhd(final_x, observed_array, full_expected_array, 
                         sparse_penalty, full_sparse_index)

    if final_lhd < prev_lhd:
        return x0, [prev_lhd,]
    #assert numpy.isnan(prev_lhd) or final_lhd - prev_lhd > -abs_tol, \
    #    "Final: %e\tPrev: %e\tDiff: %e\tTol: %e\t %s" % (
    #    final_lhd, prev_lhd, final_lhd-prev_lhd, abs_tol, lhds)
        
    return final_x, lhds

def estimate_transcript_frequencies_with_cvxopt( 
        observed_array, expected_array, sparse_penalty, sparse_index,
        verbose=False):
    from cvxpy import matrix, variable, geq, log, eq, program, maximize, \
        minimize, sum, quad_over_lin
    
    Xs = matrix( observed_array )
    ps = matrix( expected_array )
    thetas = variable( ps.shape[1] )
    constraints = [ eq(sum(thetas), 1), geq(thetas,0)]

    if sparse_penalty == None:
        p = program( maximize(Xs*log(ps*thetas)), constraints )
    else:
        p = program( maximize(Xs*log(ps*thetas) - sparse_penalty*quad_over_lin(
                    1., thetas[sparse_index,0])), constraints )
        
    
    p.options['maxiters']  = 1500
    value = p.solve(quiet=not verbose)
    thetas_values = numpy.array(thetas.value.T.tolist()[0])
    return thetas_values
    

def estimate_transcript_frequencies_sparse(  
        observed_array, full_expected_array,
        min_sparse_penalty, sparse_index ):
    if observed_array.sum() == 0:
        raise TooFewReadsError, (
            "Too few reads (%i)" % observed_array.sum() )
    
    n = full_expected_array.shape[1]
    if n == 1:
        return numpy.ones( 1, dtype=float )
    
    #x = project_onto_simplex(nnls(
    #        full_expected_array, observed_array)[0])
    x = numpy.ones(n, dtype=float)/n
    eps = 0.1
    sparse_penalty = 100
    if min_sparse_penalty == None:
        min_sparse_penalty = LHD_ABS_TOL
        # always skip the out of gene transcript
        #sparse_index = numpy.argmax(x[1:]) + 1
    
    start_time = time.time()
    #config.log_statement( "Iteration\tlog lhd\t\tchange lhd\tn 
    # iter\ttolerance\ttime (hr:min:sec)" )
    for i in xrange( MAX_NUM_ITERATIONS ):
        prev_x = x.copy()
        
        x, lhds = estimate_transcript_frequencies_line_search(  
            observed_array, full_expected_array, x, 
            sparse_penalty, sparse_index,
            dont_zero=False, abs_tol=eps )
        lhd = calc_lhd( x, observed_array, full_expected_array, 
                        sparse_penalty, sparse_index )
        prev_lhd = calc_lhd( prev_x, observed_array, full_expected_array,
                             sparse_penalty, sparse_index)
        if config.DEBUG_VERBOSE:
            config.log_statement( "Zeroing %i\t%.2f\t%.2e\t%i\t%e\t%e\t%s" % ( 
                i, lhd, (lhd - prev_lhd)/len(lhds), len(lhds ), eps, 
                sparse_penalty,
                make_time_str((time.time()-start_time)/len(lhds)) ) )
            
        start_time = time.time()
        
        if float(lhd - prev_lhd)/len(lhds) < eps:
            if eps == LHD_ABS_TOL and sparse_penalty == min_sparse_penalty:
                break
            eps = max( eps/10, LHD_ABS_TOL )
            if sparse_penalty != None:
                sparse_penalty = max( sparse_penalty/5, min_sparse_penalty )
    
    for i in xrange( 10 ):
        prev_x = x.copy()
        x, lhds = estimate_transcript_frequencies_line_search(  
            observed_array, full_expected_array, x, 
            sparse_penalty, sparse_index,
            dont_zero=True, abs_tol=LHD_ABS_TOL )
        lhd = calc_lhd( x, observed_array, full_expected_array, 
                        sparse_penalty, sparse_index)
        prev_lhd = calc_lhd( prev_x, observed_array, full_expected_array,
                             sparse_penalty, sparse_index)
        if config.DEBUG_VERBOSE:
            config.log_statement( "Non-Zeroing %i\t%.2f\t%.2e\t%i\t%e\t%s" % ( 
                i, lhd, (lhd - prev_lhd)/len(lhds), len(lhds), eps,
                make_time_str((time.time()-start_time)/len(lhds))) )
        
        start_time = time.time()
        if len( lhds ) < 500: break
    
    return x

def estimate_sparse_transcript_frequencies(observed_array, full_expected_array):
    sparse_penalty = 1e-12

    best_x = None
    best_lhd = -1e100

    for sparse_index in xrange(1, full_expected_array.shape[1]):
        #cvx_sln = estimate_transcript_frequencies_with_cvxopt( 
        #    observed_array, full_expected_array, 
        #    sparse_penalty, sparse_index, 
        #    False )
        
        x = estimate_transcript_frequencies_sparse(
            observed_array, full_expected_array, 
            sparse_penalty, sparse_index )
        
        lhd = calc_lhd( x, observed_array, full_expected_array, 
                        sparse_penalty, sparse_index)
        if lhd > best_lhd:
            best_x = x
    
    return best_x

def estimate_transcript_frequencies(observed_array, full_expected_array):
    rv = estimate_transcript_frequencies_sparse(
        observed_array, full_expected_array, 
        None, None )
    return rv
        

def estimate_confidence_bound( f_mat, 
                               num_reads_in_bams,
                               fixed_index,
                               mle_estimate,
                               bound_type,
                               alpha):
    if bound_type == 'lb': bound_type = 'LOWER'
    if bound_type == 'ub': bound_type = 'UPPER'
    assert bound_type in ('LOWER', 'UPPER'), (
        "Improper bound type '%s'" % bound_type )
    expected_array, observed_array = f_mat.expected_and_observed(
        bam_cnts=num_reads_in_bams)
    eps = 0.1
    #expected_array = expected_array[:, 0:4]    
    #mle_estimate = mle_estimate[:-1]
    #print mle_estimate
    if 1 == expected_array.shape[1]:
        return 1.0, 1.0
    
    def min_line_search( x, gradient, max_feasible_step_size ):
        def brentq_fmin(alpha):
            return calc_lhd(project_onto_simplex(x + alpha*gradient), 
                            observed_array, expected_array) - min_lhd
        
        min_step_size = 0
        max_step_size = max_feasible_step_size

        assert brentq_fmin(min_step_size) >= 0, "brent min is %e" % brentq_fmin(min_step_size)
        if brentq_fmin(max_step_size) >= 0:
            return max_step_size, False
        
        # do a line search with brent
        step_size = brentq(brentq_fmin, min_step_size, max_step_size )
        # make sure that we're above the minimum lhd
        while brentq_fmin(step_size) < 0 and step_size > 0:
            step_size = step_size/2 - FD_SS
        rv = max(0, step_size)
        
        assert calc_lhd(project_onto_simplex(x+rv*gradient),
                        observed_array,expected_array) >= min_lhd
        return rv, True
    
    def take_param_decreasing_step(x):
        gradient = numpy.zeros(n)
        gradient[fixed_index] = -1 if bound_type == 'LOWER' else 1
        gradient = project_onto_simplex( x + 100*eps*gradient ) - x
        gradient_l1_size = numpy.absolute(gradient).sum()
        # if we can't go anywhere, then dont move
        if gradient_l1_size > 0:
            gradient /= numpy.absolute(gradient).sum()        
            assert not numpy.isnan(gradient).any()
            max_feasible_step_size, max_index = \
                calc_max_feasible_step_size_and_limiting_index(x, gradient)
            alpha, is_full_step = min_line_search(
                x, gradient, max_feasible_step_size)
            x = project_onto_simplex(x + alpha*gradient)
        
        return x
    
    def take_lhd_increasing_step(x):
        # find the simple lhd gradient at this point
        lhd_gradient = calc_projected_gradient( 
            x, expected_array, observed_array, 10*eps )
        lhd_gradient /= ( numpy.absolute(lhd_gradient).sum() + 1e-12 )
        
        # find the projected gradient to minimize x[fixed_index]
        coord_gradient = numpy.zeros(n)
        coord_gradient[fixed_index] = -1 if bound_type == 'LOWER' else 1
        coord_gradient = project_onto_simplex( x + 100*eps*coord_gradient ) - x
        coord_gradient /= ( numpy.absolute(coord_gradient).sum() + 1e-12 )
        
        # if the lhd step is already moving the paramater of itnerest in the 
        # proper direction, we just take a normal lhd step
        if ( bound_type == 'LOWER' 
             and lhd_gradient[fixed_index] <= PARAM_ABS_TOL ) \
                or ( bound_type == 'UPPER' 
                     and lhd_gradient[fixed_index] >= -PARAM_ABS_TOL ):
            theta = 0
        # otherwise, find out how miuch we need to step off of the maximum step
        # to make our parameter not move in the wrong direction
        else:
            # we want lhd_gradient[fixed_index]*(1-theta) + \
            #   theta*coord_gradient[fixed_index] == 0
            # implies lhd[i] -theta*lhd[i] + theta*cg[i] == 0
            #     =>  lhd[i] = theta*lhd[i] - theta*cg[i]
            #     =>  lhd[i]/(lhd[i] - theta*cg[i]) = theta*
            theta = lhd_gradient[fixed_index]/(
                lhd_gradient[fixed_index] - coord_gradient[fixed_index])
        
        gradient = (1-theta)*lhd_gradient + theta*coord_gradient
        projection = project_onto_simplex( x + 100*eps*gradient )
        gradient = projection - x
        gradient /= ( gradient.sum() + 1e-6 )
        #print theta, x
        #print gradient
        
        max_feasible_step_size, max_index = \
            calc_max_feasible_step_size_and_limiting_index(x, gradient)
        alpha = line_search(
            x, lambda x: calc_lhd(x, observed_array, expected_array),
            gradient, max_feasible_step_size)
        new_x = project_onto_simplex(x+alpha*gradient)
        
        assert calc_lhd(new_x, observed_array, expected_array) >= min_lhd - 1e-12,\
            "Value is %e %s %e %e %e %e" % (
            alpha,
            new_x,
            calc_lhd(new_x, observed_array, expected_array), 
            calc_lhd(project_onto_simplex(x), observed_array, expected_array), 
            min_lhd - 1e-12, 
            calc_lhd(new_x, observed_array, expected_array) - min_lhd - 1e-12 )
        return new_x
    
    n = expected_array.shape[1]    
    max_lhd = calc_lhd(
        project_onto_simplex(mle_estimate), observed_array, expected_array)
    unprojected_lhd = calc_lhd(mle_estimate, observed_array, expected_array)
    assert abs(max_lhd - unprojected_lhd) < 1e-2, "Diff: %e %e %e" % (
        max_lhd, unprojected_lhd, max_lhd - unprojected_lhd)
    max_test_stat = chi2.ppf( 1 - alpha, 1 )/2.    
    min_lhd = max_lhd-max_test_stat

    x = mle_estimate.copy()
    prev_x = mle_estimate[fixed_index]
    n_successes = 0
    for i in xrange(MAX_NUM_ITERATIONS):
        # take a downhill step
        x = take_param_decreasing_step(x)
        curr_lhd = calc_lhd(x, observed_array, expected_array)
        if bound_type == 'LOWER' \
                and abs(x[fixed_index] - MIN_TRANSCRIPT_FREQ) < PARAM_ABS_TOL \
                and curr_lhd >= min_lhd:
            break
        if bound_type == 'UPPER' \
                and abs(x[fixed_index] - (1.0-MIN_TRANSCRIPT_FREQ)) < PARAM_ABS_TOL \
                and curr_lhd >= min_lhd:
            break

        x = take_lhd_increasing_step(x)
        
        if abs( prev_x - x[fixed_index] ) < eps:
            if eps > PARAM_ABS_TOL:
                eps = max(eps/2, PARAM_ABS_TOL)
            else:
                n_successes += 1
                if n_successes > 3:
                    break
        else: 
            prev_x = x[fixed_index]
            n_successes = 0

    value = x[fixed_index] 
    if value < PARAM_ABS_TOL: value = 0.
    if 1-value < PARAM_ABS_TOL: value = 1.
    lhd = calc_lhd( x, observed_array, expected_array )
    rv = chi2.sf( 2*(max_lhd-lhd), 1), value
    return rv    

def estimate_confidence_bound_with_cvx( f_mat, 
                               num_reads_in_bams,
                               fixed_i,
                               mle_estimate,
                               bound_type,
                               alpha):
    if bound_type == 'lb': bound_type = 'LOWER'
    if bound_type == 'ub': bound_type = 'UPPER'
    assert bound_type in ('LOWER', 'UPPER'), (
        "Improper bound type '%s'" % bound_type )
    expected_array, observed_array = f_mat.expected_and_observed(
        bam_cnts=num_reads_in_bams)
    
    from cvxpy import matrix, variable, geq, log, eq, program, maximize, minimize, sum
    
    mle_log_lhd = calc_lhd(mle_estimate, observed_array, expected_array)
    
    lower_lhd_bound = mle_log_lhd - chi2.ppf( 1 - alpha, 1 )/2.
    free_indices = set(range(expected_array.shape[1])) - set((fixed_i,))
    
    Xs = matrix( observed_array )
    ps = matrix( expected_array )
    thetas = variable( ps.shape[1] )
    constraints = [ geq(Xs*log(ps*thetas), lower_lhd_bound), 
                    eq(sum(thetas), 1), geq(thetas,0)]
    
    if bound_type == 'UPPER':
        p = program( maximize(thetas[fixed_i,0]), constraints )    
    else:
        p = program( minimize(thetas[fixed_i,0]), constraints )
    p.options['maxiters']  = 1500
    value = p.solve(quiet=not DEBUG_OPTIMIZATION)
    thetas_values = numpy.array(thetas.value.T.tolist()[0])
    log_lhd = calc_lhd( thetas_values, observed_array, expected_array )
    
    return chi2.sf( 2*(mle_log_lhd-log_lhd), 1), value
