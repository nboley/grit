import os, sys
import time

from itertools import izip

import numpy
numpy.seterr(all='ignore')

from scipy.linalg import svd, inv
from scipy.stats import chi2
from scipy.optimize import fminbound, brentq, bisect, line_search
from scipy.io import savemat

import time
def make_time_str(et):
    hours = et//3600
    mins = et//60 - 60*hours
    secs = et - 3660*hours - 60*mins
    return "%i:%i:%.4f" % ( hours, mins, secs )

VERBOSE = False
DEBUG_VERBOSE = False

def log_statement(message, only_log=False ):
    print >> sys.stderr, message

MIN_TRANSCRIPT_FREQ = 1e-12
# finite differences step size
FD_SS = 1e-12
NUM_ITER_FOR_CONV = 5
DEBUG_OPTIMIZATION = False
PROMOTER_SIZE = 50
ABS_TOL = 1e-5

class TooFewReadsError( ValueError ):
    pass

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
            log_statement( "%s:\t%s" % ( key.ljust(22), val ) )
        
        log_statement( "RSS: ".ljust(22) + str(rss) )
    
    return x

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

def estimate_confidence_bounds_directly( 
        observed_array, expected_array, fixed_i, 
        mle_log_lhd, upper_bound=True, alpha=0.10 ):
    assert upper_bound in ( True, False )
    from cvxpy import matrix, variable, geq, log, eq, program, maximize, minimize, sum
    lower_lhd_bound = mle_log_lhd - chi2.ppf( 1 - alpha, 1 )/2.
    free_indices = set(range(expected_array.shape[1])) - set((fixed_i,))
    
    Xs = matrix( observed_array )
    ps = matrix( expected_array )
    thetas = variable( ps.shape[1] )
    constraints = [ geq(Xs*log(ps*thetas), lower_lhd_bound), 
                    eq(sum(thetas), 1), geq(thetas,0)]
    if upper_bound:
        p = program( maximize(thetas[fixed_i,0]), constraints )    
    else:
        p = program( minimize(thetas[fixed_i,0]), constraints )
    
    p.options['maxiters']  = 1500
    value = p.solve(quiet=not DEBUG_OPTIMIZATION)
    
    thetas_values = numpy.array(thetas.value.T.tolist()[0])
    log_lhd = calc_lhd( thetas_values, observed_array, expected_array )
    
    return chi2.sf( 2*(mle_log_lhd-log_lhd), 1), value


def project_onto_simplex( x, debug=False ):
    if ( x >= MIN_TRANSCRIPT_FREQ ).all() and abs( 1-x.sum()  ) < 1e-6: return x
    sorted_x = numpy.sort(x)[::-1]
    if debug: log_statement( "sorted x: %s" % sorted_x )
    n = len(sorted_x)
    if debug: log_statement( "cumsum: %s" % sorted_x.cumsum() )
    if debug: log_statement( "arange: %s" % numpy.arange(1,n+1) )
    rhos = sorted_x - (1./numpy.arange(1,n+1))*( sorted_x.cumsum() - 1 )
    if debug: log_statement( "rhos: %s" % rhos )
    rho = (rhos > 0).nonzero()[0].max() + 1
    if debug: log_statement( "rho: %s" % rho )
    theta = (1./rho)*( sorted_x[:rho].sum()-1)
    if debug: log_statement( "theta: %s" % theta )
    x_minus_theta = x - theta
    if debug: log_statement( "x - theta: %s" % x_minus_theta )
    x_minus_theta[ x_minus_theta < 0 ] = MIN_TRANSCRIPT_FREQ
    return x_minus_theta

def estimate_transcript_frequencies_line_search(  
        observed_array, full_expected_array, x0, 
        dont_zero, abs_tol,
        fixed_indices=[], fixed_values=[] ):
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
        try:
            steps = (x0-MIN_TRANSCRIPT_FREQ)/(gradient+1e-12)
            step_size = -steps[ steps < 0 ].max()
            step_size_i = ( steps == -step_size ).nonzero()[0]
        except:
            print steps
            raise
        return step_size, step_size_i
    
    def calc_projected_gradient( x ):
        gradient = f_gradient( x )
        gradient = gradient/gradient.sum()
        x_next = project_onto_simplex( x + 1.*gradient )
        gradient = (x_next - x)
        for i, val in izip( fixed_indices, fixed_values ):
            gradient[i] = val
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
                # if there are new non-zero entries, then remove them
                # ( unless that would drop us down to 1 parameter, which
                #   doesn't make sense )
                if len( current_nonzero_entries ) > 1 \
                        and len( current_nonzero_entries ) < len(x):
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
            log_statement( "%i\t%.2f\t%.6e\t%i" % ( 
                    i, lhd, lhd - prev_lhd, len(x) ) )
    
    final_x = numpy.ones(n)*MIN_TRANSCRIPT_FREQ
    final_x[ numpy.array(sorted(set(range(n))-zeros)) ] = x
    final_lhd = calc_lhd(final_x, observed_array, full_expected_array)
    assert final_lhd >= f_lhd(x) - abs_tol
    return final_x, lhds

def estimate_transcript_frequencies(  
        observed_array, full_expected_array ):
    if observed_array.sum() == 0:
        raise TooFewReadsError, "Too few reads (%i)" % observed_array.sum()
    
    n = full_expected_array.shape[1]
    if n == 1:
        return numpy.ones( 1, dtype=float )
    
    x = numpy.array([1./n]*n)
    eps = 10.
    start_time = time.time()
    #log_statement( "Iteration\tlog lhd\t\tchange lhd\tn iter\ttolerance\ttime (hr:min:sec)" )
    for i in xrange( 500 ):
        prev_x = x.copy()
        
        x, lhds = estimate_transcript_frequencies_line_search(  
            observed_array, full_expected_array, x, 
            dont_zero=False, abs_tol=eps )
        
        lhd = calc_lhd( x, observed_array, full_expected_array )
        prev_lhd = calc_lhd( prev_x, observed_array, full_expected_array )
        if DEBUG_VERBOSE:
            log_statement( "Zeroing %i\t%.2f\t%.2e\t%i\t%e\t%s" % ( 
                i, lhd, (lhd - prev_lhd)/len(lhds), len(lhds ), eps, 
                make_time_str((time.time()-start_time)/len(lhds)) ) )
            
        start_time = time.time()
        
        if float(lhd - prev_lhd)/len(lhds) < eps:
            if eps == ABS_TOL: break
            eps /= 5
            eps = max( eps, ABS_TOL )
        
    
    for i in xrange( 5 ):
        prev_x = x.copy()
        x, lhds = estimate_transcript_frequencies_line_search(  
            observed_array, full_expected_array, x, 
            dont_zero=True, abs_tol=ABS_TOL )
        lhd = calc_lhd( x, observed_array, full_expected_array )
        prev_lhd = calc_lhd( prev_x, observed_array, full_expected_array )
        if DEBUG_VERBOSE:
            log_statement( "Non-Zeroing %i\t%.2f\t%.2e\t%i\t%e\t%s" % ( 
                i, lhd, (lhd - prev_lhd)/len(lhds), len(lhds), eps,
                make_time_str((time.time()-start_time)/len(lhds))) )
        
        start_time = time.time()
        if len( lhds ) < 500: break
    
    return x

def estimate_confidence_bound( observed_array, 
                               expected_array, 
                               fixed_index,
                               mle_estimate,
                               bound_type,
                               alpha):
    n = expected_array.shape[1]
    if n == 1:
        return 1.0, 1.0
    
    def f_lhd(x):
        log_lhd = calc_lhd(x, observed_array, expected_array)
        return log_lhd
    
    def f_gradient(x):
        return calc_gradient( x, observed_array, expected_array )

    def f_hessian(x):
        return calc_hessian( x, observed_array, expected_array )

    max_lhd = f_lhd(mle_estimate)
    max_test_stat = chi2.ppf( 1 - alpha, 1 )/2.    
    min_lhd = max_lhd-max_test_stat
        
    def calc_max_feasible_step_size_and_limiting_index( x0, gradient ):
        """Calculate the maximum step size to stay in the feasible region.
        
        solve y - x*gradient = MIN_TRANSCRIPT_FREQ for x
        x = (y - MIN_TRANSCRIPT_FREQ)/gradient
        """
        # we use minus because we return a positive step
        steps = (x0-MIN_TRANSCRIPT_FREQ)/(gradient+1e-12)
        try:
            step_size = -steps[ steps < 0 ].max()
        except:
            print steps
            raise
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
            step_size = FD_SS
            curr_lhd = f_lhd( x )
            while step_size > FD_SS and curr_lhd > f_lhd(x+step_size*gradient):
                step_size /= 1.5
            return int(step_size> FD_SS)*step_size
        
        min_step_size = FD_SS
        max_step_size = max_feasible_step_size-FD_SS
        if max_step_size < 0: 
            return 0, True        
        if brentq_fmin(max_step_size) >= 0:
            return max_step_size, False
        if brentq_fmin(min_step_size) <= 0:
            return 0, True

        # do a line search with brent
        step_size = brentq(brentq_fmin, min_step_size, max_step_size )
        if f_lhd(x) > f_lhd(x+step_size*gradient):
            step_size = downhill_search( step_size )
            return step_size, (step_size==0)
        
        return step_size, True

    def min_line_search( x, gradient, max_feasible_step_size ):
        def brentq_fmin(alpha):
            return f_lhd(x + alpha*gradient) - min_lhd
                
        min_step_size = 0
        max_step_size = max_feasible_step_size

        if brentq_fmin(min_step_size) < 0:
            print x
            print f_lhd(x), min_lhd, f_lhd(x)-min_lhd, brentq_fmin(min_step_size), max_lhd
            assert False
        
        if brentq_fmin(max_step_size) >= 0:
            return max_step_size, False
        
        # do a line search with brent
        step_size = brentq(brentq_fmin, min_step_size, max_step_size )
        while brentq_fmin(step_size) < 0 and step_size > 0:
            step_size -= 1e-12

        rv = max(0, step_size)
        assert f_lhd(x+rv*gradient) >= min_lhd
        return rv, True
    
    def calc_joint_gradient(x, theta):
        # find the simple lhd gradient at this point
        lhd_gradient = calc_projected_gradient( x )
        lhd_gradient /= numpy.absolute(lhd_gradient).sum()

        # find teh projected gradient to minimize x[fixed_index]
        coord_gradient = numpy.zeros(n)
        coord_gradient[fixed_index] = -1 if bound_type == 'LOWER' else 1
        coord_gradient = project_onto_simplex( x + 0.1*coord_gradient ) - x
        coord_gradient /= numpy.absolute(coord_gradient).sum()        

        gradient = (1-theta)*lhd_gradient + theta*coord_gradient
        gradient = project_onto_simplex( x + 0.1*gradient ) - x
        gradient /= numpy.absolute(gradient).sum()                
        
        return gradient
    
    def take_param_decreasing_step(x):
        gradient = calc_joint_gradient( x, 1.0 )
        max_feasible_step_size, max_index = \
            calc_max_feasible_step_size_and_limiting_index(x, gradient)
        alpha, is_full_step = min_line_search(
            x, gradient, max_feasible_step_size)
        x += alpha*gradient    
        lhd = f_lhd(x)
        assert lhd >= min_lhd
        if DEBUG_OPTIMIZATION:
            print "DOWNHILL", x[fixed_index], lhd - min_lhd, \
                "MAX STEP:", max_feasible_step_size, "REAL STEP", alpha
        return x, lhd
    
    def take_lhd_decreasing_step_wth_bisection(x):
        # find the maximum theta st gradient[fixed_index] >= 0
        def brentq_theta_min(theta):
            return calc_joint_gradient( x, theta )[fixed_index]

        try:
            theta = brentq(brentq_theta_min, 0., 1. )
        except:
            theta = 0

        gradient = calc_joint_gradient( x, theta )
        max_feasible_step_size, max_index = \
            calc_max_feasible_step_size_and_limiting_index(x, gradient)
        alpha, is_full_step = line_search(
            x, gradient, max_feasible_step_size)
        assert alpha >= 0
        x = x + alpha*gradient
        if DEBUG_OPTIMIZATION:
            print "MAX LHD ", x[fixed_index], f_lhd(x) - min_lhd, \
                "MAX STEP:", max_feasible_step_size, "REAL STEP", alpha
        assert f_lhd(x) >= min_lhd
        return x

    def take_lhd_decreasing_step(x):
        # find the maximum theta st gradient[fixed_index] >= 0
        def brentq_theta_min(theta):
            return calc_joint_gradient( x, theta )[fixed_index]

        try:
            theta = brentq(brentq_theta_min, 0., 1. )
        except:
            theta = 0

        # find the simple lhd gradient at this point
        lhd_gradient = calc_projected_gradient( x )
        lhd_gradient /= numpy.absolute(lhd_gradient).sum()
        
        # find teh projected gradient to minimize x[fixed_index]
        coord_gradient = numpy.zeros(n)
        coord_gradient[fixed_index] = -1 if bound_type == 'LOWER' else 1
        coord_gradient = project_onto_simplex( x + 0.1*coord_gradient ) - x
        coord_gradient /= numpy.absolute(coord_gradient).sum()        

        # find the theta that makes the paramater we are optimizing fixed
        theta = lhd_gradient[fixed_index]/(
            lhd_gradient[fixed_index] - coord_gradient[fixed_index])
        
        gradient = (1-theta)*lhd_gradient + theta*coord_gradient
        gradient = project_onto_simplex( x + 0.1*gradient ) - x
        gradient /= numpy.absolute(gradient).sum()                
        
        max_feasible_step_size, max_index = \
            calc_max_feasible_step_size_and_limiting_index(x, gradient)
        alpha, is_full_step = line_search(
            x, gradient, max_feasible_step_size)
        assert alpha >= 0
        x = x + alpha*gradient
        if DEBUG_OPTIMIZATION:
            print "MAX LHD ", x[fixed_index], f_lhd(x) - min_lhd, \
                "MAX STEP:", max_feasible_step_size, "REAL STEP", alpha
        assert f_lhd(x) >= min_lhd
        return x
    
    x = mle_estimate.copy()
    prev_lhd = 1e-10
    
    x = mle_estimate
    prev_x = mle_estimate[fixed_index]
    n_successes = 0
    while True:        
        # take a downhill step
        x, curr_lhd = take_param_decreasing_step(x)
        if bound_type == 'LOWER' \
                and abs(x[fixed_index] - MIN_TRANSCRIPT_FREQ) < 1e-13 \
                and curr_lhd >= min_lhd:
            break
        if bound_type == 'UPPER' \
                and abs(x[fixed_index] - (1.0-MIN_TRANSCRIPT_FREQ)) < 1e-13 \
                and curr_lhd >= min_lhd:
            break
        
        x = take_lhd_decreasing_step(x)
        
        if abs( prev_x - x[fixed_index] ) < 1e-9: 
            n_successes += 1
            if n_successes > 3:
                break
        else: 
            prev_x = x[fixed_index]
            n_successes = 0
        
        

    rv = chi2.sf( 2*(max_lhd-f_lhd(x)), 1), x[fixed_index]
    return rv
    
    final_x = numpy.ones(n)*MIN_TRANSCRIPT_FREQ
    final_x[ numpy.array(sorted(set(range(n))-zeros)) ] = x
    final_lhd = calc_lhd(final_x, observed_array, full_expected_array)
    assert final_lhd >= f_lhd(x) - abs_tol
    return final_x, lhds


def estimate_confidence_bound_CVXOPT( observed_array, 
                                      expected_array, 
                                      fixed_index,
                                      mle_estimate,
                                      bound_type,
                                      alpha = 0.025):
    global DEBUG_OPTIMIZATION
    DEBUG_OPTIMIZATION = True
    rv = estimate_confidence_bounds_directly( 
        observed_array,  expected_array, fixed_index,
        calc_lhd(mle_estimate, observed_array, expected_array), 
        (bound_type=='UPPER'), alpha )
    print "DIRECT APPROACH", bound_type, rv
    assert False
