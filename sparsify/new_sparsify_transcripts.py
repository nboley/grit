import sys, os
import numpy
import math
from scipy.linalg import svd
from scipy.stats import chi2
from scipy.optimize import fmin_slsqp

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

MIN_TRANSCRIPT_FREQ = 1e-6
COMPARE_TO_DCP = False
num_threads = 1

DEBUG_OPTIMIZATION = False

"""
START LOCATION STUFF

    def f_init_taylor(xs):
        x_full = [ max( MIN_TRANSCRIPT_FREQ, x ) for x in xs ]
        x_full.insert( fixed_i, MIN_TRANSCRIPT_FREQ )
        log_lhd = calc_taylor_lhd( x_full, observed_array, expected_array )
        if DEBUG_VERBOSE: print x_full, log_lhd
        assert not numpy.isnan( log_lhd )
        # we add the extra penalty to force the objective away from 0
        positive_penalty = sum([ x-y for x,y in zip(x_bounded, xs.tolist())])
        sum_to_one_penalty = abs(1-numpy.absolute(xs).sum())
        return -log_lhd + 100*( positive_penalty + sum_to_one_penalty )

    def f_init(xs):
        x_bounded = [ max( MIN_TRANSCRIPT_FREQ, x ) for x in xs ]
        x_full = x_bounded
        x_full.insert( fixed_i, MIN_TRANSCRIPT_FREQ )
        log_lhd = calc_lhd( x_full, observed_array, expected_array )
        positive_penalty = sum([ x-y for x,y in zip(x_bounded, xs.tolist())])
        sum_to_one_penalty = abs(1-numpy.absolute(xs).sum())
        return -log_lhd + 100*( positive_penalty + sum_to_one_penalty )
    


        x0 = [ x+float(mle_est[fixed_i])/(len(mle_est)-1) for i, x in enumerate(mle_est) if i != fixed_i ]
        bounds = [(MIN_TRANSCRIPT_FREQ,1)]*len(x0) 

        res = fmin_slsqp( f_init, x0=x0, bounds=bounds,
                          f_eqcons=sum_to_one_constraint, 
                          disp=DEBUG_VERBOSE, epsilon=1e-2,
                          full_output=True)    
        # check the error code
        if res[3] != 0:
            log_warning( "Invalid Optimization Return Code From Lower Bnd X0 Search Taylor Objective. Proceeding.\n" + "\n".join( map( str,res )[-2:] ) + "\n" )
        x0 = res[0]
        
        res = fmin_slsqp( f_init, x0=x0, bounds=bounds,
                          f_eqcons=sum_to_one_constraint, 
                          disp=DEBUG_VERBOSE, epsilon=1e-3,
                          full_output=True)    
        
        # check the error code
        if res[3] != 0:
            log_warning( "Invalid Optimization Return Code From Lower Bnd X0 Search Objective. Proceeding.\n" + "\n".join( map( str,res )[-2:] ) + "\n" )
        x0 = res[0]
        
        x0 = res[0]
        x0.insert( fixed_i, MIN_TRANSCRIPT_FREQ )

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
    expected_mat = expected_mat/expected_mat.sum(0)
    
    observed_mat = numpy.array( observed_mat )
    
    return expected_mat, observed_mat

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
    
    read_groups_and_read_lens =  { (RG, read_len) for RG, read_len, bin 
                                   in binned_reads.iterkeys() }
    
    fl_dists_and_read_lens = [ (fl_dists[RG], read_len) for read_len, RG  
                               in read_groups_and_read_lens ]
    
    expected_cnts = calc_expected_cnts( 
        exon_boundaries, transcripts_non_overlapping_exon_indices, 
        fl_dists_and_read_lens)
    
    observed_cnts = build_observed_cnts( binned_reads, fl_dists )
    
    expected_array, observed_array = build_expected_and_observed_arrays( 
        expected_cnts, observed_cnts )
    
    return expected_array, observed_array

def find_convex_hull( expected_array ):
    #expected_array = matrix( ((1,0,0.5), (0,1,0.5), (1,1,1)) )
    if DEBUG_VERBOSE: print "Expected Array:\n", expected_array
    
    zeros = numpy.matrix(numpy.zeros(expected_array.shape[0])).T
    expected_array = numpy.hstack( (expected_array, zeros )  )
    if DEBUG_VERBOSE: print "Augmented Expected Array:\n", expected_array
    #print svd( expected_array )

    from scipy.linalg import orth
    basis = orth(expected_array)
    if DEBUG_VERBOSE: print "Basis:\n", basis
    rotated_expected_array = (basis.T*numpy.matrix(expected_array)).T
    if DEBUG_VERBOSE: print "Rotated Expected:\n", rotated_expected_array
    
    dimension = basis.shape[1]
    num_points = expected_array.shape[1]
    qconvex_input = [ str(dimension), str(num_points)]
    for vector in rotated_expected_array.tolist():
        qconvex_input.append( "\t" + "\t".join( map(str, vector) ) )
    
    qconvex_input = "\n".join( qconvex_input )
    fp = open( "tmp.txt", "w" )
    fp.write(qconvex_input)
    fp.close()

    p = subprocess.Popen( ["qconvex", "Fx"], stdin=subprocess.PIPE, 
                          stderr=subprocess.STDOUT, stdout=subprocess.PIPE )
    output = p.communicate( input=qconvex_input  )[0].split("\n")
    if DEBUG_VERBOSE: print "Ouput:\n", "\n".join(output)
    indices = sorted(map(int, output[-(num_points+1):-1] ))[:-1]
    if DEBUG_VERBOSE: print "Indices:", indices
    
    return indices

def calc_lhd( freqs, observed_array, expected_array ):
    return float(observed_array*numpy.log( numpy.matrix( expected_array )*numpy.matrix(freqs).T ))

def calc_lhd_deriv( freqs, observed_array, expected_array ):
    rv = numpy.zeros( expected_array.shape[1] )
    denom = numpy.matrix( expected_array )*numpy.matrix(freqs).T
    for i in xrange( expected_array.shape[0] ):
        for j in xrange( expected_array.shape[1] ):
            rv[j] += observed_array[j]*expected_array[i,j]/denom[j]
    
    return rv.tolist()

def calc_taylor_lhd( freqs, observed_array, expected_array ):
    bin_prbs = numpy.matrix( expected_array )*numpy.matrix(freqs).T
    return float(observed_array*(-numpy.power(bin_prbs-1.14, 20)))

def calc_lhd_for_subprocess( args ):
    freqs, observed_array, expected_array = args
    return calc_lhd( freqs, observed_array, expected_array )

def estimate_transcript_frequencies_wth_dcp( observed_array, expected_array ):
    #convex_hull_indices = find_convex_hull( expected_array )
    convex_hull_indices = range( expected_array.shape[1] )
    
    Xs = cvxpy.matrix( observed_array )
    ps = cvxpy.matrix( expected_array[:,convex_hull_indices] )
    num_transcripts = ps.shape[1]    
    
    # add the sum to one condition implicitly
    thetas = cvxpy.variable( num_transcripts-1 )
    loss_fn = Xs*(cvxpy.log( ps*cvxpy.vstack((thetas, 1-cvxpy.sum(thetas))) ))
    
    p = cvxpy.program( cvxpy.maximize(loss_fn), 
                       [cvxpy.geq(thetas,0), cvxpy.leq(cvxpy.sum(thetas),1)] )
    p.options['maxiters']  = 1500
    log_lhd = p.solve(quiet=not DEBUG_VERBOSE)
    
    freq_estimates = [0]*num_transcripts
    for index, value in zip(convex_hull_indices, thetas.value.T.tolist()[0]):
        freq_estimates[index] = value
    
    freq_estimates[-1] = 1 - sum(freq_estimates)
    if VERBOSE: print log_lhd, freq_estimates, calc_lhd( 
        freq_estimates, observed_array, expected_array )
        
    return log_lhd, freq_estimates

def estimate_transcript_frequencies( observed_array, expected_array, 
                                     fixed_indices=[], fixed_values=[],
                                     minimum_eps=1e-8, x0=None ):
    assert len( fixed_values ) == len( fixed_indices )
    assert sum( fixed_values ) < 1
    
    num_transcripts = expected_array.shape[1] - len( fixed_indices )
    bounds = [(MIN_TRANSCRIPT_FREQ,1)]*num_transcripts
    indices_and_values = sorted(zip(fixed_indices, fixed_values))
    eq_value = 1. - sum( fixed_values )
    uniform_x0 = [eq_value/num_transcripts]*num_transcripts
    
    def adjust_estimate(x):
        """Ensure the estimate is within the bounds

        """
        x_full = [ max( MIN_TRANSCRIPT_FREQ, x ) for x in x ]
        x_full = [ x*(eq_value/sum(x_full)) for x in x_full ]
        for i, val in indices_and_values:
            x_full.insert( i, val )
        
        assert abs(1. - sum(x_full)) < 1e-9
        return x_full
        
    def boundary_loss( xs ):
        loss_0 = sum( int(x<MIN_TRANSCRIPT_FREQ)*(MIN_TRANSCRIPT_FREQ-x)
                      for x in xs )
        loss_1 = sum( int(x > 1.)*(1-x)
                      for x in xs )
        loss_sum = abs(1-sum(xs))
        return 100*(loss_0 + loss_1 + loss_sum)
    
    def f_taylor(xs):
        x_full = adjust_estimate( xs )
        log_lhd = calc_taylor_lhd( x_full, observed_array, expected_array )
        return -log_lhd + boundary_loss(xs)
    
    def f(xs):
        x_full = adjust_estimate( xs )
        log_lhd = calc_lhd( x_full, observed_array, expected_array )
        return -log_lhd + boundary_loss(xs)
    
    def eq_const( x ):
        return numpy.array( numpy.matrix( eq_value - x.sum() ) )
    
    def minimize( x0, f, epsilon ):
        res = fmin_slsqp( f, x0=x0, 
                          bounds=bounds,
                          f_eqcons=eq_const, 
                          disp=2*int(DEBUG_OPTIMIZATION),
                          full_output=True, 
                          epsilon=epsilon,
                          iter=100)
        
        return res[0], res[-2]

    def find_x0():
        """Find x0 using the polynomial approximation to log.

        """
        epsilon_init = 1e-1
        while epsilon_init > 1e-8:
            x0, rv = minimize( uniform_x0, f_taylor, epsilon_init )
            if rv == 0: 
                return x0
            epsilon_init /= 10
        
        return None

    def find_x(x0):
        freq_estimates = None
        estimate_epsilon = None
        estimate_lhd = None
        for epsilon in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8 ]:
            if epsilon < minimum_eps: break
            x0, rv = minimize( x0, f, epsilon )
            x0_adj = adjust_estimate(x0)
            lhd = calc_lhd( x0_adj, observed_array, expected_array )
            if rv == 0 or lhd < estimate_lhd: 
                freq_estimates = x0_adj
                estimate_lhd = lhd
                estimate_epsilon = epsilon

        if freq_estimates == None:
            raise ValueError, "Couldn't find an optimum."
        
        return freq_estimates

    if x0 == None:
        x0 = find_x0()
        if x0 == None:
            print "WARNING: Couldn't find a starting location."
            x0 = uniform_x0
    
    freq_estimates = find_x( x0 )
    
    if abs(1 - sum(freq_estimates)) > 1e-2:
        print freq_estimates, sum(freq_estimates)
        raise ValueError, "Constraints not satisfied."
    
    freq_estimates = [ x/sum(freq_estimates) for x in freq_estimates ]
    
    log_lhd = calc_lhd( freq_estimates, observed_array, expected_array )
    
    if COMPARE_TO_DCP:
        dcp_lhd, dcp_est = estimate_transcript_frequencies_wth_dcp(
            observed_array, expected_array)

        print "Lhd:", log_lhd, dcp_lhd
        for x, y in zip( freq_estimates, dcp_est ):
            print x, y
    
    return log_lhd, freq_estimates

def estimate_confidence_bound_wth_dcp( observed_array, expected_array, 
                                       mle_log_lhd, fixed_i, upper_bound=True, 
                                       alpha=0.10 ):
    lower_lhd_bound = mle_log_lhd - chi2.ppf( 1 - alpha, 1 )/2.
    free_indices = set(range(expected_array.shape[1])) - set((fixed_i,))
    
    Xs = cvxpy.matrix( observed_array )
    ps = cvxpy.matrix( expected_array )
    thetas = cvxpy.variable( ps.shape[1] )
    
    constraints = [ cvxpy.geq(Xs*cvxpy.log(ps*thetas), lower_lhd_bound), 
                    cvxpy.eq(cvxpy.sum(thetas), 1), 
                    cvxpy.geq(thetas, MIN_TRANSCRIPT_FREQ) ]
    
    if upper_bound:
        p = cvxpy.program( cvxpy.maximize(thetas[fixed_i,0]), constraints )    
    else:
        p = cvxpy.program( cvxpy.minimize(thetas[fixed_i,0]), constraints )
    
    p.options['maxiters']  = 1500
    value = p.solve(quiet=not DEBUG_VERBOSE)
    
    thetas_values = thetas.value.T.tolist()[0]
    log_lhd = calc_lhd( thetas_values, observed_array, expected_array )
    
    return chi2.sf( 2*(mle_log_lhd-log_lhd), 1), value, thetas_values

def estimate_confidence_bound_by_bisection( observed_array, expected_array,
                                            fixed_i, mle_est, 
                                            bound_type,
                                            alpha=0.10 ):
    def etf_wrapped( x ):
        return estimate_transcript_frequencies( 
            observed_array, expected_array, 
            [fixed_i,], fixed_values=[x,], 
            minimum_eps=1e-7,
            x0=[ x for i, x in enumerate(mle_est) if i != fixed_i] )[0]
    
    def bisect_for_bnd( lower_bnd, upper_bnd, MAX_ITER=100, PRECISION=1e-4 ):
        for i in xrange( MAX_ITER ):
            curr = (upper_bnd + lower_bnd)/2
            curr_lhd = etf_wrapped( curr )
            try:
                assert curr_lhd <= mle_log_lhd
            except:
                print curr_lhd, mle_log_lhd
                raise
            
            if DEBUG_VERBOSE:
                print curr, curr_lhd, upper_lhd_bound-curr_lhd, \
                    bound_type, mle_log_lhd
            
            if upper_bnd - lower_bnd < PRECISION:
                return (upper_bnd + lower_bnd)/2., curr_lhd
            
            if curr_lhd < upper_lhd_bound:
                if bound_type=='UPPER': upper_bnd = curr
                else: lower_bnd = curr
            else:
                if bound_type=='UPPER': lower_bnd = curr
                else: upper_bnd = curr

        
        raise ValueError, "Bisection failed to find confidence bound."
    
    mle_bnd = mle_est[fixed_i]
    mle_log_lhd = calc_lhd( mle_est, observed_array, expected_array )
    upper_lhd_bound = mle_log_lhd - chi2.ppf( 1 - alpha, 1 )/2.
    
    if bound_type=='UPPER': other_bnd = 1.0 - MIN_TRANSCRIPT_FREQ
    else: other_bnd = MIN_TRANSCRIPT_FREQ
    
    lower_bnd, upper_bnd = min( mle_bnd, other_bnd ), max( mle_bnd, other_bnd )
    bnd, log_lhd = bisect_for_bnd( lower_bnd, upper_bnd )
    
    if COMPARE_TO_DCP:
        dcp_est = estimate_confidence_bound_wth_dcp( 
            observed_array, expected_array,
            mle_log_lhd, fixed_i, upper_bound, alpha )
        dcp_lhd = calc_lhd( dcp_est[2], observed_array, expected_array )
        
        print "DCP_COMPARISON", upper_bound, freq_estimates[fixed_i], \
            dcp_est[1], log_lhd, dcp_lhd, upper_lhd_bound
    
    assert log_lhd >= upper_lhd_bound - 0.1
    
    return chi2.sf( 2*(mle_log_lhd-log_lhd), 1), bnd

class MaxIterError( ValueError ):
    pass

def estimate_confidence_bounds_directly( observed_array, 
                                         expected_array, 
                                         fixed_index,
                                         mle_estimate,
                                         bound_type,
                                         alpha = 0.025):
    def adjust_estimate(x):
        """Ensure the estimate is within the bounds

        """
        x_full = [ max( MIN_TRANSCRIPT_FREQ, x ) for x in x ]
        x_full = [ x*(eq_value/sum(x_full)) for x in x_full ]
        assert abs(1. - sum(x_full)) < 1e-9
        return x_full
        
    def boundary_loss( xs ):
        loss_0 = sum( int(x<MIN_TRANSCRIPT_FREQ)*(MIN_TRANSCRIPT_FREQ-x)
                      for x in xs )
        loss_1 = sum( int(x > 1.)*(1-x)
                      for x in xs )
        loss_sum = abs(1-sum(xs))
        return 100*(loss_0 + loss_1 + loss_sum)
        
    def lhd_ratio_constraint(xs):
        x_full = adjust_estimate( xs )
        log_lhd = calc_lhd( x_full, observed_array, expected_array )
        return numpy.array(numpy.matrix(log_lhd + boundary_loss(xs) - upper_lhd_bound))
    
    def eq_const( x ):
        return numpy.array( numpy.matrix( eq_value - x.sum() ) )
    
    def f( xs ):
        return (2*int(bound_type == "LOWER")-1)*xs[fixed_index]
    
    def minimize( x0, f, epsilon, maxiter ):
        bounds = [(MIN_TRANSCRIPT_FREQ,1)]*num_transcripts
        res = fmin_slsqp( f, x0=x0, 
                          bounds=bounds,
                          f_eqcons=eq_const, 
                          f_ieqcons=lhd_ratio_constraint,
                          disp=2*int(DEBUG_OPTIMIZATION),
                          full_output=True, 
                          epsilon=epsilon,
                          iter=maxiter)
        
        return res[0], res[-2]
    
    def find_x(x0, maxiter):
        maxiter_reached = False
        freq_estimates = None
        estimate_epsilon = None
        estimate_lhd = None
        for epsilon in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8 ]:
            x0, rv = minimize( x0, f, epsilon, maxiter )
            x0_adj = adjust_estimate(x0)
            lhd = calc_lhd( x0_adj, observed_array, expected_array )
            if rv == 9:
                maxiter_reached = True
            if rv == 0: 
                freq_estimates = x0_adj
                estimate_lhd = lhd
                estimate_epsilon = epsilon

        if freq_estimates == None:
            if maxiter_reached:
                raise MaxIterError, "Couldn't find an optimum."
            else:
                raise ValueError, "Couldn't find an optimum."
        
        return freq_estimates

    num_transcripts = expected_array.shape[1]
    eq_value = 1.
    x0 = mle_estimate
    mle_log_lhd = calc_lhd( mle_estimate, observed_array, expected_array )
    upper_lhd_bound = mle_log_lhd - chi2.ppf( 1 - alpha, 1 )/2.

    # special case estiamtes already at the boundaries
    if bound_type == 'LOWER' and \
            mle_estimate[fixed_index]-MIN_TRANSCRIPT_FREQ < MIN_TRANSCRIPT_FREQ:
        return 1, MIN_TRANSCRIPT_FREQ
    elif bound_type == 'UPPER' and \
            1-mle_estimate[fixed_index] < MIN_TRANSCRIPT_FREQ:
        return 1, 1.0
    
    for maxiter in ( 100, 1000 ):
        try:
            freq_estimates = find_x( x0, maxiter )
        except MaxIterError:
            if DEBUG_VERBOSE:
                print "WARNING: max iter reached. Increasing to 1000."
            continue
        else:
            break
    
    if abs(1 - sum(freq_estimates)) > 1e-2:
        print freq_estimates, sum(freq_estimates)
        raise ValueError, "Constraints not satisfied."
    
    log_lhd = calc_lhd( freq_estimates, observed_array, expected_array )
    estimate = freq_estimates[fixed_index]
    
    if DEBUG_VERBOSE:
        print estimate, log_lhd, upper_lhd_bound - log_lhd, mle_log_lhd
    
    return chi2.sf( 2*(mle_log_lhd-log_lhd), 1), estimate

def estimate_confidence_bound( observed_array, 
                               expected_array, 
                               fixed_index,
                               mle_estimate,
                               bound_type,
                               alpha = 0.025):
    try:
        return estimate_confidence_bounds_directly( 
            observed_array,  expected_array, fixed_index,
            mle_estimate, bound_type, alpha )
    except ValueError:
        print "WARNING: Couldn't find bound for transcript %i %s directly: trying bisection." % (fixed_index+1, bound_type)
        return estimate_confidence_bound_by_bisection( 
            observed_array,  expected_array, fixed_index,
            mle_estimate, bound_type )
                    
    pass

def build_grid( expected_array, observed_array ):
    grid = []
    n = 10
    for i in range( n+1 ):
        for j in range( n+1 ):
            for k in range( n+1 ):
                entry = [i/float(n), j/float(n), k/float(n)]
                if sum( entry ) > 1: continue
                entry.append( round(1 - sum( entry ),1) )
                grid.append( entry  )
    
    from multiprocessing import Pool
    p = Pool( 50 )
    res = p.map(calc_lhd, [ (x, observed_array, expected_array) for x in grid])

    for entry, lhd in zip( grid, res ):
        print "\t".join(map(str, entry)) + "\t" + str(lhd)
        
    sys.exit()


def estimate_gene_expression( gene, bam_fname, fl_dists, estimate_confidence_bounds=False ):
    expected_array, observed_array = build_design_matrices( 
        gene, bam_fname, fl_dists )
    
    log_lhd, mle_estimate = estimate_transcript_frequencies( 
        observed_array, expected_array )
    if VERBOSE: print gene.id, log_lhd, [ "%.2e" % x for x in mle_estimate ]
        
    if estimate_confidence_bounds:
        bnds = []
        for index, mle_value in enumerate( mle_estimate ):
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
        ofps[bam_fn] = ThreadSafeFile("%s.%s.gtf" % (args.ofprefix, os.path.basename(bam_fn)), "w")
        ofps[bam_fn].write("track name=transcripts.%s useScore=1\n" % os.path.basename(bam_fn))
        
    return args.gtf, bam_fns, ofps, fl_dists, read_group_mappings, args.estimate_confidence_bounds

def estimate_gene_expression_worker(input_queue, fl_dists, ofps, 
                                    estimate_confidence_bounds=False, 
                                    filter_value=1e-3):
    while not input_queue.empty():
        gene, bam_fn = input_queue.get()
        if VERBOSE: print "Processing gene %s sample %s." % ( gene.id, bam_fn )
        mles, lbs, ubs = estimate_gene_expression( gene, bam_fn, fl_dists )
        for mle, lb, ub, transcript in zip(mles, lbs, ubs,gene.transcripts):
            if mle/max( mles ) < filter_value: continue
            meta_data = { "frac": "%.2e" % mle }
            transcript.score = max( 1, int(1000*mle/max(mles)) )
            if estimate_confidence_bounds:
                meta_data["conf_lo"] = "%.2e" % lb
                meta_data["conf_hi"] = "%.2e" % ub
            
            ofps[bam_fn].write( transcript.build_gtf_lines(
                    gene.id, meta_data, source="grit") + "\n" )
        if VERBOSE: print "Finished gene %s sample %s." % ( gene.id, bam_fn )

def main():
    # Get file objects from command line
    gtf_fp, bam_fns, ofps, fl_dists, rg_mappings, estimate_confidence_bounds \
        = parse_arguments()

    manager = multiprocessing.Manager()
    input_queue = manager.Queue()    
    genes = load_gtf( gtf_fp.name )
    for gene in genes:
        for bam_fn in bam_fns:
            input_queue.put( (gene, bam_fn) )
    
    ps = []
    for i in xrange( min( num_threads, len(genes)*len(bam_fns) ) ):
        p = Process(target=estimate_gene_expression_worker, 
                    args=(input_queue, fl_dists, ofps, estimate_confidence_bounds))
        p.start()
        ps.append( p )
    
    for p in ps:
        p.join()

    for ofp in ofps.values():
        ofp.close()
    
    return

if __name__ == "__main__":
    main()
