import sys, os
import numpy
from scipy.linalg import svd
from scipy.stats import chi2

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

from cvxpy import maximize, minimize, geq, eq, variable, matrix, program, log, \
    sum, quad_form, square, quad_over_lin, geo_mean

num_threads = 1

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
    return (0,1,2,3)
    sys.exit()
    
    dimension = expected_array.shape[0]
    num_points = expected_array.shape[1]
    qconvex_input = [ str(dimension), str(num_points)]
    for vector in expected_array.T.tolist():
        for i in xrange( 8 ):
            qconvex_input.append( "\t" + "\t".join( map(str, vector) ) )
    
    qconvex_input = "\n".join( qconvex_input )
    fp = open( "tmp.txt", "w" )
    fp.write(qconvex_input)
    fp.close()
    sys.exit()
    p = subprocess.Popen( ["qconvex", "Fx"], stdin=subprocess.PIPE, 
                          stderr=subprocess.STDOUT, stdout=subprocess.PIPE )
    indices = p.communicate( input=qconvex_input  )[0]
    
    sys.exit()

def calc_lhd( freqs, observed_array, expected_array ):
    return float(observed_array*log( matrix( expected_array )*matrix(freqs).T ))

def calc_lhd_for_subprocess( args ):
    freqs, observed_array, expected_array = args
    return calc_lhd( freqs, observed_array, expected_array )

def estimate_transcript_frequencies( observed_array, expected_array ):
    convex_hull_indices = range(expected_array.shape[1]) 
    # find_convex_hull( expected_array )
    
    Xs = matrix( observed_array )
    ps = matrix( expected_array[:,convex_hull_indices] )
    thetas = variable( ps.shape[1] )
    
    p = program( maximize(Xs*log(ps*thetas)), 
                 [eq(sum(thetas), 1), geq(thetas,0)])
    
    p.options['maxiters']  = 1500
    log_lhd = p.solve(quiet=not EXTRA_VERBOSE)
    
    freq_estimates = [0]*expected_array.shape[1]
    for index, value in zip(convex_hull_indices, thetas.value.T.tolist()[0]):
        freq_estimates[index] = value
    
    return log_lhd, freq_estimates

def estimate_confidence_bound( observed_array, expected_array, mle_log_lhd, 
                               fixed_i, upper_bound=True, alpha=0.10 ):
    lower_lhd_bound = mle_log_lhd - chi2.ppf( 1 - alpha, 1 )/2.
    
    free_indices = set(range(expected_array.shape[1])) - set((fixed_i,))
    
    Xs = matrix( observed_array )
    ps = matrix( expected_array )
    thetas = variable( ps.shape[1] )
    
    constraints = [ geq(Xs*log(ps*thetas), lower_lhd_bound), 
                    eq(sum(thetas), 1), geq(thetas, 1e-6)]
    
    if upper_bound:
        p = program( maximize(thetas[fixed_i,0]), constraints )    
    else:
        p = program( minimize(thetas[fixed_i,0]), constraints )
    
    p.options['maxiters']  = 1500
    value = p.solve(quiet=not EXTRA_VERBOSE)
    
    thetas_values = thetas.value.T.tolist()[0]
    log_lhd = calc_lhd( thetas_values, observed_array, expected_array )
    
    return chi2.sf( 2*(mle_log_lhd-log_lhd), 1), value

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


def estimate_gene_expression( gene, bam_fname, fl_dists ):
    expected_array, observed_array = build_design_matrices( 
        gene, bam_fname, fl_dists )
    
    log_lhd, mle_estimate = estimate_transcript_frequencies( 
        observed_array, expected_array )
    if VERBOSE: print gene.id, log_lhd, [ "%.2e" % x for x in mle_estimate ]

    bnds = []
    for index, mle_value in enumerate( mle_estimate ):
        bnds.append( [] )
        for upper in ( False, True ):
            p_value, bnd = estimate_confidence_bound( 
                observed_array, expected_array, log_lhd, index, upper )
            bnds[-1].append( round( bnd, 6 ) )
        if VERBOSE: 
            print "Gene %s\tTranscript %i\t\tBnds %s" % (gene.id, index+1, bnds[-1])
    
    lower_bnds, upper_bnds = zip( *bnds )
        
    return mle_estimate, lower_bnds, upper_bnds

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='Determine valid transcripts and estimate frequencies.')
    parser.add_argument( 'ofname', \
        help='Output file name')
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

    parser.add_argument( '--write-meta-data', '-m', default=False, 
        action='store_true', help='Whether or not to write out meta data.')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
                             help='Whether or not to print status information.')
    parser.add_argument( '--extra-verbose', default=False, action='store_true',\
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
            raise ValueError, "Mean and SD for normal fl_dist are not properly formatted."
        
        mean, sd = fl_dist_norm
        fl_min = max( 0, mean - (4 * sd) )
        fl_max = mean + (4 * sd)
        fl_dists = { 'mean': build_normal_density( fl_min, fl_max, mean, sd ) }
        read_group_mappings = []
    else:
        fl_dists, read_group_mappings = load_fl_dists( args.fl_dists )
    
    global EXTRA_VERBOSE
    EXTRA_VERBOSE = args.extra_verbose
    
    global VERBOSE
    VERBOSE = ( args.verbose or EXTRA_VERBOSE )

    # we change to the output directory later, and these files need to opened in
    # each sub-process for thread safety, so we get the absokute path while we can.
    bam_fns = [ os.path.abspath( bam_fn ) for bam_fn in args.bam_fns ]

    global PROCESS_SEQUENTIALLY
    if args.threads == 1:
        PROCESS_SEQUENTIALLY = True
    
    global WRITE_META_DATA
    WRITE_META_DATA = args.write_meta_data
    
    global num_threads
    num_threads = args.threads
    
    ofp = ThreadSafeFile( args.ofname, "w" )
    
    return args.gtf, bam_fns, ofp, fl_dists, read_group_mappings

def estimate_gene_expression_worker(input_queue, fl_dists, ofp):
    while not input_queue.empty():
        gene, bam_fn = input_queue.get()
        if VERBOSE: print "Processing gene %s." % gene.id
        mles, lbs, ubs = estimate_gene_expression( gene, bam_fn, fl_dists )
        for mle, lb, ub, transcript in zip(mles, lbs, ubs,gene.transcripts):
            meta_data = { "frac": mle, "conf_lo": "%.2f" % lb, 
                          "conf_hi" : "%.2f" % ub }
            ofp.write( transcript.build_gtf_lines(
                    gene.id, meta_data, source="grit") + "\n" )
        if VERBOSE: print "Finished gene %s." % gene.id

def main():
    # Get file objects from command line
    gtf_fp, bam_fns, ofp, fl_dists, rg_mappings = parse_arguments()

    manager = multiprocessing.Manager()
    input_queue = manager.Queue()    
    genes = load_gtf( gtf_fp.name )
    for gene in genes:
        for bam_fn in bam_fns:
            input_queue.put( (gene, bam_fn) )
    
    ps = []
    for i in xrange( min( num_threads, len(genes) ) ):
        p = Process(target=estimate_gene_expression_worker, 
                    args=(input_queue, fl_dists, ofp))
        p.start()
        ps.append( p )
    
    for p in ps:
        p.join()
    
    ofp.close()
    
    return

if __name__ == "__main__":
    main()
