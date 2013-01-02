import sys, os
import numpy
from itertools import izip

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

def estimate_transcript_frequencies( observed_array, expected_array ):
    ps = matrix(expected_array)
    thetas = variable( expected_array.shape[1] )
    uniform_theta_value = 1.0/expected_array.shape[1]
    Xs = matrix( observed_array )

    # ridge_lambda = Xs.sum()/100
    # -ridge_lambda*quad_form(thetas,1)
    # -ridge_lambda*sum(square(thetas))),
    p = program( maximize(Xs*log(ps*thetas)),
                 [eq(sum(thetas), 1), geq( thetas, 1e-10 )] )
    # +1*quad_over_lin(1,thetas[1, 0]) ), 
    # p = program( minimize( -Xs*log(ps*thetas) ),
    #             [eq( sum(thetas), 1), geq( thetas, 0 )] )

    p.options['maxiters']  = 500
    log_lhd = -p.solve(quiet=not VERBOSE)
    thetas_values = thetas.value
    
    return thetas_values.T.tolist()[0]

def estimate_gene_expression( gene, bam_fname, fl_dists ):
    expected_array, observed_array = build_design_matrices( 
        gene, bam_fname, fl_dists )
    
    transcript_frequencies = estimate_transcript_frequencies( 
        observed_array, expected_array )
    
    return 

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
        fl_dists = build_normal_density( fl_min, fl_max, mean, sd )
        read_group_mappings = []
    else:
        fl_dists, read_group_mappings = load_fl_dists( args.fl_dists )
        
    global VERBOSE
    VERBOSE = args.verbose
    
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
    
    return args.gtf, bam_fns, args.ofname, fl_dists, read_group_mappings

if __name__ == "__main__":
    # Get file objects from command line
    gtf_fp, bam_fns, ofname, fl_dists, rg_mappings = parse_arguments()
    
    genes = load_gtf( gtf_fp.name )
    for gene in genes:
        for bam_fn in bam_fns:
            estimate_gene_expression( gene, bam_fn, fl_dists )
    

    sys.exit()























print "done"
