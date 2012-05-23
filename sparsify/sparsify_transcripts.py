#!/usr/bin/python
import sys
import os

sys.path.append( os.path.join( os.path.dirname(__file__), "..", "scikit-learn" ) )
from sklearn import linear_model

import subprocess
import shutil
import re
import numpy
import multiprocessing
import time
import itertools
import traceback
import datetime
import signal
import Queue

from math import log, sqrt, exp
from collections import defaultdict

from frag_len import FlDist

MAX_NUM_EXONS = 100

FILTER_BY_IMPROBABILITY = False
FREQ_FILTER = 0.01
FREQ_FILTER_THRESH = 0.05

PRINT_READ_PAIRS = False
PAUSE_ON_ERROR = False
DO_PROFILE = False
VERBOSE = False
MINIMAL_VERBOSE = False

DEBUG = False
DEBUG_VERBOSE = False
WRITE_META_DATA = False

# when debugging, we sometimes want to do everything from 
# the main thread. This builds all of the output, and then
# processes it all sequentially, so that we can get debugging 
# info.
PROCESS_SEQUENTIALLY = ( DO_PROFILE or False )

# if an exon boundry is also a junction edge,
# then it can't start or stop a transcript. ( BC 
# splice sites are very different from transcription start sites/
# transcription stop sites )
FILTER_EXONS_THAT_ARE_ALSO_JNS = True
MINIMUM_INTRON_SIZE = 20
LONG_SINGLE_EXON_GENES_FNAME = "high_read_depth_exons.gtf"


import reads
reads.PRINT_READ_PAIRS = PRINT_READ_PAIRS
reads.DEBUG_VERBOSE = DEBUG_VERBOSE
reads.MINIMUM_INTRON_SIZE = MINIMUM_INTRON_SIZE
reads.VERBOSE = VERBOSE
import transcripts
transcripts.VERBOSE = VERBOSE

from frag_len import build_normal_density, load_fl_dists
from reads import Reads, BinnedReads
from gene_models import parse_gff_line, GeneBoundaries
from transcripts import *
from f_matrix import *

class ThreadSafeFile( file ):
    def __init__( *args ):
        file.__init__( *args )
        self.lock = multiprocessing.Lock()

    def write( self, string ):
        self.lock.acquire()
        file.write( self, string )
        self.lock.release()


def get_memory_usage( pid ):
    cmd = "ps -o%%mem --pid %i | awk '{sum+=$1} END {print sum}'" % pid 
    process = subprocess.Popen( cmd, shell=True, stdout=subprocess.PIPE )
    return float( process.communicate()[0].strip() )

def get_tracking_filename( bam_fn, gene_name ):
    # strip off the last four characters to get rid of the .bam
    assert bam_fn[-4:] == '.bam'
    source = os.path.basename( bam_fn )[:-4]
    return os.path.join( "./tracking_files", source, \
                             source + "." + gene_name + ".tracking" )


def estimate_freqs_with_nnls( expected_cnts, observed_cnts, indices ):
    """re-estimate using unbiased estimator
    
    """
    filtered_expected_cnts = expected_cnts[ :, indices ]
    
    coefs, residuals = optimize.nnls( filtered_expected_cnts, observed_cnts )
    if not any( coefs ):
        raise ValueError
    dp = numpy.dot( filtered_expected_cnts, coefs )
    assert abs(  numpy.sqrt(( ( observed_cnts - dp )**2 ).sum()) -  residuals ) < 1e-6
    freqs = list( coefs/coefs.sum() )
    # insert the 0 freqs where appropriate
    freqs_mapping = dict( zip( indices, freqs ) )
    freqs = []
    for transcript_index in xrange( expected_cnts.shape[1] ):
        if transcript_index in freqs_mapping:
            freqs.append( freqs_mapping[ transcript_index ] )
        else:
            freqs.append( 0.0 )
    
    return numpy.array( freqs ), residuals

def choose_lambda_via_cv( expected_cnts, observed_cnts, lambdas ):
    """Perform k-fold CV on the lasso path.

    """
    # find the bin prbs
    num_obs = observed_cnts.sum()
    bin_freqs = numpy.array( observed_cnts, dtype=float )/num_obs
    
    def decide_to_stop( mses ):
        """Returns the minimum acceptable lambda index if one exists, None otherwise.
        
        """
        # we need at least three values of lambda
        if len( mses ) < 3: return None
        # return lambda if the MSE increases for 2 consecutive steps 
        if mses[-1] > mses[-3] and mses[-2] > mses[-3]:
            return mses.index( min( mses ) )
        else:
            return None
        assert False
    
    def bootstrap_freqs( indices ):    
        training_size = int(num_obs)/2
        training_set = numpy.random.multinomial( training_size, bin_freqs )
        test_set = observed_cnts - training_set
        freqs, residuals = estimate_freqs_with_nnls( \
            expected_cnts, training_set, indices )
        freq_scale = num_obs - training_size
        MSE = sqrt( \
            ( ( test_set - numpy.dot( expected_cnts, freqs*freq_scale ) )**2 ).sum() )
        MSE_training = sqrt( \
            ( ( training_set - \
                    numpy.dot( expected_cnts, freqs*freq_scale ) )**2 ).sum() )
        return freqs, ( MSE, MSE_training, residuals )

    ordered_lambdas = sorted( lambdas, reverse=True )
    
    mse_path = []
    
    if VERBOSE:
        print "\nEstimating Optimal Lambda via Cross Validation:"
        print "Lambda".ljust(20), "MSE ( min, mean, max )"

    prev_lambda = -1.0
    min_mse = 1e100
    # remove zeros from ordered_lambdas
    if 0 in ordered_lambdas:
        first_zero_index = ordered_lambdas.index(0)
        ordered_lambdas = ordered_lambdas[:first_zero_index]
    for last_lambda in ordered_lambdas:
        if  (1 - last_lambda/prev_lambda) < 1e-2:
            continue
        
        prev_lambda = last_lambda
        
        indices = (lambdas >= last_lambda).nonzero()[0]
        
        NUM_SAMPLES = 5
        mses = []
        try:
            for sample_num in xrange( NUM_SAMPLES ):
                quals = bootstrap_freqs( indices )[1]
                mses.append( quals[0] )
        except ValueError:
            if last_lambda == ordered_lambdas[-1]:
                raise ValueError, 'Cannot estimate frequencies when all ' + \
                    'coefficients are zero.'

            continue
        mses = numpy.array( mses )
        min_mse = min( min_mse, mses.mean() )                       
        mse_path.append( ( mses.mean() - mses.std(), mses.mean(), \
                           mses.mean() + mses.std() ) )
        if VERBOSE:
            print ("%e" % last_lambda).ljust( 20 ), \
                  "[ %e, %e, %e ]" % \
                  ( mses.mean() - 1.0*mses.std(), mses.mean(), \
                    mses.mean() + 1.0*mses.std() )
        
        # add an early stop condition
        if mses.mean() > 1.5*min_mse: break

    # typically we would choose the min mean MSE as the value of the 
    # lambda cutoff, however, because there is significant noise in 
    # the MSE estiamtes and we prefer a larger transcript set than a 
    # smaller one ( due to downstream filtering ) we expand the mse to 
    # include any transcripts that are within a std of the min.
    min_mean_mse = min( item[1] for item in mse_path )
    min_mean_mse_index = max( index for index, mse in enumerate( mse_path ) \
                              if mse[1] == min_mean_mse )
    
    min_mean_mse_max_value = mse_path[ min_mean_mse_index ][2]
    for min_mse, mean_mse, max_mse in mse_path[min_mean_mse_index+1:]:
        if min_mse < min_mean_mse_max_value:
            min_mean_mse_index += 1
        else: break
    
    
    if VERBOSE:
        print "Chose the optimal lamda: %e" % ordered_lambdas[ min_mean_mse_index ]
    
    return ordered_lambdas[ min_mean_mse_index ]

def build_lars_path( expected_cnts, observed_cnts ):
    lasso_lambdas, added_predictors, coefs = linear_model.lars_path( \
        expected_cnts, observed_cnts, \
        non_negative=True, verbose=VERBOSE )

    ### reorder the lambdas so that they correspond with the trasncripts' order
    # we do this to keep track of the transcripts with lambda 0.
    lambda_mapping = dict( zip( added_predictors, lasso_lambdas )  )
    lasso_lambdas = []
    for transcript_index in xrange( expected_cnts.shape[1] ):
        if transcript_index in lambda_mapping:
            lasso_lambdas.append( lambda_mapping[ transcript_index ] )
        else:
            lasso_lambdas.append( 0.0 )
    
    return numpy.array( lasso_lambdas )
    

def estimate_stability( expected_cnts, observed_cnts, lasso_lambda, num_samples=100 ):
    all_change_points = [ [] for loop in xrange( expected_cnts.shape[1] ) ]
    cnts = numpy.zeros( expected_cnts.shape[1]  )
    for loop in xrange( num_samples ):
        num_obs = observed_cnts.sum()
        bin_freqs = numpy.array( observed_cnts, dtype=float )/num_obs
        training_set = numpy.random.multinomial( num_obs, bin_freqs )

        """
        lambdas = build_lars_path( expected_cnts, training_set )
        for index in xrange(len(lambdas)):
            all_change_points[ index ].append( lambdas )
        """

        clf = linear_model.Lasso( alpha=lasso_lambda, fit_intercept=False )
        clf.fit( expected_cnts, training_set )
        cnts += numpy.array( clf.coef_, dtype=bool )

    """
    for index in xrange(len(all_change_points)):
        all_change_points[index] = numpy.array( all_change_points[index] )
        #print all_change_points[index]
    """

    return cnts/num_samples

def calc_improbability( f_mats, candidate_transcripts, freqs = None ):
    freqs = None
    if freqs == None:
        freqs = numpy.array([FREQ_FILTER]*len(candidate_transcripts))
    
    #### do more filtering
    # find the expected number of reads for each transcript
    single_bins_fmats = reads.build_nonoverlapping_bin_observations( f_mats )
    
    expected_cnts, observed_cnts = convert_f_matrices_into_arrays( single_bins_fmats )
    
    weights = expected_cnts.sum(0)
    weights = numpy.dot( weights, FREQ_FILTER )
    weights = weights/weights.sum()

    expected_read_cnts = weights*observed_cnts.sum()
    
    improbability_measure = []
    for transcript_index in xrange( expected_cnts.shape[1] ):
        exp_bin_fracs = expected_cnts[ :, transcript_index ] \
                        /expected_cnts[ :, transcript_index ].sum() 

        exp_num_reads = expected_read_cnts[ transcript_index ]
        log_prbs = [ 0, ]
        for cnt, exp_frac in zip( observed_cnts, exp_bin_fracs ):
            if cnt > 0: continue
            if exp_frac == 0.0 or exp_frac == 1.0: continue
            # calculate the probability of observing 0 reads 
            log_prb = exp_num_reads*log( 1 - exp_frac )
            log_prbs.append( log_prb )
        
        improbability_measure.append( exp(min(log_prbs)) )
    
    return improbability_measure
    
    # filter the transcripts by their metadata
    good_indices = set([ i for i in xrange(len(improbability_measure)) \
                             if improbability_measure[i] > FREQ_FILTER_THRESH ])
    improbabilites = [ x for x in improbability_measure if x > FREQ_FILTER_THRESH ]
    trans = list( candidate_transcripts.iter_transcripts_and_metadata() )
    new_trans = [ trans for index, trans in enumerate( trans ) if index in good_indices ]
    return Transcripts( candidate_transcripts.gene, new_trans ), improbabilites

def perform_optimization( f_mats, gene, binned_reads, candidate_transcripts ):
    """

    """
    expected_cnts, observed_cnts = convert_f_matrices_into_arrays( f_mats )
    if observed_cnts.sum() == 0:
        raise GeneProcessingError( gene, binned_reads.reads, "No observations." )
    
    #### filter transcripts    
    lasso_lambdas = build_lars_path( expected_cnts, observed_cnts )    
    optimal_lambda = choose_lambda_via_cv( expected_cnts, observed_cnts, lasso_lambdas )
    
    if False and FILTER_BY_IMPROBABILITY:
        good_indices = numpy.array([ i for i in xrange(len(lasso_lambdas)) \
                                 if lasso_lambdas[i] >= ( optimal_lambda - 1e-12)\
                                    and improbabilities[i] > FREQ_FILTER_THRESH ])
    else:
        good_indices = numpy.array([ i for i in xrange(len(lasso_lambdas)) \
                                 if lasso_lambdas[i] >= ( optimal_lambda - 1e-12)])

    
    if VERBOSE:
        print "\nRe-estimating frequencies on filtered transcripts using " + \
            "an unbiased procedure..."  
    
    assert len( good_indices ) > 0
    freqs, residuals = estimate_freqs_with_nnls( \
        expected_cnts, observed_cnts, good_indices )    
    
    improbabilities = numpy.array( calc_improbability( \
            f_mats, candidate_transcripts, freqs ) )
    
    if False and FILTER_BY_IMPROBABILITY:
        good_indices = (improbabilities > FREQ_FILTER_THRESH).nonzero()[0]

    transcripts = Transcripts( gene )
    for index, ( transcript, lasso_lambda, freq, log_prb ) \
            in enumerate( zip( candidate_transcripts, lasso_lambdas, \
                                   freqs, improbabilities ) ):
        transcripts.add_transcript( transcript, binned_reads.sourcefile, \
                                        lasso_lambda, freq, log_prb=log_prb )
    
    transcripts.sort(order='freq')
    
    return transcripts, {"good_indices": good_indices,     \
                             "optimal_lambda": optimal_lambda, \
                             'improbs': improbabilities }

def estimate_gene_expression( gene, candidate_transcripts, binned_reads, fl_dists, \
                                  read_group_mappings ):
    """ as we filter out the different transcripts, append the objects
    to this list. The first entry is always the un lasso fit transcripts:
    ie, they are only fit using the filtering heuristics. Then we start 
    paring the transcript lists down
    """
    if VERBOSE:
        print "==================NEW GENE %s ====================" % gene.name
    
    # Make sure that we get at least one transcript, 1 exon and that we have valid bins.
    assert len( candidate_transcripts ) > 0
    assert len( gene.exon_bndrys ) > 0
    if len( binned_reads.binned_reads.keys() ) == 0:
        raise GeneProcessingError( gene, binned_reads.reads, "Zero valid bins." )

    ### Estimate transcript freqs
    # build the f matrix
    f_mats = build_f_matrix( candidate_transcripts, binned_reads, gene, fl_dists )
        
    transcripts, meta_data = perform_optimization( \
        f_mats, gene, binned_reads, candidate_transcripts )
    
    # build the new transcripts object
    if VERBOSE:
        print str("\nTranscript:").ljust( 50 ), \
              str(" Lambda:").ljust(15), \
              " Freq:".ljust(15), \
              " Improb:"
        
        for transcript, md in transcripts.iter_transcripts_and_metadata():
            print str(transcript).ljust( 50 ), \
                  ("%e" % md.lasso_lambda).ljust(15), \
                  ( "%s" % md.freq ).ljust( 15 ), \
                  md.improbability
    
    if VERBOSE:
        print "Finished estimating transcript frequencies for", gene.name
        
    return transcripts, meta_data

def build_reads_objs( bam_fns, fl_dists, read_grp_mappings ):
    def some_reads_are_paired(reads):
        for cnt, read in enumerate(reads.fetch()):
            if read.is_paired:
                return True
            if cnt > 1000:
                return False
    
    # create reads objects and store in dict
    reads_objs = {}
    for bam_fn in bam_fns:
        reads = Reads( bam_fn )        
        
        try:
            reads.fetch("X", 0, 1)
        except ValueError:
            print "Warning: %s does not have associated bam index [*.bai] file." % bam_fn
            print "Skipping: ", bam_fn
            continue
        
        if not some_reads_are_paired( reads ):
            print "Warning: SLIDE is not compatible with single end reads."
            print "Skipping: ", bam_fn
            continue
        
        reads.fl_dists = fl_dists
        reads.read_group_mappings = read_group_mappings
        
        reads_objs[ bam_fn ] = reads

    return reads_objs

class GeneProcessingError( Exception ):
    def make_string( self ):
        error_string = '-'*60
        error_string += "\nGene " + self.gene.name + " which has " + \
            str( len( self.gene.exon_bndrys ) ) + \
            " exons has produced an error and will not be processed:\n"
        error_string += "Error occured while processing " + \
            os.path.basename( self.reads.filename ) + " reads file.\n\n"
        error_string += traceback.format_exc()
        error_string += '-'*60
        error_string += "\nError Detail:\n"
        error_string += self.detail
        error_string += "\n\n"
        return error_string
    
    def __init__( self, gene, reads, detail="" ):
        self.gene = gene
        self.reads = reads
        self.detail = detail
        self.error_string = self.make_string()
        if VERBOSE:
            print self.error_string

    def __str__( self ):
        return self.error_string


def make_error_log_string( gene, reads_filename, error_inst ):
    error_fields = []
    error_fields.append( gene.name )
    error_fields.append( str( len( gene.exon_bndrys ) ) )
    error_fields.append( os.path.basename( reads_filename ) )
    if type( error_inst ) == GeneProcessingError:
        error_fields.append( error_inst.detail )
    else:
        error_fields.append( str( error_inst ) )
    
    for index, field in enumerate(error_fields):
        field = field.replace( ',', ';' )
        field = field.replace( '\n', '\t' )
        error_fields[ index ] = field

    error_string = ','.join( error_fields )
    
    if VERBOSE:
        print error_string
        print traceback.format_exc()
    
    return error_string

def process_gene_and_reads( gene, candidate_transcripts, reads, op_transcripts):
    binned_reads = BinnedReads( gene, reads, reads.read_group_mappings )
    
    fl_dists = reads.fl_dists
    read_group_mappings = reads.read_group_mappings
    
    try:
        transcripts, meta_data = estimate_gene_expression( \
            gene, candidate_transcripts, binned_reads, fl_dists, read_group_mappings )
    except Exception, inst:
        return make_error_log_string( gene, reads.filename, inst )
    
    op_transcripts.add_transcripts( transcripts, gene, WRITE_META_DATA )
    
    return

def estimate_genes_expression_worker( genes_source_input_queue, output_queue, \
                                      transcripts_ofp,                        \
                                      fl_dists, read_group_mappings, \
                                      processed_genes = []):
    """Estimate transcript frequencies for genes.
    
    Load genes and source bam_fns from genes_source_input_queue, 
        process them, and write them to the output_queue.
    
    We need to initialise the bam files first
    """    
    def process_queue():
        # load the bam files
        reads_objs = build_reads_objs( bam_fns, fl_dists, read_group_mappings )
        
        while not genes_source_input_queue.empty():
            try:
                gene, candidate_transcripts, bam_fn = \
                    genes_source_input_queue.get(block=False)
                processed_genes.append( ( gene, bam_fn ) )
            except Queue.Empty:
                break

            # if the reads don't exist, we pbly had trobule 
            # with the fl dist ( or something )
            try:
                reads = reads_objs[ bam_fn ]
            except KeyError: 
                continue

            item = process_gene_and_reads( \
                gene, candidate_transcripts, reads, transcripts_ofp )

            # typically process_gene_and_reads will just return a transcripts object,
            # but it's possible for it to return a GeneProcessingError as well. 
            output_queue.put( (item, gene.name) )
    
    def handle_error( x, frame ):
        gene = processed_genes[-1][0]
        reads_fname = processed_genes[-1][1]
        output_queue.put( (make_error_log_string( \
                    gene, reads_fname, "Out of memory error."), gene.name) )
        sys.exit( 0 )
    
    signal.signal(signal.SIGUSR1, handle_error)
    
    process_queue()
        
    return

def process_output_queue( output_queue, log_fp, num_sources ):
    """Process the output queue.
    
    """
    num_genes_processed = 0
    while not output_queue.empty():
        try:
            item, gene_name = output_queue.get(block=False)
        except Queue.Empty:
            return num_genes_processed
        
        num_genes_processed += 1
        
        # if this is an exception, write it to the log imemdiately. 
        if isinstance( item, str ):
            #if isinstance( item, str ):
            log_fp.write( item + '\n' )
            log_fp.flush()
        else:
            assert type(None) == type(item)
        
    return num_genes_processed

def estimate_genes_expression( genes, gene_transcripts, bam_fns, \
                               fl_dists, read_group_mappings, \
                               num_threads, ofname ):
    """

    """
    ofp = TranscriptsFile( ofname, "w" )

    log_fname = ofname + ".log"
    log_fp = open( log_fname, 'w' )
    
    # create queues to store input and output data
    manager = multiprocessing.Manager()
    input_queue = manager.Queue()
    output_queue = manager.Queue()
    
    ## populate the input queue with each gene, candidate_transcripts, reads combination
    # find how many genes to get 100000 into the queue
    NUM = 100000
    min_num_genes = ( NUM/len( bam_fns ) ) + 1
    cnt = 0
    for gene in genes[:min_num_genes]:
        for bam_fn in bam_fns:
            if VERBOSE and cnt%10000 == 0:
                print "{0:d} / {1:d} ( {2:d} {3:d} )".format( \
                    cnt, len(genes)*len(bam_fns), len(genes), len(bam_fns) )
            input_queue.put( ( gene, gene_transcripts[gene.name], bam_fn ) )
            cnt += 1
    
    # start the gene processing
    args = ( input_queue, output_queue, ofp, fl_dists, read_group_mappings )
    
    if PROCESS_SEQUENTIALLY:
        for gene in genes[min_num_genes:]:
            for bam_fn in bam_fns:
                if VERBOSE and cnt%10000 == 0:
                    print "{0:d} / {1:d} ( {2:d} {3:d} )".format( \
                        cnt, len(genes)*len(bam_fns), len(genes), len(bam_fns) )
                input_queue.put( ( gene, gene_transcripts[gene.name], bam_fn ) )
                cnt += 1

        estimate_genes_expression_worker( *args )
    else:
        # spawn threads to estimate genes expression
        processes = []
        for thread_id in xrange( num_threads ):
            p = multiprocessing.Process(\
                target=estimate_genes_expression_worker, args=args)
            p.start()
            processes.append( p )

    for gene in genes[min_num_genes:]:
        for bam_fn in bam_fns:
            if VERBOSE and cnt%10000 == 0:
                print "{0:d} / {1:d} ( {2:d} {3:d} )".format( \
                    cnt, len(genes)*len(bam_fns), len(genes), len(bam_fns) )
            input_queue.put( ( gene, gene_transcripts[gene.name], bam_fn ) )
            cnt += 1
        
    ### process the output queue. 
    # Keep a dictionary of all of the output transcripts objects, grouped by gene.
    # we do this because after all of the read source files have been processed, we
    # need to create a merged gtf file.
    num_genes_processed = 0
    while (not PROCESS_SEQUENTIALLY) and any( p.is_alive() for p in processes ):
        num_new_genes_processed = process_output_queue( \
            output_queue, log_fp, len(bam_fns) )
        
        for index, process in enumerate( processes ):
            memory_used = get_memory_usage( process.pid )
            
            if memory_used > 100.0/len( processes ):
                print "WARNING: terminating process for using {0:.2%} memory".format( \
                    memory_used/100 )
                os.kill( process.pid, signal.SIGUSR1 )
                p = multiprocessing.Process(\
                    target=estimate_genes_expression_worker, args=args)
                p.start()
                processes[ index ] = p
        
        num_genes_processed += num_new_genes_processed
        if num_new_genes_processed > 0 :
            print '{0:.0%} of genes completed ({1:d}/{2:d})'.format( \
                num_genes_processed / float( len(bam_fns) * len(genes) ), \
                    num_genes_processed, ( len(bam_fns) * len(genes) ) )
        time.sleep( 0.1 )
    
    # get any remaining entries in the ouput queue
    num_genes_processed += process_output_queue( \
        output_queue, log_fp, len(bam_fns) )
    print "100%% of genes completed ({0:d}/{1:d})".format( \
        num_genes_processed, num_genes_processed )
    
    log_fp.close()
    
    return num_genes_processed

def get_raw_transcripts( gtf_fp ):
    """ Get gene transcript structure from gtf file
    Also verify that transcripts are valid.
    
    This function is used in get_elements scripts so it is not a contained function 
    of build_gene_transcripts
    """
    def is_valid_transcript( exons ):
        # check that there are no repeated exons
        if sorted( set( exons ) ) != sorted( exons ):
            return False
        # check that consecutive exons are not overlapping and are on the 
        # same strand and chrm
        for i, exon in enumerate( exons[:-1] ):
            if exon.stop >= exons[i+1].start or \
                    exon.strand != exons[i+1].strand or \
                    exon.chr != exons[i+1].chr:
                return False
        
        return True
    
    # struct of raw_gene_transcripts: gene_name -> trans_name -> exon_regions
    # fill raw_gene_transcripts with raw exon inforamtion grouped by gene and trans
    raw_gene_transcripts = defaultdict(lambda : defaultdict(set))
    for line in gtf_fp:
        gene_name, trans_name, feature_type, region = parse_gff_line( line )
        if feature_type == 'exon':
            raw_gene_transcripts[gene_name][trans_name].add( region )
    
    # verify transcripts
    verified_raw_gene_trans = defaultdict(dict)
    for gene_name, transcripts in raw_gene_transcripts.iteritems():
        for trans_name, exons in transcripts.iteritems():
            exons = sorted( exons )
            if not is_valid_transcript( exons ):
                if VERBOSE:
                    print 'Transcript ' + trans_name + ' contained invalid exons ' + \
                        'and was removed.'
                continue
            
            verified_raw_gene_trans[gene_name][trans_name] = exons
        
        if len( verified_raw_gene_trans[gene_name] ) == 0:
            if VERBOSE:
                print 'Gene ' + gene_name + \
                    ' contained no valid transcripts and will be removed.'
            del verified_raw_gene_trans[gene_name]
    
    return verified_raw_gene_trans

def build_gene_transcripts( gtf_fp, genes ):
    """ Build a transcripts object from each cooresponding entry in gtf_fp
    returned object has structure dict[gene_name]->Transcripts object
    """
    # struct of raw_gene_transcripts: gene_name -> trans_name -> exon_regions
    raw_gene_transcripts = get_raw_transcripts( gtf_fp )
    
    # initialize transcripts objects for each gene
    gene_transcripts = {}
    for gene_name, gene in genes.iteritems():
        gene_transcripts[ gene_name ] = Transcripts( gene )
    
    # create transcripts objects for each gene
    for gene_name, transcripts in raw_gene_transcripts.iteritems():
        for trans_name, exons in transcripts.iteritems():
            # get exon indices from gene object for each transcript
            exon_indices = []
            for exon in exons:
                exon_indices.append( \
                    genes[gene_name].exon_index( exon.start, exon.stop ) )

            # create transcript object and add it to its transcripts object
            transcript = Transcript( \
                [exon_indices[0]], [exon_indices[-1]], exon_indices[1:-1], trans_name )
            gene_transcripts[gene_name].add_transcript( transcript )
    
    return gene_transcripts

def build_objects( gtf_fp ):
    """Build objects which can be passed to processes simutaniously

    Read objects must be contained in each process separately b/c 
    they are linked to open files
    """
    # create gene object and close gtf_file
    genes = GeneBoundaries( gtf_fp )
    if MINIMAL_VERBOSE:
        print "Built gene objects from gtf file."
    # move file position back to beginning of file to be read for creating transcripts
    gtf_fp.seek(0)
    
    # create gene_transcripts dict
    gene_transcripts = build_gene_transcripts( gtf_fp, genes )
    gtf_fp.close()
    if MINIMAL_VERBOSE:
        print "Built transcripts objects from gtf file."
        
    return genes, gene_transcripts

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
        
    global VERBOSE
    VERBOSE = args.verbose
    
    if args.fl_dist_norm:
        try:
            mean, sd = args.fl_dist_norm.split(':')
            mean = int(mean)
            sd = int(sd)
            fl_dist_norm = (mean, sd)
        except ValueError:
            print "WARNING: User entered mean and sd for normal fl_dist " \
                + "are not properly formatted.\n" + \
                "\tUsing default produced from input BAM."
            fl_dist_norm = None
    else:
        fl_dist_norm = None

    # we change to the output directory later, and these files need to opened in
    # each sub-process for thread safety, so we get the absokute path while we can.
    bam_fns = [ os.path.abspath( bam_fn ) for bam_fn in args.bam_fns ]

    global PROCESS_SEQUENTIALLY
    if args.threads == 1:
        PROCESS_SEQUENTIALLY = True
    
    global WRITE_META_DATA
    WRITE_META_DATA = args.write_meta_data
        
    return args.gtf, bam_fns, args.ofname, \
        args.fl_dists, fl_dist_norm, args.threads

if __name__ == "__main__":
    # Get file objects from command line
    gtf_fp, bam_fns, ofname, fl_dist_fns, fl_dist_norm, threads=parse_arguments()
    
    # build objects from file objects
    genes, gene_transcripts = build_objects( gtf_fp )
    
    if fl_dist_norm:
        mean, sd = fl_dist_norm
        fl_min = max( 0, mean - (4 * sd) )
        fl_max = mean + (4 * sd)
        fl_dists = build_normal_density( fl_min, fl_max, mean, sd )
    else:
        fl_dists, read_group_mappings = load_fl_dists( fl_dist_fns )


    def foo():
        estimate_genes_expression( \
            genes.values(), gene_transcripts, bam_fns, \
            fl_dists, read_group_mappings, threads, ofname )
    
    if DO_PROFILE:
        import cProfile
        cProfile.run('foo()')
    else:
        foo()
        
