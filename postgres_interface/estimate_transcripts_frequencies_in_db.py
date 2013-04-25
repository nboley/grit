import os, sys
import pickle
import time
from itertools import izip

from multiprocessing import BoundedSemaphore, Process, cpu_count

VERBOSE = False
DEBUG = False

import psycopg2
import pysam

from load_gene_from_db import load_gene_from_db
from build_fl_dist_from_db import build_fl_dist
from db_conn import open_conn

import s3_data

# for fl dist object
sys.path.append( os.path.join(os.path.dirname(__file__), "../sparsify/") )
import new_sparsify_transcripts
from new_sparsify_transcripts import build_design_matrices
from new_sparsify_transcripts import estimate_transcript_frequencies
from new_sparsify_transcripts import calc_fpkm, TooFewReadsError

def add_freq_estimates_to_DB( cursor, gene, bam_fn, freq_estimates, fpkms ):
    for transcript, freq_estimate, fpkm in zip( 
            gene.transcripts, freq_estimates, fpkms ):
        query_template  = "INSERT INTO transcript_expression "
        query_template += "( transcript, reads_fn, frequency, rpkm ) VALUES "
        query_template += "( %s, %s, %s, %s );"
        cursor.execute( query_template, 
                        (transcript.id, bam_fn, freq_estimate, fpkm) )
    
    return

def add_design_matrices_to_DB( cursor, gene, bam_fn, 
                               expected_array, observed_array, 
                               unobservable_transcripts ):
    # insert the observed array
    query_template = """INSERT INTO observed_arrays 
                        ( gene, reads_fn, observed_bin_cnts ) 
                        VALUES ( %s, %s, %s );"""
    #cursor.execute( query_template, ( gene.id, bam_fn, observed_array ) )
    query = cursor.mogrify( query_template, 
                            (gene.id, bam_fn, observed_array.tolist()) )
    cursor.execute( query )
    
    for trans_i, (transcript, expected) in enumerate( 
            izip(gene.transcripts, expected_array.tolist()) ):
        if trans_i in unobservable_transcripts: continue
        query_template = """INSERT INTO expected_arrays 
                            ( transcript, reads_fn, expected_bin_fracs ) 
                            VALUES ( %s, %s, %s );"""
        query = cursor.mogrify( query_template, 
                                ( transcript.id, bam_fn, expected ) )
        cursor.execute( query )
    
    return

def estimate_transcript_frequencies(conn, gene_id, reads_key, bam_fn, fl_dists):
    try:
        cursor = conn.cursor()
        
        # load the genes
        gene = load_gene_from_db( gene_id, conn )
        gene.chrm = "chr" + gene.chrm[0]

        try:
            # build the design matrices
            expected_array, observed_array, unobservable_transcripts \
                = new_sparsify_transcripts.build_design_matrices( 
                    gene, bam_fn, fl_dists, None, reverse_strand=False )
            if DEBUG_VERBOSE: print "Sum counts:", observed_array.sum()
            add_design_matrices_to_DB( cursor, gene, reads_key, 
                                       expected_array, observed_array, 
                                       unobservable_transcripts )

            # estimate the transcript frequencies, and add them into the DB
            mle_estimates = new_sparsify_transcripts.estimate_transcript_frequencies( 
                observed_array, expected_array, 1e-5)
            bam_file = pysam.Samfile(bam_fn)
            fpkms = calc_fpkm( gene, fl_dists, mle_estimates, bam_file.mapped, 
                               observed_array.sum() )
            bam_file.close()
        except TooFewReadsError:
            mle_estimates = [None]*len(gene.transcripts)
            fpkms = [0.]*len(gene.transcripts)
        
        if DEBUG_VERBOSE:
            print "Estimates:", mle_estimates
        add_freq_estimates_to_DB(cursor, gene, reads_key, mle_estimates, fpkms)
        
        # updae the queue
        query = "UPDATE gene_expression_queue " \
              + "SET processing_status = 'FINISHED' " \
              + "WHERE gene = '%s'" % gene_id \
              + "  AND reads_fn = '%s';" % reads_key
        cursor.execute( query )
    except Exception, inst:
        if VERBOSE: print "ERROR in %s:" % gene_id, inst
        if DEBUG: raise
        inst = str(inst).replace( "'", "" )
        query = "UPDATE gene_expression_queue " \
              + "SET processing_status = 'FAILED', " \
              + "error_log = '%s' " % inst \
              + "WHERE gene = %s" % gene_id \
              + "  AND reads_fn = '%s';" % reads_key
        cursor.execute( query )
        
    cursor.close()
    conn.commit()
    
    return 

def get_queue_items(conn, max_number=1 ):
    """Get an item from the queue, return None if it's empty

    """
    cursor = conn.cursor()
    query = """
    UPDATE gene_expression_queue 
       SET processing_status = 'PROCESSING' 
     WHERE (gene, reads_fn) in (
            SELECT gene, reads_fn
            FROM gene_expression_queue 
            WHERE processing_status = 'UNPROCESSED' 
            FOR UPDATE
            LIMIT %i
    ) RETURNING gene, reads_fn;
    """ % max_number
    cursor.execute( query )
    res = cursor.fetchall()
    cursor.close()
    conn.commit()
        
    return res

def gene_queue_is_empty(cursor):
    query = "SELECT * FROM gene_expression_queue" \
            + " WHERE processing_status = 'UNPROCESSED' LIMIT 1;"
    cursor.execute( query )
    return len( cursor.fetchall() ) == 0

def load_fl_dists( bam_fn ):
    with open( bam_fn + ".fldist" ) as fl_dists_fp:
        fl_dists = pickle.load( fl_dists_fp )
    return fl_dists

def get_reads_fn_and_fl_dist(cache, reads_url):
    reads_fp = cache.open_local_copy_from_s3_url(reads_url)
    reads_fn = reads_fp.name
    reads_fp.close()
    fl_dist = load_fl_dists( reads_fn )
    return reads_fn, fl_dist

def spawn_process( conn_info, cache ):
    # open a new connection ( which we need to do because of the fork )
    conn = open_conn(*conn_info)
    
    # get and lock a gene to process
    rv = get_queue_items(conn)
    if len(rv) == 0: return
    
    for gene_id, reads_location in rv:
        try:
            if VERBOSE: print "Processing ", gene_id
          
            if reads_location.startswith( "s3" ):
                reads_fn, fl_dist \
                    = get_reads_fn_and_fl_dist( cache, reads_location )
            else:
                reads_fn, fl_dist = \
                    reads_location, load_fl_dists( reads_location )
            
            estimate_transcript_frequencies( 
                conn, gene_id, reads_location, reads_fn, fl_dist  )
        except ValueError, inst:
            inst = str(inst).replace( "'", "" )
            cursor = conn.cursor()
            query = "UPDATE gene_expression_queue " \
                  + "SET processing_status = 'FAILED', " \
                  + "error_log = '%s' " % inst \
                  + "WHERE gene = %s" % gene_id \
                  + "  AND reads_fn = '%s';" % reads_location
            cursor.execute( query )
            cursor.close()
            conn.commit()
            if DEBUG_VERBOSE: print "CLEANED UP ", inst, reads_location
            continue


    conn.close()
    return

def main_loop(conn_info, s3_info, nthreads, is_daemon):
    parent_conn = open_conn(*conn_info)
    cursor = parent_conn.cursor()
    
    cache = s3_data.S3Cache( s3_info[0], s3_info[1], parent_conn ) \
        if s3_info != None else None
    processes = []
    while True:
        running_ps = []
        for p in processes:
            if p.is_alive():
                running_ps.append( p )
        processes = running_ps
        if len( processes ) == nthreads:
            time.sleep( 1.0 )
            continue
        
        if gene_queue_is_empty(cursor): 
            if is_daemon:
                time.sleep(1)
                continue
            else:
                break                
        
        # release the thread, and stop the process
        p = Process( target=spawn_process, args=[conn_info, cache] )
        p.start()
        processes.append( p )
    
    # wait until all of the threads have terminated
    for p in processes:
        p.join()

    return

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='Estimate transcript frequencies from DB.')

    parser.add_argument( '--db-name', default='rnaseq', 
                         help='Database to insert the data into. ' )
    parser.add_argument( '--db-host', 
                         help='Database host. default: socket connection' )
    parser.add_argument( '--db-user', 
                         help='Database connection user. Default: unix user' )
    parser.add_argument( '--db-pass', help='DB connection password.' )

    parser.add_argument( '--access-key', default=None, help='Amazon access key.' )
    parser.add_argument( '--secret-key', default=None, help='Amazon secret key.' )

    parser.add_argument( '--threads', '-t', default='MAX',
                         help='Number of threads to spawn ( MAX for the number of cores ).')    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',
                         help='Print status information.')
    parser.add_argument( '--debug-verbose', default=False, action='store_true',
                         help='Print additional status information.')
    parser.add_argument( '--daemon', default=False, action='store_true',
                         help='Whether or not to run this as a daemon.')
    
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose or args.debug_verbose
    new_sparsify_transcripts
    s3_data.VERBOSE = VERBOSE
    
    global DEBUG_VERBOSE
    DEBUG_VERBOSE = args.debug_verbose
    new_sparsify_transcripts.DEBUG_VERBOSE = DEBUG_VERBOSE
    s3_data.DEBUG_VERBOSE = DEBUG_VERBOSE
    
    if args.threads == 'MAX':
        args.threads = cpu_count()
    
    s3_info = None if args.access_key == None or args.secret_key == None \
        else (args.access_key, args.secret_key )
    return ( args.db_name, args.db_host, args.db_user, args.db_pass), \
        int(args.threads), args.daemon, s3_info

def main():
    conn_info, nthreads, is_daemon, s3_info  = parse_arguments()
    if is_daemon:
        import daemon
        with daemon.DaemonContext():
            main_loop(conn_info, s3_info, nthreads, is_daemon)
    else:
        main_loop(conn_info, s3_info, nthreads, is_daemon)

if __name__ == '__main__':
    main()
