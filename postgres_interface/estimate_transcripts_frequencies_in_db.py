import os, sys
import pickle
import time
from itertools import izip

from multiprocessing import BoundedSemaphore, Process, cpu_count

VERBOSE = False
DEBUG = False

import psycopg2
import pysam
import boto
import boto.s3.bucket
import boto.s3.key

from load_gene_from_db import load_gene_from_db
from build_fl_dist_from_db import build_fl_dist

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
                    gene, bam_fn, fl_dists, None, reverse_strand=True )
            if DEBUG_VERBOSE: print "Sum counts:", observed_array.sum()
            #add_design_matrices_to_DB( cursor, gene, reads_key, 
            #                           expected_array, observed_array, 
            #                           unobservable_transcripts )

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

def get_queue_item(conn):
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
            LIMIT 1
    ) RETURNING gene, reads_fn;
    """
    cursor.execute( query )
    res = cursor.fetchall()
    cursor.close()
    conn.commit()
    
    # if we can't get an item, then just exit
    if len( res ) == 0: os._exit(os.EX_OK)
    return res[0]

def gene_queue_is_empty(cursor):
    query = "SELECT * FROM gene_expression_queue" \
            + " WHERE processing_status = 'UNPROCESSED' LIMIT 1;"
    cursor.execute( query )
    return len( cursor.fetchall() ) == 0

def load_fl_dists( bam_fn ):
    with open( bam_fn + ".fldist" ) as fl_dists_fp:
        fl_dists = pickle.load( fl_dists_fp )
    return fl_dists

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='Estimate transcript frequencies from DB.')

    parser.add_argument( '--host', default='localhost', help='Database host.' )
    parser.add_argument( '--db-name', default='rnaseq_data', 
                         help='Database name. ' )

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
    new_sparsify_transcripts.VERBOSE = VERBOSE
    
    global DEBUG_VERBOSE
    DEBUG_VERBOSE = args.debug_verbose
    new_sparsify_transcripts.DEBUG_VERBOSE = DEBUG_VERBOSE
    
    if args.threads == 'MAX':
        args.threads = cpu_count()
    
    return ( args.db_name, args.host), int(args.threads), args.daemon

def get_reads_fn_and_fl_dist( manager_inst, reads_location, conn ):
    # check to see if the filename is already in the manager
    # if it is, then find the corresponding filenames and return them
    if reads_location in manager_inst.location_to_fn_mapping:
        while True:
            reads_fn, frag_len_fn, status \
                = manager_inst.location_to_fn_mapping[reads_location]
            # if they are finished ( ie downloaded, and the fl dist
            # is built ) then return them
            if status == 'FINISHED':
                return reads_fn, load_fl_dists(reads_fn)
            elif status == 'ERROR':
                raise ValueError, "Cant get files"
            # otherwise, wait a second and try again
            else:
                if DEBUG_VERBOSE: 
                    print "Waiting on %s: status %s" % ( reads_fn, status )
                time.sleep(1.)
    # otherwise, we need to get them
    else:
        try:
            # choose a filename
            op_fname = "/scratch/" + ''.join( 
                x for x in reads_location.replace("/", "_") 
                if x.isalnum() or x == '_' ) + ".bam"
            manager_inst.location_to_fn_mapping[reads_location] = (
                op_fname, op_fname + ".fldist", 'downloading' )

            # download the file
            if VERBOSE: print "Downloading ", reads_location, " to ", op_fname
            
            # get the bucket
            conn = boto.connect_s3(
                "", 
                "")
            bucket = boto.s3.bucket.Bucket( conn, reads_location[5:].split("/")[0] )
            key = boto.s3.key.Key(bucket)
            # get the bam file
            key.key = "/".join(reads_location[5:].split("/")[1:])
            key.get_contents_to_filename(op_fname)
            
            # get the bai file
            key.key = "/".join(reads_location[5:].split("/")[1:]) + '.bai'
            key.get_contents_to_filename(op_fname)
            
            # build the fl dist
            if not os.path.exists( op_fname + ".fldist" ):
                build_fl_dist(op_fname, conn)
            
            fl_dists = load_fl_dists(op_fname)

            # release the lock
            manager_inst.location_to_fn_mapping[reads_location] = (
                op_fname, op_fname + ".fldist", 'FINISHED' )        
        except Exception, inst:
            if VERBOSE: print "ERROR loading files. %s" % inst
            # release the lock
            manager_inst.location_to_fn_mapping[reads_location] = (
                op_fname, op_fname + ".fldist", 'ERROR' )        
            
            raise ValueError, "Cant get files" + str(inst)
        
        return op_fname, fl_dists

    assert False

class FilesManager(object):
    def __init__(self):
        import multiprocessing
        self.manager = multiprocessing.Manager()
        self.lock = multiprocessing.Lock()
        self.location_to_fn_mapping = self.manager.dict()
        return

def spawn_process( conn_info, manager_inst ):
    # open a new connection ( which we need to do because of the fork )
    conn = psycopg2.connect("dbname=%s host=%s user=nboley" % conn_info)

    # get and lock a gene to process
    gene_id, reads_location = get_queue_item(conn)

    try:
        reads_fn, fl_dist \
            = get_reads_fn_and_fl_dist( manager_inst, reads_location, conn )
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
        conn.close()
        if DEBUG_VERBOSE: print "CLEANED UP ", inst, reads_location
        return

    if VERBOSE: print "Processing ", gene_id

    estimate_transcript_frequencies( 
        conn, gene_id, reads_location, reads_fn, fl_dist  )

    conn.close()
    return

def main_loop(conn_info, nthreads, is_daemon):
    parent_conn = psycopg2.connect("dbname=%s host=%s user=nboley" % conn_info)
    cursor = parent_conn.cursor()
    
    manager_inst = FilesManager()
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
            if not is_daemon:
                break
            else:
                time.sleep(30)
        
        # release the thread, and stop the process
        p = Process( target=spawn_process, args=[conn_info,manager_inst] )
        p.start()
        processes.append( p )
    
    # wait until all of the threads have terminated
    for p in processes:
        p.join()

    return

def main():
    conn_info, nthreads, is_daemon  = parse_arguments()
    if is_daemon:
        import daemon
        with daemon.DaemonContext():
            main_loop(conn_info, nthreads, is_daemon)
    else:
        main_loop(conn_info, nthreads, is_daemon)

if __name__ == '__main__':
    main()
