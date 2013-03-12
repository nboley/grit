import os, sys
import pickle
import time
from itertools import izip

from multiprocessing import BoundedSemaphore, Process

VERBOSE = False
import psycopg2

from load_gene_from_db import load_gene_from_db
from build_fl_dist_from_db import build_fl_dist

# for fl dist object
sys.path.append( os.path.join(os.path.dirname(__file__), "../sparsify/") )
import new_sparsify_transcripts
from new_sparsify_transcripts import build_design_matrices
from new_sparsify_transcripts import estimate_transcript_frequencies

def add_freq_estimates_to_DB( cursor, gene, bam_fn, freq_estimates ):
    for transcript, freq_estimate in zip( gene.transcripts, freq_estimates ):
        query_template  = "INSERT INTO transcript_expression "
        query_template += "( transcript, reads_fn, frequency ) VALUES "
        query_template += "( %s, %s, %s );"
        cursor.execute( query_template, 
                        (transcript.id, bam_fn, freq_estimate) )
    
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

def estimate_transcript_frequencies(conn, gene_id, bam_fn, fl_dists):
    try:
        cursor = conn.cursor()
        
        # load the genes
        gene = load_gene_from_db( gene_id, conn )
        gene.chrm = gene.chrm[0]

        # build the design matrices
        expected_array, observed_array, unobservable_transcripts \
            = new_sparsify_transcripts.build_design_matrices( 
                gene, bam_fn, fl_dists, None )

        add_design_matrices_to_DB( cursor, gene, bam_fn, 
                                   expected_array, observed_array, 
                                   unobservable_transcripts )
        
        # estimate the transcript frequencies, and add them into the DB
        ABS_TOL = 1e-4
        mle_estimates = new_sparsify_transcripts.estimate_transcript_frequencies( 
            observed_array, expected_array, ABS_TOL)
        if VERBOSE:
            print "Estimates:", mle_estimates
        add_freq_estimates_to_DB( cursor, gene, bam_fn, mle_estimates )
        
        # updae the queue
        query = "UPDATE gene_expression_queue " \
              + "SET processing_status = 'FINISHED' " \
              + "WHERE gene = '%s'" % gene_id \
              + "  AND reads_fn = '%s';" % bam_fn

        cursor.execute( query )
    except Exception, inst:
        if VERBOSE: print "ERROR in %s:" % gene_id, inst
        raise
        inst = str(inst).replace( "'", "" )
        query = "UPDATE gene_expression_queue " \
              + "SET processing_status = 'FAILED', " \
              + "error_log = '%s' " % inst \
              + "WHERE gene = %s" % gene_id \
              + "  AND reads_fn = '%s';" % bam_fn
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
            ORDER BY reads_fn
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

    parser.add_argument( '--threads', '-t', default=1,
                         help='Number of threads to spawn.')    
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
    
    return ( args.db_name, args.host), int(args.threads), args.daemon

def get_reads_fn_and_fl_dist( manager_inst, reads_location, conn ):
    # check to see if the filename is already in the manager
    manager_inst.lock.acquire()
    # if it is, then find the corresponding filenames and return them
    if reads_location in manager_inst.location_to_fn_mapping:
        while True:
            reads_fn, frag_len_fn, status \
                = manager_inst.location_to_fn_mapping[reads_location]
            # if they are finished ( ie downloaded, and the fl dist
            # is built ) then return them
            if status == 'FINISHED':
                manager_inst.lock.release()
                return reads_fn, load_fl_dists(reads_fn)
            # otherwise, wait a second and try again
            else:
                if DEBUG_VERBOSE: 
                    print "Waiting on %s: status %s" % ( reads_fn, status )
                manager_inst.lock.release()
                time.sleep(1.)
                manager_inst.lock.acquire()
    # otherwise, we need to get them
    else:
        # choose a filename
        op_fname = "/scratch/" + ''.join( 
            x for x in reads_location.replace("/", "_") 
            if x.isalnum() or x == '_' ) + ".bam"
        manager_inst.location_to_fn_mapping[reads_location] = (
            op_fname, op_fname + ".fldist", 'downloading' )
        manager_inst.lock.release()
        
        # download the file
        if VERBOSE: print "Downloading ", reads_location
        res = os.system( "s3cmd get %s %s --continue" 
                         % (reads_location, op_fname) )
        res = os.system( "s3cmd get %s.bai %s.bai --continue" 
                         % (reads_location, op_fname) )
        
        # build the fl dist
        if not os.path.exists( op_fname + ".fldist" ):
            build_fl_dist(op_fname, conn)
        
        # release the lock
        manager_inst.lock.acquire()
        manager_inst.location_to_fn_mapping[reads_location] = (
            op_fname, op_fname + ".fldist", 'FINISHED' )        
        manager_inst.lock.release()

        return op_fname, load_fl_dists(op_fname)

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

    reads_fn, fl_dist \
        = get_reads_fn_and_fl_dist( manager_inst, reads_location, conn )

    if VERBOSE: print "Processing ", gene_id

    estimate_transcript_frequencies( 
        conn, gene_id, reads_fn, fl_dist  )

    conn.close()
    return

def main():
    conn_info, nthreads, daemon  = parse_arguments()
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
            if not daemon:
                break
            else:
                time.sleep(1)
        
        # release the thread, and stop the process
        p = Process( target=spawn_process, args=[conn_info,manager_inst] )
        p.start()
        processes.append( p )
    
    # wait until all of the threads have terminated
    for loop in range(nthreads): thread_cntr.acquire()

if __name__ == '__main__':
    main()
