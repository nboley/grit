import os, sys
import pickle

from multiprocessing import BoundedSemaphore

VERBOSE = False
import psycopg2

from load_gene_from_db import load_gene_from_db

# for fl dist object
sys.path.append( os.path.join(os.path.dirname(__file__), "../sparsify/") )
import new_sparsify_transcripts
from new_sparsify_transcripts import build_design_matrices
from new_sparsify_transcripts import estimate_transcript_frequencies

def estimate_transcript_frequencies(conn, gene_id, bam_fn, fl_dists):
    """

    """
    try:
        cursor = conn.cursor()
        
        # load the genes
        gene = load_gene_from_db( gene_id, conn )
        gene.chrm = gene.chrm[0]

        # build the design matrices
        expected_array, observed_array, unobservable_transcripts \
            = new_sparsify_transcripts.build_design_matrices( 
                gene, bam_fn, fl_dists, None )

        # estimate the transcript frequencies
        ABS_TOL = 1e-4
        mle_estimates = new_sparsify_transcripts.estimate_transcript_frequencies( 
            observed_array, expected_array, ABS_TOL)
        if VERBOSE:
            print "Estimates:", mle_estimates
    
        # add the transcript frequencies into the DB
        for transcript, mle_estimate in zip( gene.transcripts, mle_estimates ):
            query_template  = "INSERT INTO transcript_expression "
            query_template += "( transcript, reads_fn, frequency ) VALUES "
            query_template += "( %s, %s, %s );"
            cursor.execute( query_template, 
                            (transcript.id, bam_fn, mle_estimate) )

        query = "UPDATE gene_expression_queue " \
              + "SET processing_status = 'FINISHED' " \
              + "WHERE gene = '%s'" % gene_id \
              + "  AND reads_fn = '%s';" % bam_fn

        cursor.execute( query )
    except Exception, inst:
        if VERBOSE: print "ERROR in %s:" % gene_id, inst
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

def load_fl_dists( bam_fn ):
    with open( bam_fn + ".fldist" ) as fl_dists_fp:
        fl_dists = pickle.load( fl_dists_fp )
    return fl_dists

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='Estimate transcript frequencies from DB.')

    # TODO add code to determine this automatically
    parser.add_argument( '--annotation-name', required=True, type=int,
                         help='Annotation for which to estimate transcript frequencies.')
    
    parser.add_argument( '--host', default='localhost', help='Database host.' )
    parser.add_argument( '--db-name', default='rnaseq_data', 
                         help='Database name. ' )

    parser.add_argument( '--threads', '-t', default=1,
                         help='Number of threads to spawn.')    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',
                         help='Print status information.')
    parser.add_argument( '--debug-verbose', default=False, action='store_true',
                         help='Print additional status information.')

    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose or args.debug_verbose
    new_sparsify_transcripts.VERBOSE = VERBOSE
    
    global DEBUG_VERBOSE
    DEBUG_VERBOSE = args.debug_verbose
    new_sparsify_transcripts.DEBUG_VERBOSE = DEBUG_VERBOSE
    
    return ( args.db_name, args.host), args.annotation_name, int(args.threads)

def gene_queue_is_empty(cursor):
    query = "SELECT * FROM gene_expression_queue" \
            + " WHERE processing_status = 'UNPROCESSED' LIMIT 1;"
    cursor.execute( query )
    return len( cursor.fetchall() ) == 0

def clear_zombies( pids ):
    for pid in sorted(pids):
        with open("/proc/%d/stat" %(pid)) as procfile:
            status = procfile.readline().split()[2]
        if status == 'Z': 
            os.waitpid( pid, 0 )
            pids.remove(pid)
    
    return

def main():
    conn_info, ann_name, nthreads  = parse_arguments()
    parent_conn = psycopg2.connect("dbname=%s host=%s" % conn_info)
    cursor = parent_conn.cursor()
    thread_cntr = BoundedSemaphore( nthreads )
    processes = set()
    while True:
        # wait4 zombied processes
        clear_zombies( processes )
        
        if gene_queue_is_empty(cursor): 
            break
        
        # acquire a thread, and for. We us ethe double fork
        # trick to avoid zombie processes
        thread_cntr.acquire()
        pid = os.fork()
        if pid != 0: 
            processes.add( pid )
            continue
        
        # open a new connection ( which we need to do because of the fork )
        conn = psycopg2.connect("dbname=%s host=%s" % conn_info)
        
        # get and lock a gene to process
        gene_id, reads_fn = get_queue_item(conn)
        fl_dist = load_fl_dists(reads_fn)
        
        if VERBOSE: print "Processing ", gene_id
        
        estimate_transcript_frequencies( 
            conn, gene_id, reads_fn, fl_dist  )

        conn.close()
        # release the thread, and stop the process
        thread_cntr.release()
        os._exit(os.EX_OK)
    
    # wait until all of the threads have terminated
    for loop in range(nthreads): thread_cntr.acquire()

if __name__ == '__main__':
    main()
