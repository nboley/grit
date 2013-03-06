import os, sys
import pickle

VERBOSE = False
import psycopg2

from load_gene_from_db import load_gene_from_db

# for fl dist object
sys.path.append( os.path.join(os.path.dirname(__file__), "../sparsify/") )
import new_sparsify_transcripts
from new_sparsify_transcripts import build_design_matrices
from new_sparsify_transcripts import estimate_transcript_frequencies

def estimate_transcript_frequencies(conn, ann_id, gene_name, bam_fn, fl_dists):
    """

    """
    try:
        cursor = conn.cursor()
        
        # load the genes
        gene = load_gene_from_db( (gene_name, ann_id), conn )
        gene.id = gene.id[0]
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
              + "WHERE gene = '%s'" % gene_name \
              + "  AND annotation = '%s'" % ann_id \
              + "  AND reads_fn = '%s';" % bam_fn

        cursor.execute( query )
    except Exception, inst:
        print "ERROR in %s:" % gene_name, inst
        inst = str(inst)
        inst.replace( "'", "" )
        query = "UPDATE gene_expression_queue " \
              + "SET processing_status = 'FAILED', " \
              + "error_log = '%s' " % inst \
              + "WHERE gene = '%s'" % gene_name \
              + "  AND annotation = '%s'" % ann_id \
              + "  AND reads_fn = '%s';" % bam_fn
        
        cursor.execute( query )
        
    
    cursor.close()
    conn.commit()
    
    return 

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
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
                             help='Print status information.')
    parser.add_argument( '--debug-verbose', default=False, action='store_true',\
                             help='Print additional status information.')

    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    new_sparsify_transcripts.VERBOSE = VERBOSE
    
    global DEBUG_VERBOSE
    DEBUG_VERBOSE = args.debug_verbose
    new_sparsify_transcripts.DEBUG_VERBOSE = DEBUG_VERBOSE
        
    conn = psycopg2.connect("dbname=%s host=%s" % ( args.db_name, args.host) )
    return conn, args.annotation_name

def load_fl_dists( bam_fn ):
    with open( bam_fn + ".fldist" ) as fl_dists_fp:
        fl_dists = pickle.load( fl_dists_fp )
    return fl_dists
    

def main():
    all_fl_dists = {}
    
    conn, ann_name  = parse_arguments()
    while True:
        cursor = conn.cursor()
        query = "SELECT gene, annotation, reads_fn  " \
                + "FROM gene_expression_queue " \
                + "WHERE processing_status = 'UNPROCESSED' LIMIT 1 FOR UPDATE;"
        cursor.execute( query )
        
        # get and parse the arguments
        res = cursor.fetchall()
        if len( res ) == 0: break
        gene_id, annotation, reads_fn = res[0]
        if reads_fn not in all_fl_dists:
            all_fl_dists[reads_fn] = load_fl_dists(reads_fn)
            
        if VERBOSE:
            print "Processing ", gene_id
        
        # update the queue so that we are processing it
        query = "UPDATE gene_expression_queue " \
              + "SET processing_status = 'PROCESSING' " \
              + "WHERE gene = '%s'" % gene_id \
              + "  AND annotation = '%s'" % annotation \
              + "  AND reads_fn = '%s';" % reads_fn
        
        cursor.execute( query )
        cursor.close()

        # commit so that the FOR UPDATE unlocks.
        conn.commit()
        
        estimate_transcript_frequencies( 
            conn, annotation, gene_id, reads_fn, all_fl_dists[ reads_fn ] )
    
    conn.close()

if __name__ == '__main__':
    main()
