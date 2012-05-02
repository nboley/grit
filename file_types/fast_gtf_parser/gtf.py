import sys, os
import numpy
import signal
from ctypes import *

VERBOSE = False

# load the gtf library
gtf_o = cdll.LoadLibrary( os.path.join( ".", \
        os.path.dirname( __file__), "libgtf.so" ) )

class c_transcript_t(Structure):
    """
struct transcript {
    char* trans_id;
    int num_exon_bnds;
    int* exon_bnds;
};
"""
    _fields_ = [("trans_id", c_char_p),

                ("num_exon_bnds", c_int),
                ("exon_bnds", POINTER(c_int))
               ]

class c_gene_t(Structure):
    """
struct gene {
    char* gene_id;
    char* chrm;
    char strand;
    int min_loc;
    int max_loc;
    int num_transcripts;
    struct transcript** transcripts;
};
"""
    _fields_ = [("gene_id", c_char_p),

                ("chrm", c_char_p),
                ("strand", c_char),

                ("min_loc", c_int),
                ("max_loc", c_int),
                
                ("num_transcripts", c_int),
                ("transcripts", POINTER(POINTER(c_transcript_t)))
               ]

def load_gtf( fname ):
    c_genes_p = c_void_p()   
    num_genes = c_int()
    gtf_o.get_gtf_data( c_char_p(fname), byref(c_genes_p), byref( num_genes) )
    # check for a NULL pointer
    if not bool(c_genes_p):
        return []
        # raise IOError, "Couldn't load '%s'" % fname
    
    c_genes_p = cast( c_genes_p, POINTER(POINTER(c_gene_t)) )

    genes = []
    for i in xrange( num_genes.value ):
        g = c_genes_p[i].contents
        # convert the meta data into standard python types
        chrm = g.chrm
        if chrm.startswith( "chr" ):
            chrm = chrm[3:]
        data = [g.gene_id, chrm, g.strand, g.min_loc, g.max_loc]
        transcripts = []
        for j in xrange( g.num_transcripts ):
            trans_id = g.transcripts[j].contents.trans_id
            num_exon_bnds = g.transcripts[j].contents.num_exon_bnds
            exon_bnds = g.transcripts[j].contents.exon_bnds[:num_exon_bnds]
            transcripts.append( (trans_id, exon_bnds) )
        
        data.append( transcripts )
        genes.append( data )
    
    gtf_o.free_gtf_data( c_genes_p, num_genes )    
    
    return genes

def load_gtfs( fnames, num_threads ):
    """Multi threaded gtf parser.
    """
    from multiprocessing import Process, Queue, Manager
    # build a thread safe queue of filenames. WE add items to
    # the queue in the order of smallest file to largest,
    # because items are popped off of the top and we want the
    # largest files to be processed first. 
    input_q = Queue()
    fnames_and_sizes = [ (os.path.getsize(fname), fname ) \
                            for fname in fnames ]
    fnames_and_sizes.sort(reverse = True)
    for size, fname in fnames_and_sizes:
        input_q.put( fname )
    
    # initialize a data structure to hold the parsed gtf objects
    manager = Manager()
    output = manager.dict()

    # a wrapper that reads fromt he input queue, and
    # write to the output
    def foo( input_q, output ):
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        while not input_q.empty():
            try:
                fname = input_q.get()
            except Queue.Empty, inst:
                break
            
            output[fname] = load_gtf( fname )

            if VERBOSE:
                print >> sys.stderr, "Finished parsing ", fname
            
        return
    
    # initialize the list of processes
    procs = []
    for p_index in xrange( min( len(fnames), num_threads ) ):
        p = Process( target=foo, args=(input_q, output) )
        procs.append( p )
        p.start()
    
    # join all of the processes
    try:
        for p in procs:
            p.join()
    except KeyboardInterrupt:
        for p in procs:
            print "Terminating " + str(p) + "due to SIGINT"
            p.terminate()
    
    # return the data in the same order we got it
    rv = []
    for fname in fnames:
        rv.append( output[ fname ]  )
    
    return rv

if __name__ == "__main__":
    VERBOSE = True
    gene_grps = load_gtfs( sys.argv[1:], 4 )
    for genes in gene_grps:
        for gene in genes:
            print gene[0]
            #print "\t", [ i[0] for i in gene[-1] ]
