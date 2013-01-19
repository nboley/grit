# Copyright (c) 2011-2012 Nathan Boley

import sys, os
import numpy
import signal
from ctypes import *
from itertools import izip, repeat, chain

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      ".." ) )
from gtf_file import create_gtf_line, GenomicInterval

VERBOSE = False

# load the gtf library
gtf_o = cdll.LoadLibrary( os.path.join( ".", \
        os.path.dirname( __file__), "libgtf.so" ) )

class c_transcript_t(Structure):
    """
struct transcript {
    char* trans_id;
    int cds_start;
    int cds_stop;
    int num_exon_bnds;
    int* exon_bnds;
    int score;
    int rpkm;
    int rpk;
};
"""
    _fields_ = [("trans_id", c_char_p),
                
                ("cds_start", c_int),
                ("cds_stop", c_int),
                
                ("num_exon_bnds", c_int),
                ("exon_bnds", POINTER(c_int)),

                ("score", c_int),
                ("rpkm", c_double),
                ("rpk", c_double)
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

class Transcript( list ):
    def __init__(self, trans_id, chrm, strand, exons, cds_region,
                 gene_id=None, score=None, rpkm=None, rpk=None ):
        self.gene_id = gene_id
        self.id = trans_id
        self.chrm = chrm
        self.strand = strand

        self.score = score
        self.rpkm = rpkm
        self.rpk = rpk
        
        exon_bnds = list( chain( *exons ) )
        self.exon_bnds = exon_bnds
        self.exons = tuple(zip(exon_bnds[:-1:2], exon_bnds[1::2]))
        assert list(self.exons) == list(exons)
        self.introns = tuple([ (x+1, y-1) for x, y in 
                               izip(exon_bnds[1:-2:2], exon_bnds[2:-1:2]) ])
        
        self.is_protein_coding = ( cds_region != None )
        self.cds_region = cds_region
        self.cds_exons = None
        self.fp_utr_exons = None
        self.tp_utr_exons = None

        if cds_region != None:
            # find the start index of the coding sequence
            cds_start_bndry_index = exon_bnds.index( cds_region[0] )
            cds_stop_bndry_index = exon_bnds.index( cds_region[1] )
            
            cds_exon_bnds = exon_bnds[cds_start_bndry_index:cds_stop_bndry_index+1]
            fp_bnds = exon_bnds[:cds_start_bndry_index]
            tp_bnds = exon_bnds[cds_stop_bndry_index+1:]

            self.cds_exons = tuple(zip(cds_exon_bnds[:-1:2], cds_exon_bnds[1::2]))
            self.us_exons = tuple(zip(fp_bnds[:-1:2], fp_bnds[1::2]))
            self.ds_exons = tuple(zip(tp_bnds[:-1:2], tp_bnds[1::2]))

            # if this is a reverse strand transcript, rev 5' and 3' ends
            if self.strand == '+':
                self.fp_utr_exons, self.tp_utr_exons \
                    = self.us_exons, self.ds_exons
            else:
                self.fp_utr_exons, self.tp_utr_exons \
                    = self.ds_exons, self.us_exons
            
            # now, remove CDS boundaries if they aren't real ( ie, if they just
            # differentiate between coding and non coding sequence )
            
            # start with the stop so we don't have to re-search for the index
            if cds_stop_bndry_index+1 < len(exon_bnds) \
                    and exon_bnds[ cds_stop_bndry_index+1 ] == cds_region[1]+1:
                del exon_bnds[ cds_stop_bndry_index:cds_stop_bndry_index+2 ]
            
            if cds_start_bndry_index > 0 \
                    and exon_bnds[ cds_start_bndry_index-1]+1 == cds_region[0]:
                del exon_bnds[ cds_start_bndry_index-1:cds_start_bndry_index+1 ]
                
        # add these for compatability
        self.append( trans_id )
        self.append( exon_bnds )        
    
    def __hash__( self ):
        if self.cds_region != None:
            return hash(( self.chrm, self.strand, 
                          self.exons, tuple(self.cds_region) ))
        else:
            return hash( (self.chrm, self.strand, self.exons, None) )
    
    def IB_key( self ):
        """Return a key for matching transcripts on their external bounds.
        
        """
        if len( self.exon_bnds ) == 2:
            return ( self.chrm, self.strand, "SE_GENE", self.cds_region )
        else:
            return (self.chrm, self.strand, 
                    tuple(self.exon_bnds[1:-1]), self.cds_region)
    
    def build_gtf_lines( self, gene_id, meta_data, source='.'):
        ret_lines = []
        def build_lines_for_feature( exons, feature, is_CDS=False ):
            current_frame = 0
            score = str(self.score) if self.score != None else '.'
            for start, stop in exons:
                region = GenomicInterval( self.chrm, self.strand, start, stop )
                frame = current_frame if is_CDS else '.'
                yield create_gtf_line( region, gene_id, self.id, meta_data,
                                       score, feature=feature, frame=str(frame),
                                       source=source)
                current_frame = ( current_frame + stop - start + 1 )%3
            return
        
        if self.cds_region == None:
            ret_lines.extend( build_lines_for_feature( 
                    self.exons, 'exon', False ) )
        else:
            us_exons, ds_exons = self.fp_utr_exons, self.tp_utr_exons
            us_label, ds_label = 'five_prime_UTR', 'three_prime_UTR'
            us_label, ds_label = 'exon', 'exon'
            if self.strand == '-': 
                us_exons, ds_exons = ds_exons, us_exons
                us_label, ds_label = ds_label, us_label
            
            ret_lines.extend( build_lines_for_feature( 
                    us_exons, us_label, False ) )

            ret_lines.extend( build_lines_for_feature( 
                    self.cds_exons, 'CDS', True ) )
            
            ret_lines.extend( build_lines_for_feature( 
                    ds_exons, ds_label, False ) )
        
        return "\n".join( ret_lines )
    
class Gene( list ):
    def __init__( self, id, chrm, strand, start, stop, transcripts ):
        self.id = id
        self.chrm = chrm
        self.strand = strand
        self.start = start
        self.stop = stop
        self.transcripts = transcripts
        list.extend( self, ( id, chrm, strand, start, stop, transcripts ) )
        return    

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
        transcripts = []
        for j in xrange( g.num_transcripts ):
            trans_id = g.transcripts[j].contents.trans_id
            cds_start = g.transcripts[j].contents.cds_start
            cds_stop = g.transcripts[j].contents.cds_stop
            score = None if g.transcripts[j].contents.score == -1 \
                else int(g.transcripts[j].contents.score)
            rpkm = None if g.transcripts[j].contents.rpkm < -1e-10 \
                else float(g.transcripts[j].contents.rpkm)
            rpk = None if g.transcripts[j].contents.rpk < -1e-10 \
                else float(g.transcripts[j].contents.rpk)
            
            num_exon_bnds = g.transcripts[j].contents.num_exon_bnds
            exon_bnds = g.transcripts[j].contents.exon_bnds[:num_exon_bnds]
            if cds_start == -1 or cds_stop == -1:
                assert cds_start == -1 and cds_stop == -1
                cds_region = None
            else:
                cds_region = ( cds_start, cds_stop )
            exons = tuple(zip(exon_bnds[:-1:2], exon_bnds[1::2]))
            transcripts.append( 
                Transcript(trans_id, chrm, g.strand, 
                           exons, cds_region, g.gene_id, 
                           score=score, rpkm=rpkm, rpk=rpk ) 
                )
        
        gene =  Gene( g.gene_id, chrm, g.strand, 
                      g.min_loc, g.max_loc, transcripts )
        genes.append( gene )
    
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

    # a wrapper that reads from the input queue, and
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
    gene_grps = load_gtfs( sys.argv[1:], 1 )
    for genes in gene_grps:
        for gene in genes:
            print gene[0]
            for trans in gene[-1]:
                print trans.chrm, trans.strand, trans.id, trans.cds_region
                print trans.fp_utr_exons
                print trans.cds_exons
                print trans.tp_utr_exons
                print trans.score
                print trans.rpkm
                print trans.rpk
                print
                print trans.exons
                print trans.introns
                print
