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

class c_exon_t(Structure):
    """
struct exon {
    int start;
    int stop;
};
"""
    _fields_ = [
                ("start", c_int),
                ("stop", c_int),
               ]

class c_transcript_t(Structure):
    """
struct transcript {
    char* trans_id;
    int cds_start;
    int cds_stop;
    int num_exons;
    struct* exons;
    int score;
    int rpkm;
    int rpk;
};
"""
    _fields_ = [("trans_id", c_char_p),
                
                ("cds_start", c_int),
                ("cds_stop", c_int),
                
                ("num_exons", c_int),
                ("exons", POINTER(c_exon_t)),

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

def partition_coding_and_utr_segments( exons, cds_start, cds_stop ):
    """Split the exons into UTR and CDS exons.

    """
    # find the exon index the the start codon intersects
    cds_start_i = [ i for i, (start, stop) in enumerate(exons) 
                    if cds_start >= start and cds_start <= stop ]
    assert len( cds_start_i ) == 1
    cds_start_i = cds_start_i[0]
    
    # we start at the cds_start exon because we know this index must be >= 
    assert cds_stop >= cds_start
    cds_stop_i = [ i for i, (start, stop) in enumerate(exons[cds_start_i:]) 
                   if cds_stop >= start and cds_stop <= stop ]
    assert len( cds_stop_i ) == 1
    cds_stop_i = cds_stop_i[0] + cds_start_i
    
    def mod_external_bndrys( exons, lower_bnd, upper_bnd ):
        """If necessary, shrink the external boundaries"""
        # if it's an empty set, there is nothing to be done
        if len( exons ) == 0: return exons
        exons[0] = ( max(exons[0][0], lower_bnd), exons[0][1] )
        exons[-1] = ( exons[-1][0], min(upper_bnd, exons[-1][1] ) )
        return exons
    
    cds_exons = mod_external_bndrys( 
        list(exons[cds_start_i:cds_stop_i+1]), cds_start, cds_stop )
    
    us_utr_stop_i = cds_start_i if exons[cds_start_i][0] < cds_start \
        else cds_start_i - 1
    us_utr_exons = mod_external_bndrys(
        list(exons[:us_utr_stop_i+1]), 1, cds_start-1)
    
    ds_utr_stop_i = cds_stop_i if exons[cds_stop_i][1] > cds_stop \
        else cds_stop_i + 1
    ds_utr_exons = mod_external_bndrys(
        list(exons[ds_utr_stop_i:]), cds_stop+1, 1e100)
    
    return us_utr_exons, cds_exons, ds_utr_exons

class Transcript( object ):
    def __init__(self, trans_id, chrm, strand, exons, cds_region,
                 gene_id=None, score=None, rpkm=None, rpk=None, promoter=None ):
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
        self.introns = tuple([ (x+1, y-1) for x, y in 
                               izip(exon_bnds[1:-2:2], exon_bnds[2:-1:2]) ])
        
        self.is_protein_coding = ( cds_region != None )
        
        self.cds_region = cds_region
        self.start_codon = None
        self.stop_codon = None
        
        self.cds_exons = None
        self.fp_utr_exons = None
        self.tp_utr_exons = None
        
        self.promoter = promoter
        
        if cds_region != None:
            self.us_exons, self.cds_exons, self.ds_exons = \
                partition_coding_and_utr_segments( 
                self.exons, self.cds_region[0], self.cds_region[1] )
        
            us_codon, ds_codon = self.cds_region[0], self.cds_region[1]
            # if this is a reverse strand transcript, rev 5' and 3' ends
            if self.strand == '+':
                self.fp_utr_exons, self.tp_utr_exons \
                    = self.us_exons, self.ds_exons
                self.start_codon, self.stop_codon = us_codon, ds_codon
            else:
                self.fp_utr_exons, self.tp_utr_exons \
                    = self.ds_exons, self.us_exons
                self.stop_codon, self.start_codon = us_codon, ds_codon
        
        # add these for compatability
        #self.append( trans_id )
        #self.append( exon_bnds )        
    
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
        
        ret_lines.extend( build_lines_for_feature( 
                self.exons, 'exon', False ) )
        
        if self.cds_region != None:
            us_exons, ds_exons = self.fp_utr_exons, self.tp_utr_exons
            us_label, ds_label = 'five_prime_UTR', 'three_prime_UTR'
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

def flatten( regions ):
    new_regions = []
    curr_start = regions[0][0]
    curr_end = regions[0][1]
    for i,(start,end) in enumerate(regions):
        if curr_end > end:
            end = curr_end
        if i+1 == len( regions ): 
            if len(new_regions) == 0:
                new_regions = [curr_start, curr_end]
                break
            if new_regions[-1][1] == end:
                break
            else:
                new_regions.append( [ curr_start, curr_end] )
                break
        if regions[i+1][0]-end <= 1:
            curr_end = max( regions[i+1][1], end ) 
        else:
            new_regions.append( [curr_start, curr_end] )
            curr_start = regions[i+1][0]
            curr_end = regions[i+1][1]
    if type(new_regions[0]) == int:
        return [new_regions]
    else:
        return new_regions

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
            
            num_exons = g.transcripts[j].contents.num_exons
            exons = g.transcripts[j].contents.exons[:num_exons]
            if cds_start == -1 or cds_stop == -1:
                assert cds_start == -1 and cds_stop == -1
                cds_region = None
            else:
                cds_region = ( cds_start, cds_stop )
            
            exons = [(x.start, x.stop) for x in exons]
            exons = flatten( exons )
            
            transcripts.append( 
                Transcript(trans_id, chrm, g.strand, 
                           exons, cds_region, g.gene_id, 
                           score=score, rpkm=rpkm, rpk=rpk ) 
                )

        min_loc = min( min(t.exon_bnds) for t in transcripts )
        gene =  Gene( g.gene_id, chrm, g.strand, 
                      min_loc, g.max_loc, transcripts )
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
            print list( gene )
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
