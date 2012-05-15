import sys, os
import numpy
import signal
from ctypes import *
import numpy.ctypeslib

VERBOSE = False

# load the gtf library
bedgraph_o = cdll.LoadLibrary( os.path.join( ".", \
        os.path.dirname( __file__), "libbedgraph.so" ) )

class c_contig_t(Structure):
    """
struct contig_t {
    char *name;
    int size;
    double* values;
};
"""
    _fields_ = [("name", c_char_p),
                ("size", c_int),
                ("values", POINTER(c_double))
               ]

class c_contigs_t(Structure):
    """
struct contigs_t {
    struct contig_t* contigs;
    int size;
};
"""
    _fields_ = [
                 ("contigs", POINTER(c_contig_t)),
                 ("size", c_int),
               ]

def iter_bedgraph_tracks( fname ):
    c_contigs_p = c_void_p()
    bedgraph_o.load_bedgraph( c_char_p(fname), byref(c_contigs_p) )

    # check for a NULL pointer
    if not bool(c_contigs_p):
        raise IOError, "Couldn't load '%s'" % fname
    
    c_contigs = cast( c_contigs_p, POINTER(c_contigs_t) ).contents
    for i in xrange( c_contigs.size ):
        values = cast( c_contigs.contigs[i].values, POINTER(c_double) )
        array = numpy.ctypeslib.as_array( values, (c_contigs.contigs[i].size,))
        yield c_contigs.contigs[i].name, array
    
    return

if __name__ == "__main__":
    for name, track in iter_bedgraph_tracks( sys.argv[1] ):
        print name, track

# WAT?
