# Copyright (c) 2011-2012 Nathan Boley

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
    
    res = []
    
    c_contigs = cast( c_contigs_p, POINTER(c_contigs_t) ).contents
    for i in xrange( c_contigs.size ):
        values = c_contigs.contigs[i].values
        size = c_contigs.contigs[i].size
        name = c_contigs.contigs[i].name
        
        if 0 == size:
            print >> sys.stderr, "WARNING: found a 0 length chrm '%s'. Skipping it." % name
            continue
        
        array = numpy.ctypeslib.as_array( values, (size,))

        """
        try:
            assert( not numpy.isnan( array.sum() ) )
        except:
            print fname, size, name
            print numpy.where( numpy.isnan(array) )
            print array[ 1710:1720 ]
            print values[ 1710:1720 ]
            print c_contigs.contigs[i].values[ 1710:1720 ]
            raise
        """
        
        res.append( ( name, array ) )
    
    return res

if __name__ == "__main__":
    for name, track in iter_bedgraph_tracks( sys.argv[1] ):
        print name, track, track.min(), track.sum()

# WAT?
