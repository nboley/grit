import numpy as np
cimport numpy as np

cdef extern from "math.h":
    double log(double)

cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    void* calloc( size_t n, size_t size )
    void free(void *ptr)

cimport cython
@cython.boundscheck(False)
def calc_lhd( np.ndarray[np.double_t, ndim=1] freqs not None, 
              np.ndarray[np.int_t, ndim=1] observed_array not None, 
              np.ndarray[np.double_t, ndim=2] expected_array not None ):
    cdef int num_transcripts = freqs.shape[0]
    cdef int num_bins = expected_array.shape[0]
    # build the expected bin frequencies
    cdef double lhd = 0
    cdef double freq = 0
    cdef int i = 0
    cdef int j = 0
    for i in range(num_bins):
        # calculate this bin's frequency
        freq = 0
        for j in range(num_transcripts):
            freq += freqs[j]*expected_array[i,j]
        
        lhd += observed_array[i]*log(freq)
    
    return lhd

@cython.boundscheck(False)
@cython.cdivision(True)
def calc_lhd_deriv( np.ndarray[np.double_t, ndim=1] freqs not None, 
                    np.ndarray[np.int_t, ndim=1] observed_array not None, 
                    np.ndarray[np.double_t, ndim=2] expected_array not None ):
    cdef int num_transcripts = freqs.shape[0]
    cdef int num_bins = expected_array.shape[0]

    cdef int i = 0
    cdef int j = 0
    cdef double* weights = <double *>calloc( num_bins, sizeof( double ) )
    cdef double freq
    for i in range(num_bins):
        # calculate this bin's frequency
        freq = 0
        for j in range(num_transcripts):
            freq += freqs[j]*expected_array[i,j]
        weights[i] = observed_array[i]/freq
    
    gradient = np.zeros( num_transcripts, dtype=np.double )
    cdef double curr_grad_value
    for i in range(num_transcripts):
        curr_grad_value = 0
        for j in range(num_bins):
            curr_grad_value += weights[j]*expected_array[j,i]
        gradient[i] = curr_grad_value
    
    return -gradient
