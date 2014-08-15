from libc.math cimport log, exp, lgamma

def calc_moments(double p, int n, double TOL=1e-6):
    cdef double mean = 0
    cdef double var = 0
    cdef int x
    cdef double prob, value, mean_update, var_update
    
    cdef double log_n_fac, log_x_fac, log_n_mx_fac
    log_n_fac = lgamma(n+1)
    log_n_mx_fac = log_n_fac
    log_x_fac = 0
    
    cdef double log_p, log_1m_p
    log_p = log(p)
    log_1m_p = log(1-p)
        
    for x in xrange(1,n+1):
        log_x_fac += log(x)
        log_n_mx_fac -= log(n-x+1)
        
        # calculate the binomial probability
        prob = log_n_fac - log_x_fac - log_n_mx_fac
        prob += x*log_p
        prob += (n-x)*log_1m_p
        prob = exp(prob)
        
        # calculate the value of the statistic for this x
        value =  x*log_p - log_x_fac
        mean_update = prob*value
        mean += mean_update
        var_update = mean_update*value
        var += var_update
        if x > p*n + 1 and abs(mean_update)*(n-x) < TOL and abs(var_update)*(n-x) < TOL: 
            break
    
    return mean, var - mean**2
