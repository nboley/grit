import sys, os
import re
from time import time, sleep
from urllib import urlencode
import urllib2

# NCBI Batch - Conserved Domain Database Search url
URL = "http://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"

# jobState values from url source
JOB_STATUS = {0:'DONE', 1:'INVALIDRID', 2:'NULLINPUT', 3:'RUNNING', 4:'QMERROR', 5:'DATAERROR', 6:'ABUSEINPUT'}

#vvvvvvvvvvvv VARS FOR OUTPUT FORMATING vvvvvvvvvvvvvv
# should be in [ 'hits', 'aligns', 'feats' ]
# if TDATA = 'feats' only QUERY_DEF_LINE will be used
TDATA_MODE = 'hits'
# should be in [ 'asn', 'xml', 'json', 'text' ]
# only used in TDATA_MODE = 'aligns'
ALN_FORMT = None

# should be in [ 'rep', 'all' ] equivalent to [ 'Concise', 'Full' ] from webpage
MODE = 'rep'

# evuivalent to 'Superfamiy Only' on batch cdd webpage
SUPER_ONLY = False
# evuivalent to 'Include query defline' on batch cdd webpage
QUERY_DEF_LINE = True
# evuivalent to 'Include domain defline' on batch cdd webpage
CDD_DEF_LINE = False
#^^^^^^^^^^^^^^ VARS FOR OUTPUT FORMATING ^^^^^^^^^^^^^^^

def wait_for_server( cdsid, sleep_secs=10 ):
    """wait for NCBI server to complete a given request with cdsid
    """
    wait_data = { 'cdsid': cdsid, 'tdata': 'none', 'callback': 'jsonp' }
    while( True ):
        # update timestamp attribute (not sure why its needed)
        wait_data['_'] = int(time())
        # format GET request url
        waiting_url = URL + "?" + urlencode( wait_data )
        # create GET request for waiting status
        req = urllib2.Request( waiting_url, headers={ 'Content-Type': 'application/json' } )
        
        # parse status from jsonp output
        status = int(re.search( "'status': (\d)", urllib2.urlopen( req ).read() ).group(1))
        # determine what should be done given current job status
        if JOB_STATUS[status] not in ('RUNNING', 'DONE'):
            print "Encountered unexpected status code {0}. Exiting.".format(JOB_STATUS[status])
            sys.exit(1)
        # Server is still running request, so wait for it to finish
        elif JOB_STATUS[status] == 'RUNNING':
            sleep(sleep_secs)
        # once job is done return
        else:
            break
    
    return

def write_output( out_fp, download_cdsid ):
    # Send final POST request for ultimate response
    final_data = { 'cdsid': download_cdsid }
    final_response = urllib2.urlopen( URL, urlencode(final_data) )
    
    # write final output to out_fp
    out_fp.write( final_response.read() )
    out_fp.close()
    
    return

def request_output( initial_cdsid ):
    # server has finished getting cds correctly so download the output
    # build final results download form from 
    results_data = {'cdsid': initial_cdsid, 'callback':'jsonp', '_': int(time()), 'tdata': TDATA_MODE, 'mode': MODE}
    if SUPER_ONLY:
        results_data[ 'clonly' ] = 'true'
    if QUERY_DEF_LINE:
        results_data[ 'qdefl' ] = 'true'
    if CDD_DEF_LINE:
        results_data[ 'cddefl' ] = 'true'
    if TDATA_MODE == 'aligns':
        results_data[ 'alnfmt' ] = ALN_FORMT
    download_url = URL + "?" + urlencode(  results_data )
    
    # format download GET request
    download_req = urllib2.Request( download_url, headers={ 'Content-Type': 'application/json' } )
    
    # get response from server
    download_respose = urllib2.urlopen( download_req ).read()
    
    # parse cdsid from download post request for lookups
    return re.search( "'cdsid': '(.*)'", download_respose ).group(1)

def submit_query( fasta_input ):
    # format initial query data
    query_data = urlencode( {"queries": fasta_input, "enctype": "multipart/form-data" } )
    
    # submit post request to server
    query_response = urllib2.urlopen( URL, query_data )
    
    # parse cdsid for lookups out of initial response body
    return re.search( 'var ctrlHandle = "(.*)";', query_response.read()).group(1)

def run_cdd_search( fasta_fp, out_fp ):
    fasta_input = fasta_fp.read()
    
    if VERBOSE: print 'submitting initial query...'
    initial_cdsid = submit_query( fasta_input )
    
    if VERBOSE: print 'waiting for initial query...'
    # wait for the searver to finish processing request
    # this may take some time for larger requeasts
    wait_for_server( initial_cdsid )
    
    if VERBOSE: print 'requesting output...'
    download_cdsid = request_output( initial_cdsid )
    
    if VERBOSE: print 'waiting for output...'
    # wait for the server to format output and return final data
    wait_for_server( download_cdsid )
    
    if VERBOSE: print 'writing output to file...'
    write_output( out_fp, download_cdsid )
    
    return

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(\
        description='Run Concerverd Domain Database (CDD) searches on a fasta file.' )
    parser.add_argument( 'fasta', type=file, \
            help='Fasta file from which to search for conserved domains.')
    parser.add_argument( '--out_fname', '-o', \
            help='Prefix of output file. (default: <in_fasta>.cdd.txt)')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
            help='Whether or not to print status information.')
    args = parser.parse_args()
    
    # create out_fp (default as described above)
    out_fn = args.out_fname if args.out_fname != None else os.path.basename( args.fasta.name ) + '.cdd.txt'
    out_fp = open( out_fn, 'w' )
    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.fasta, out_fp

def main():
    fasta_fp, out_fp = parse_arguments()
    
    run_cdd_search( fasta_fp, out_fp )
    
    return

if __name__ == '__main__':
    main()

