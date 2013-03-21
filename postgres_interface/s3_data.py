import os, sys
import multiprocessing
import boto

from build_fl_dist_from_db import build_fl_dist

VERBOSE = False
DEBUG_VERBOSE = False

# time to sleep before cehcking download
TIME_TO_SLEEP_BEFORE_RETRIES = 1.0

VALID_STATUSES = [ 'FINISHED', 'DOWNLOADING', 'ERROR' ]

class S3Cache(object):
    """A process safe mapping from S3 files to local files.

    This returns the filename of a local copy of an s3 file, downloading and
    verifying if necessary.
    """
    def __init__( self, access_key, private_key, db_conn, scratch_directory="/scratch/" ):
        self._access_key = access_key
        self._private_key = private_key
        self.conn = boto.connect_s3( self._access_key, self._private_key )
        self.db_conn = db_conn
        
        self.scratch_directory = scratch_directory
        
        self._manager = multiprocessing.Manager()
        self._lock = multiprocessing.Lock()
        self._location_to_fn_mapping = self._manager.dict()
    
    @staticmethod
    def get_bucket_and_key_from_s3_url( url ):
        bucket = url[5:].split("/")[0]
        key = "/".join(url[5:].split("/")[1:])
        return bucket, key
    
    def find_local_filename_from_s3_location( self, bucket, key ):
        return os.path.join( 
            self.scratch_directory, 
            ''.join( x for x in bucket.replace("/", "_") 
                     if x.isalnum() or x == '_' ) \
            + ''.join( x for x in key.replace("/", "_") 
                     if x.isalnum() or x == '_' ) + ".bam"
            )

    def _get_previously_cached_file( self, bucket_name, file_key ):
        while True:
            local_fn, status \
                = self._location_to_fn_mapping[(bucket_name, file_key)]
            assert status in VALID_STATUSES
            # if it does, and it's finished, awesome, return it
            if status == 'FINISHED':
                return local_fn
            # if it does, but there's an error associated with it, then 
            # raise an exception
            elif status == 'ERROR':
                raise ValueError, "Cant get file %s %s" % ( bucket_name, file_key )
            # otherwise, wait for it to download
            else:
                assert status == 'DOWNLOADING'
                if DEBUG_VERBOSE: 
                    print "Waiting on %s: status %s" % ( local_fn, status )
                time.sleep(TIME_TO_SLEEP_BEFORE_RETRIES)

        assert False
    
    def _add_file_to_cache( self, bucket_name, file_key, local_fn ):
        bucket = boto.s3.bucket.Bucket( self.conn, bucket_name )
        key = bucket.lookup( file_key )
        
        # get the bam file
        if DEBUG_VERBOSE: print "Downloading ", file_key
        key.key = file_key
        key.get_contents_to_filename(local_fn)
        
        # get the bai file
        if DEBUG_VERBOSE: print "Downloading ", file_key + '.bai'
        key.key = file_key + '.bai'
        key.get_contents_to_filename(local_fn + '.bai')
        
        # build the fl dist
        if not os.path.exists( local_fn + ".fldist" ):
            if DEBUG_VERBOSE: print "Building ", file_key + '.fldist'
            build_fl_dist(local_fn, self.db_conn)
        
        # release the lock
        self._location_to_fn_mapping[(bucket_name, file_key)] = (
            local_fn, 'FINISHED' )        
        
        return
    
    def open_local_copy( self, bucket_name, file_key ):
        """Returns the filename of the cached object.
        
        """
        local_fn = self.find_local_filename_from_s3_location( bucket_name, file_key )
        
        # open a lock to make sure a file doesn't get added between deciding we need
        # to add it, and actually starting to add it
        self._lock.acquire()        
        
        # see if the file exists in the manager. If it does, then return it or raise
        # an exception
        if local_fn in self._location_to_fn_mapping:
            self._lock.release()
            local_fn = self._get_previously_cached_file( self, bucket_name, file_key )
        # otherwise, see if the filename exists locally, and verify it's complete with a checksum
        # we assume the bam index and fragment length file exist
        elif os.path.exists( local_fn ):
            bucket = boto.s3.bucket.Bucket( self.conn, bucket_name )
            key = bucket.lookup( file_key )
            
            # calculate the checksum
            local_size = os.stat( local_fn ).st_size
            amazon_size = key.size
            
            # if they are the same, then add it into the db
            if DEBUG_VERBOSE: 
                print "Local size: %s\tAmazon size: %s\tDifference: %s" % (
                    local_size, amazon_size, amazon_size - local_size )
            
            if local_size == amazon_size:
                self._location_to_fn_mapping[(bucket_name, file_key)] = (
                    local_fn, 'FINISHED' )        
                self._lock.release()
            else:
                # remove the local file
                #os.remove(local_fn)
                #os.remove(local_fn + '.bami')
                #os.remove(local_fn + 'fldist')
                self._location_to_fn_mapping[(bucket_name, file_key)] = (
                    local_fn, 'DOWNLOADING' )
                self._lock.release()
                self._add_file_to_cache( bucket_name, file_key, local_fn )
            
        # if it's not in the mapping and doesn't seem to exists locally,
        # then we need to download it. So do that.
        else:
            self._location_to_fn_mapping[(bucket_name, file_key)] = (
                local_fn, 'DOWNLOADING' )
            self._lock.release()
            self._add_file_to_cache( bucket_name, file_key, local_fn )
        
        return open( local_fn )

    def open_local_copy_from_s3_url( self, url ):
        bucket, key = self.get_bucket_and_key_from_s3_url( url )
        return self.open_local_copy( bucket, key )
        
    
def main():
    global VERBOSE
    VERBOSE = True
    global DEBUG_VERBOSE
    DEBUG_VERBOSE = True
    
    access_key = sys.argv[1]
    private_key = sys.argv[2]
    s3_url = sys.argv[3]
    
    cache = S3Cache( access_key, private_key, None )
    fp = cache.open_local_copy_from_s3_url( s3_url )
    print fp

if __name__ == '__main__':
    main()
