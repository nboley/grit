import multiprocessing

class ProcessSafeOPStream( object ):
    def __init__( self, writeable_obj ):
        self.writeable_obj = writeable_obj
        self.lock = multiprocessing.Lock()
        self.name = self.writeable_obj.name
        return
    
    def write( self, data ):
        self.lock.acquire()
        self.writeable_obj.write( data )
        self.writeable_obj.flush()
        self.lock.release()
        return
    
    def close( self ):
        self.writeable_obj.close()

class Pool(object):
    """A working version of the multiprocessing's Pool.
    
    """
    def __init__(self, nthreads):
        self.nthreads = nthreads
        self.processes = [None]*nthreads
    
    def apply( self, fn, all_args ):
        while len(all_args) > 0:
            if all( p != None and p.is_alive() for p in self.processes ):
                time.sleep(0.1)
                continue
            else:
                proc_i = min( i for i, p in enumerate(self.processes) 
                              if p == None or not p.is_alive() )
                args = all_args.pop()
                self.processes[proc_i] = multiprocessing.Process(
                    target=fn, args=args)
                self.processes[proc_i].start()
            
        for i, p in enumerate(self.processes):
            if p != None: 
                p.join()
                self.processes[i] = None
        
        return
