import os
import time
import signal
import multiprocessing

class ThreadSafeFile( file ):
    def __init__( *args ):
        file.__init__( *args )
        args[0].lock = multiprocessing.Lock()

    def write( self, string ):
        with self.lock:
            file.write( self, string )
            self.flush()

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

def handle_interrupt_signal(signum, frame):
    os._exit(0)

def fork_and_wait(n_proc, target, args):
    """Fork n_proc processes, run target(*args) in each, and wait to finish.
    
    """
    if n_proc == 1:
        target(*args)
        return
    else:
        pids = []
        for i in xrange(n_proc):
            pid = os.fork()
            if pid == 0:
                signal.signal(signal.SIGINT, handle_interrupt_signal)
                target(*args)
                os._exit(0)
            else:
                pids.append(pid)
        try:
            for pid in pids: 
                ret_pid, error_code = os.waitpid(pid, 0)
                assert pid == ret_pid
                if error_code != 0: 
                    raise OSError, "Process '{pid}' returned error code '%i'".format(
                        pid, error_code) 
                else:
                    pids.remove(pid)
        except KeyboardInterrupt:
            for pid in pids:
                os.kill(pid, signal.SIGINT)
            raise
        except OSError:
            for pid in pids:
                os.kill(pid, SIGINT)
            raise            
        return
