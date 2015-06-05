"""
Copyright (c) 2011-2015 Nathan Boley

This file is part of GRIT.

GRIT is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GRIT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GRIT.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import sys
import time
import signal
import multiprocessing
import traceback

from grit import config

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
    os._exit(os.EX_TEMPFAIL)

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
                try:
                    signal.signal(signal.SIGINT, handle_interrupt_signal)
                    target(*args)
                    os._exit(os.EX_OK)
                except Exception, inst:
                    config.log_statement( "Uncaught exception in subprocess\n" 
                                          + traceback.format_exc(), log=True)
                    os._exit(os.EX_SOFTWARE)
            else:
                pids.append(pid)
        try:
            while len(pids) > 0:
                ret_pid, error_code = os.wait()
                pids.remove(ret_pid)
                if error_code != os.EX_OK: 
                    raise OSError, "Process '{}' returned error code '{}'".format(
                        ret_pid, error_code) 
        except KeyboardInterrupt:
            for pid in pids:
                try: os.kill(pid, signal.SIGHUP)
                except: pass
            raise
        except OSError:
            for pid in pids:
                try: os.kill(pid, signal.SIGHUP)
                except: pass
            raise
        return
