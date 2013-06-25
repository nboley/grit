import os
import time
import curses
import multiprocessing

def manage_curses_display(stdscr, msg_queue, msg_queue_lock):
    curses.curs_set(0)    
    stdscr.timeout(0)
    header = stdscr.derwin(2, 80, 0, 0)

    thread_data = stdscr.derwin(NTHREADS+2, 80, 3, 0)
    thread_data.insstr( 0, 0, "Main:".ljust(11) )
    for i in xrange(NTHREADS):
        thread_data.insstr( i+1, 0, ("Thread %i:" % (i+1)).ljust(11) )
    
    stdscr.addstr(NTHREADS+2+2+1, 0, "Log:" )
    
    log_border = stdscr.derwin(20+2, 80, NTHREADS+2+2+2, 0)
    log_border.border()
    log = log_border.derwin( 20, 78, 1, 1 )
    
    header.addstr(0, 0, "GRIT (version %s)" % (__version__, ) )
    stdscr.refresh()    
    while True:
        try:
            msg_queue_lock.acquire()
            thread_index, msg = msg_queue.pop()
        except IndexError, inst:
            msg_queue_lock.release()
            time.sleep(0.1)
            continue

        if thread_index in (0, None ):
            log.insertln()
            log.insstr( msg )
            log.refresh()
        
        if 0 == thread_index:
            thread_data.addstr(0, 0, "Main:".ljust(13) + msg.ljust(80-13) )
            thread_data.refresh()            
        elif thread_index != None:
            thread_data.addstr(thread_index, 0, ("Thread %i:" % (thread_index)).ljust(13) + msg.ljust(80-13) )
            thread_data.refresh()

        stdscr.refresh()        
        msg_queue_lock.release()
        if msg == 'BREAK': break
    
    return

class Logger( object ):
    def __init__(self, nthreads, log_fname="log.txt"):
        self.manager = multiprocessing.Manager()
        self.ofp = open( log_fname, "w" )
        self.msgs_lock = self.manager.Lock()
        self.msgs = self.manager.list()
        p = multiprocessing.Process( target=curses.wrapper, 
                     args=(manage_curses_display, self.msgs, self.msgs_lock) )
        p.start()
        self.curses_p = p
        self.main_pid = os.getpid()

        self.pid_to_index_mapping = self.manager.list()
        self.pid_to_index_mapping.append( self.main_pid )
        for loop in xrange(nthreads):
            self.pid_to_index_mapping.append( None )
        
        return
    
    def __call__( self, message, ONLY_LOG=False ):
        self.ofp.write(message.strip() + "\n" )
        self.ofp.flush()
        
        self.msgs_lock.acquire()
        if ONLY_LOG: 
            p_index = None
        else:
            try: 
                p_index = self.pid_to_index_mapping.index( os.getpid() )
            except ValueError:
                p_index = min( i for i, pid in enumerate(self.pid_to_index_mapping) 
                               if pid == None or not os.path.exists("/proc/%i" % pid) )
                self.pid_to_index_mapping[p_index] = os.getpid()
        
        self.msgs.append( (p_index, message) )
        self.msgs_lock.release()
        # make sure that the message has time to be displayed
        time.sleep(0.2)
    
    def close(self):
        self.ofp.close()
        self.msgs.append( (None, 'BREAK') )
        self.curses_p.join()

