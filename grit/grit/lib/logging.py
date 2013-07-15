import os, sys
import time
import curses
import multiprocessing

__version__ = "0.1.1"
MAX_REFRESH_TIME = 1e-2
MAX_NCOL = 120

def manage_curses_display(stdscr, msg_queue, msg_queue_lock, nthreads=1):
    curses.curs_set(0)    
    stdscr.timeout(0)
    header = stdscr.derwin(2, MAX_NCOL, 0, 0)

    thread_data_windows = []
    thread_data_windows.append( stdscr.derwin(1, MAX_NCOL, 3, 0) )
    thread_data_windows[-1].insstr( 0, 0, "Thread 0:".ljust(11) )
    for i in xrange(nthreads):
        thread_data_windows.append( stdscr.derwin(1, MAX_NCOL+i, 3+i+1, 0))
        thread_data_windows[i+1].insstr(0, 0, ("Thread %i:" % (i+1)).ljust(11))

    stdscr.addstr(nthreads+2+2+1, 0, "Log:" )

    log_border = stdscr.derwin(20+2, MAX_NCOL, nthreads+2+2+2, 0)
    log_border.border()
    log = log_border.derwin( 20, 78, 1, 1 )

    header.addstr(0, 0, "GRIT (version %s)" % (__version__, ) )
    stdscr.refresh()

    while True:
        start_time = time.time()
        while True:
            # make sure that we refresh every one in a while
            if start_time - time.time() > MAX_REFRESH_TIME:
                break

            # try to acquire a message. If none exists, brak
            # to refresh and sleep
            try:
                msg_queue_lock.acquire()
                if len( msg_queue ) == 0:
                    msg_queue_lock.release()
                    break
                else:
                    thread_index, do_log, msg = msg_queue.pop()
            except IOError, inst:
                break
            
            # release the message queue lock
            msg_queue_lock.release()

            # if the message is BREAK, then we are done so exit the thread
            if msg == 'BREAK': 
                return

            if do_log:
                log.insertln()
                log.insstr( msg )
                log.refresh()
            
            # truncate the message so that it doesnt extend past 80 charcters
            msg = msg[:MAX_NCOL-11]
            if thread_index != None:
                line = ("Thread %i:" % (thread_index)).ljust(11) \
                    + msg.ljust(MAX_NCOL-11)
                thread_data_windows[thread_index].erase()
                thread_data_windows[thread_index].insstr(0, 0, line )
                thread_data_windows[thread_index].refresh()
        
        stdscr.refresh()
    
    return

class Logger( object ):
    def _init_ncurses_manager(self):
        self.manager = multiprocessing.Manager()
        self.msgs_lock = self.manager.Lock()
        self.msgs = self.manager.list()
        p = multiprocessing.Process( target=curses.wrapper, 
                     args=(manage_curses_display, 
                           self.msgs, self.msgs_lock, self.nthreads) )
        p.start()
        self.curses_p = p
        self.main_pid = os.getpid()

        self.pid_to_index_mapping = self.manager.list()
        self.pid_to_index_mapping.append( self.main_pid )
        for loop in xrange(self.nthreads):
            self.pid_to_index_mapping.append( None )


    def __init__(self, nthreads, use_ncurses=False, log_ofstream=None):
        self.use_ncurses = use_ncurses
        self.nthreads = nthreads
        self.log_ofstream = log_ofstream
        
        # if we're not using ncurses and the log_ofstream is not set,
        # then by default write messages to stderr
        if not self.use_ncurses and self.log_ofstream == None:
            self.log_ofstream = sys.stderr
        
        if self.use_ncurses:
            self._init_ncurses_manager()
        
        return
    
    def __call__( self, message, only_log=False, do_log=None ):
        if self.log_ofstream != None:
            self.log_ofstream.write(message.strip() + "\n" )
            self.log_ofstream.flush()

        if self.use_ncurses:
            self.msgs_lock.acquire()
            try: 
                p_index = self.pid_to_index_mapping.index( os.getpid() )
            except ValueError:
                p_index = min( i for i, pid in enumerate(self.pid_to_index_mapping) 
                               if pid == None or not os.path.exists("/proc/%i"%pid))
                self.pid_to_index_mapping[p_index] = os.getpid()

            # only log message from main
            if do_log == None: do_log = True if p_index == 0 else False
            self.msgs.append( (p_index, only_log or do_log, message) )
            self.msgs_lock.release()
        
        time.sleep(1e-2)
    
    def close(self):
        if self.use_ncurses:
            self.msgs.append( (None, False, 'BREAK') )
            self.curses_p.join()
