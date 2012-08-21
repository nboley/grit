# Copyright (c) 2011-2012 Nathan Boley

import sys
import os

import numpy
import re
from collections import defaultdict, namedtuple
from itertools import izip, product
from copy import deepcopy

from chrm_sizes import ChrmSizes

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "./fast_wiggle_parser/" ) )
from bedgraph import iter_bedgraph_tracks

VERBOSE = False


GenomicInterval = namedtuple('GenomicInterval', \
                                 ['chr', 'strand', 'start', 'stop'])


def threshold_and_scale_read_cvg_from_nc_density( \
        cvg, nc_wig, smoothing_window_size, ratio_threshold, \
        intermediate_wiggles_prefix=None ):
    if VERBOSE: print >> sys.stderr, '\tSmoothing coverage wiggle.'        
    # get a smoothed copy of the tes coverage wiggle
    smoothed_cvg = cvg.get_smoothed_copy( smoothing_window_size )

    if VERBOSE: print >> sys.stderr, '\tSmoothing negative control wiggle.'
    
    # smooth the nc wiggle
    smoothed_nc_wig = nc_wig.get_smoothed_copy( smoothing_window_size )

    for key in smoothed_cvg.iterkeys():
        # wiggle arrays were made from same genome_data file so thier sizes 
        # should be the same
        assert len( nc_wig[key] ) == len( smoothed_cvg[ key ] )
        
        # rescale the NC array to be on the same scale, noting that it
        # has been previously normalized 
        smoothed_nc_wig[key] = smoothed_nc_wig[key] \
            * smoothed_cvg[ key ].sum()
    
    cvg.threshold_cvrg( smoothed_cvg, smoothed_nc_wig, ratio_threshold )
    if intermediate_wiggles_prefix != None:
        if VERBOSE: print >> sys.stderr, '\tWriting thresholded cvg wiggles'
        track_name_prefix = os.path.basename( intermediate_wiggles_prefix ) \
            + ".thresholded"
        track_name_prefix = track_name_prefix.replace( ".", "_" )
        cvg.write_wiggles( \
            intermediate_wiggles_prefix + '.thresholded.plus.wig',  \
            intermediate_wiggles_prefix + '.thresholded.minus.wig', \
            ignore_zeros=True, track_name_prefix=track_name_prefix )

    cvg.scale_cvrg( smoothed_nc_wig )
    if intermediate_wiggles_prefix != None:
        if VERBOSE: print >> sys.stderr, '\tWriting scaled tes wiggles...'
        track_name_prefix = os.path.basename( intermediate_wiggles_prefix ) \
            + ".thresholded"
        track_name_prefix = track_name_prefix.replace( ".", "_" )
        cvg.write_wiggles( intermediate_wiggles_prefix + '.scaled.plus.wig', \
                           intermediate_wiggles_prefix + '.scaled.minus.wig',\
                           ignore_zeros=True, \
                           track_name_prefix=track_name_prefix )
    
    return cvg

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 
           'bartlett', 'blackman'. 
        flat window will produce a moving average smoothing.
            
    output:
        the smoothed signal
        
    example:
    
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    """
    
    # make window_len odd so that convolution does not change the 
    # length of the array
    if window_len % 2 == 0:
        window_len -= 1
    
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', " \
            + "'bartlett', 'blackman'"
    
    
    # mirror the first and last elements around the beginning and the end 
    # of the array to minimize the fringe effects
    # Note that the sum of the smoothed array is not garunteed to be equal to the
    # input array, but its length will be the same
    s=numpy.r_[x[int(window_len/2):0:-1],x,x[-1:-(int(window_len/2)+1):-1]]
    
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')
    
    y=numpy.convolve(w/w.sum(),s,mode='valid')
    
    return y

def parse_bed_like_wig_line( line, strand, chrm=None ):
    data = line.split()
    if len(data) != 4: return None, None
    chrm = data[0]
    if chrm.startswith( 'chr' ): chrm = chrm[ 3:]
    try:
        start = int(data[1])
        stop = int(data[2])
        # convert bed style 0-based closed-open to gtf style 1-based closed-closed
        start += 1
        value = float(data[3])
    except ValueError:
        return None, None
    
    return GenomicInterval( chrm, strand, start, stop ), value

def parse_var_step_wig_line( line, strand, chrm ):
    # TODO:::add span varible as per variableStep formatting
    data = line.split()
    if len(data) != 2: return None, None
    try:
        pos = int(data[0])
        value = float(data[1])
    except ValueError:
        return None, None
    
    return GenomicInterval( chrm, strand, pos, pos ), value

class Wiggle( dict ):
    @staticmethod
    def _infer_strand_from_fname( fname ):
        if fname.lower().find( "plus" ) >= 0:
            return '+'
        elif fname.lower().find( "+" ) >= 0:
            return '+'
        elif fname.lower().find( "minus" ) >= 0:
            return '-'
        elif fname.lower().find( "-" ) >= 0:
            return '-'
        else:
            raise ValueError, "Couldn't infer strand from filename '%s'" % fname
        
        return
    
    def write_wiggles( self, plus_out, minus_out, \
                           ignore_zeros=True, track_name_prefix=None, \
                           filter_region=None ):
        if filter_region != None:
            assert isinstance( filter_region, GenomicInterval )
            if filter_region.strand == '+':
                minus_out = None
            if filter_region.strand == '-':
                plus_out = None
        
        def open_ofp( out ):
            if out == None: return None
            
            if isinstance( out, str ):
                return open( out, "w" )
            elif hasattr( plus_out, "write" ):
                return open( out, "w" )
            else:
                raise TypeError, "Write wiggles takes either a filename, " \
                    + "or a stream with a write method."

        # open the output file pointers, if necessary
        plus_out_fp = open_ofp( plus_out )
        minus_out_fp = open_ofp( minus_out )
        
        if track_name_prefix != None:
            # make sure the track name is valid
            track_name_prefix = track_name_prefix.replace( ".", "_" )
            if plus_out != None:
                plus_out_fp.write( "track type=bedGraph name={0}_plus\n".format( \
                        track_name_prefix ) )
            if minus_out != None:
                minus_out_fp.write( "track type=bedGraph name={0}_minus\n".format( \
                        track_name_prefix ) )
        
        if self.intervals == None:
            self.calc_all_intervals()
        
        for key in self.iterkeys():
            chrm, strand = key
            if filter_region != None:
                # we check both the tail and chr to compare properly
                # between chr4 and 4 ( or 'chr2L' and '2L' )
                if filter_region.chr != chrm and filter_region.chr[3:] != chrm:
                    continue
                if filter_region.strand not in ( strand, '.' ):
                    continue
            
            fp = plus_out_fp if strand == '+' else minus_out_fp
            for start, stop, val in self.intervals[key]:
                if ignore_zeros and val == 0:
                    continue
                if filter_region != None:
                    # TODO - use a bisect for this
                    if stop < filter_region.start:
                        continue
                    if start > filter_region.stop:
                        break
                
                # subtract one from start to make 0-based closed-open (bed format)
                fp.write( '\t'.join( \
                        ('chr' + chrm, str(int(start-1)), str(int(stop)), \
                             str(val)) ) + '\n' )
        
        # if we opened the file pointers, close them
        if isinstance( plus_out, str ):
            plus_out_fp.close()
        if isinstance( minus_out, str ):
            minus_out_fp.close()
        
        return
    
    def calc_all_intervals( self ):
        """Calculate the intervals from the coverage arrays
        """
        self.intervals = {}
        
        for key, coverage in self.iteritems():
            chrm_intervals = []
            curr_start = 0
            curr_val = coverage[1]
            for pos, val in enumerate(coverage[1:]):
                if val != curr_val:
                    # store 1-based closed-closed interval and value
                    chrm_intervals.append( (curr_start+1, pos, curr_val) )
                    curr_start = pos
                    curr_val = val
            
            if curr_start != pos:
                # store the final interval
                chrm_intervals.append( (curr_start+1, pos, curr_val) )
            
            # fix first interval to start at 1 instead of 0
            if chrm_intervals[0][1] == 0:
                del chrm_intervals[0]
            else:
                chrm_intervals[0] = ( 1, chrm_intervals[0][1], chrm_intervals[0][2] )
            
            self.intervals[key] = numpy.array( chrm_intervals )
        
        return
    
    def calc_zero_intervals( self, min_size=1 ):
        """Calculate all intervals of zero coverage over min_size
        """
        self.zero_intervals = {}
        
        for key, coverage in self.iteritems():
            chrm_zero_intervals = []
            non_zero_indices = coverage.nonzero()[0]
            
            # if there is no coverage in this chrom, strand
            if len( non_zero_indices ) == 0:
                self.zero_intervals[key] = numpy.array( \
                    [( 1, self.chrm_sizes[key[0]] )] )
                continue
            
            # if first interval is valid add it
            if non_zero_indices[0] > min_size:
                chrm_zero_intervals.append( ( 1, non_zero_indices[0]-1 ) )
            
            for index, next_index in \
                    zip(non_zero_indices[:-1], non_zero_indices[1:]):
                if next_index - index > min_size:
                    chrm_zero_intervals.append( (index+1, next_index-1) )
            
            # add interval up to the end of the chrom if valid
            if self.chrm_sizes[key[0]] - non_zero_indices[-1] >= min_size:
                chrm_zero_intervals.append( \
                    (non_zero_indices[-1]+1, self.chrm_sizes[key[0]] ) )
            
            self.zero_intervals[key] = numpy.array(chrm_zero_intervals)
        
        return
    
    def get_region_w_fraction_cvrg( self, key, start, stop, fraction, 
                                    from_upstrm=True ):
        """Find the subregion with the specified fraction of coverage

        Return the full region if there is zero coverage.
        """
        region = self[key][start:stop+1]
        # if we want to sum from the downstream position, then 
        # reverse the reigon. 
        tot_cvg = region.sum() 
        if tot_cvg == 0:
            return ( start, stop )

        thresh = tot_cvg * fraction
        
        if not from_upstrm:
            region = region[::-1]
        
        pos_from_start = numpy.cumsum( region ).searchsorted( thresh )
        
        if from_upstrm:
            return ( start, start + pos_from_start )
        
        return ( stop - pos_from_start, stop )
    
    def reset_intervals( self ):
        # reset itervals to None after coverage values have changed
        self.intervals = None
        self.zero_intervals = None
        
        return
    
    def apply( self, fn, args=[] ):
        rv = {}
        args.insert( 0, None )
        for key in self:
            args[0] = self[key]
            rv[key] = fn( *args )
        
        del args[:]
        
        return rv
    
    def update( self, fn, args=[] ):
        for key, val in self.apply( fn, args ).iteritems():
            self[key] = val
    
    def smooth( self, window_len, normalize=True ):
        self.update( smooth, [window_len,] )
        
        def normalize_fn( array ):
            array_sum = array.sum()
            if array_sum < 1e-6:
                return array
            return array / array_sum

        if normalize:
            self.update( normalize_fn )
        
        self.reset_intervals()
        
        return
    
    def get_smoothed_copy( self, window_len, normalize=True ):
        smoothed_self = deepcopy( self )
        smoothed_self.smooth( window_len, normalize=normalize )

        return smoothed_self
    
    def threshold_cvrg( self, smoothed_self, other_wig, ratio_thresh ):
        """threshold coverage with another wiggle object and a ratio_threshold
        other wiggle nust have same keys and be the same size
        i.e. created from the same chrm_sizes_file
        """
        for key in self:
            assert len( other_wig[key] ) == len( self[ key ] ) and \
                len( smoothed_self[key] ) == len( self[ key ] )
            self[ key ][ smoothed_self[key] < other_wig[key]*ratio_thresh ] = 0
        
        self.reset_intervals()
        return
    
    def scale_cvrg( self, other_wig ):
        """scale coverage with another wiggle object 
        other wiggle nust have same keys and be the same size
        i.e. created from the same chrm_sizes_file
        """
        for key in self:
            # add one to avoid divide by zero
            self[key] = self[key] / (other_wig[key] + 1)
        
        self.reset_intervals()
        return

    def add_cvg( self, chrm, strand, start, stop, value ):
        try:
            self[(chrm, strand)][ start : stop + 1 ] += value
        except IndexError:
            print >> sys.stderr, \
                'WARNING: Region of coverage past end of chromosome.'
            print >> sys.stderr, \
                '\tchr{0}:{1:d}-{2:d}:{3}'.format( chrm, strand, start, stop )
        return
    
    def add_cvg_from_wiggle( self, wig_fp, strand ):
        parse_wig_line = parse_bed_like_wig_line
        chrm = None
        get_chr_pattern = re.compile( 'variableStep\s+chrom=chr(\w+)' )
        for line in wig_fp:
            if line.startswith( 'track' ): continue
            # handle both bed style and variableStep type wig files
            if line.startswith( 'variableStep' ):
                parse_wig_line = parse_var_step_wig_line
                try:
                    # TODO:::Handle span argument here (uses default val of 1 now)
                    chrm = get_chr_pattern.match( line ).group(1)
                except AttributeError:
                    print >> sys.stderr, "WARNING: Invalid variableStep wig line " \
                        + "(file may not be parsed correctly): '" + line + "'"
                    continue
            
            region, value = parse_wig_line( line, strand, chrm )
            if region == None: continue
            self.add_cvg( region.chr, region.strand, \
                          region.start, region.stop, \
                          value )
        
        return

    def add_cvg_from_array( self, array, chrm_name, strand ):
        if numpy.isnan(numpy.sum(array)):
            print >> sys.stderr,"WARNING: Detected NaNs in 'add_cvg_from_array'"
            print >> sys.stderr,"NaN locs:", numpy.where( numpy.isnan(array) )
            assert not numpy.isnan(numpy.sum(array))
        
        # strip the leading 'chr' for UCSC compatability
        if chrm_name.startswith( 'chr' ):
            chrm_name = chrm_name[3:]

        if chrm_name not in self.chrm_sizes:
            print >> sys.stderr, "WARNING: '%s' does not" % chrm_name \
                + "exist in the chrm lens file. Skipping this contig." 
            return

        key = ( chrm_name, strand )
        contig_size = self.chrm_sizes[ chrm_name ]
        if len( array ) > contig_size:
            print >> sys.stderr, "WARNING: Values extend past the end of " \
            + "the '%s'. Truncating the array from %i to %i." \
            % ( chrm_name, len(array), contig_size )
            array = array[:contig_size]

        self[ key ][ :len(array) ] += array
        

    def add_cvg_from_bedgraph( self, fname, strand ):
        for chrm_name, array in iter_bedgraph_tracks( fname ):
            if numpy.isnan(numpy.sum(array)):
                print >> sys.stderr,"WARNING: Detected NaNs in 'add_cvg_from_bedgraph'"
                print >> sys.stderr,"NaN locs:", numpy.where( numpy.isnan(array) )
                array[ numpy.isnan(array) ] = 0

            if array.min() < -1e50:
                print >> sys.stderr,"WARNING: Detected Very Small Numbers ( %e ) in 'add_cvg_from_bedgraph'" % array.min()
                print >> sys.stderr,"Small num locs:", numpy.where( array.min() < -1e50 )
                assert False

            if array.min() > 1e50:
                print >> sys.stderr,"WARNING: Detected Very Large Numbers ( %e ) in 'add_cvg_from_bedgraph'" % array.max()
                print >> sys.stderr,"Big Num locs:", numpy.where( array.min() > 1e50 )
                assert False
            
            self.add_cvg_from_array( array, chrm_name, strand )
        
        return

    @staticmethod
    def _fp_is_bedgraph( wig_fp ):
        if 'bedgraph' == wig_fp.name.split( "." )[-1].lower():
            return True
        
        wig_fp.seek(0)
        data = wig_fp.readline().split()
        if len( data ) == 0: 
            return False
        
        if not data[0].lower() == 'track':
            wig_fp.seek(0)
            return False
        
        meta_data = dict( item.split("=") for item in data[1:] )
        if 'type' not in meta_data:
            wig_fp.seek(0)
            return False
        
        if meta_data['type'].lower() != 'bedgraph':
            wig_fp.seek(0)
            return False
        
        wig_fp.seek(0)
        return True
    
    def load_data_from_fp( self, wig_fp, strand="infer_from_fname" ):
        # if we don't know the strand, try and infer it from
        # the filename
        if strand == "infer_from_fname":
            strand = self._infer_strand_from_fname( wig_fp.name )
        else:
            if strand not in "+-":
                raise ValueError, "Unrecognized strand '{0}'".format( strand )
            
        if self._fp_is_bedgraph( wig_fp ):        
            fname = wig_fp.name
            self.add_cvg_from_bedgraph( fname, strand )
        else:
            self.add_cvg_from_wiggle( wig_fp, strand )
        
        self.reset_intervals()
        
        return
    
    def load_data_from_positions( self, positions ):
        """data should be a dict of arrays of positions with coverage
        """
        keys = set( product( self.chrm_sizes, "+-" ) ).intersection( positions )
        for key in keys:
            # load the unsmoothed read positions into the wiggle object
            for pos in positions[ key ]:
                try:
                    self[ key ][ pos ] += 1
                except IndexError:
                    print >> sys.stderr, \
                        'WARNING: Region of coverage past end of chromosome.'
                    print >> sys.stderr, \
                        '\tchr{0[0]}:{1:d}-{1:d}  {0[1]}'.format( key, pos )
                    continue
        
        self.reset_intervals()
        return
    
    def __getitem__( self, key ):
        """If this contig hasn't been accessed yet, initialize it. 

        """
        if key not in self:
            if len( key ) != 2: 
                return dict.__getitem__( self, key )
            chrm, strand = key
            
            if chrm not in self.chrm_sizes or strand not in "+-":
                return dict.__getitem__( self, key )
            
            contig_size = self.chrm_sizes[ chrm ]
            self[ key ] = numpy.zeros( contig_size + 1)
        
        return dict.__getitem__( self, key )
    
    def __init__( self, chrm_sizes_fp, fps=[], strands=[ 'infer_from_fname', ]):
        self.chrm_sizes = ChrmSizes( chrm_sizes_fp )
        
        # initialize wiggle arrays
        """
        for chrm, size in self.chrm_sizes.iteritems():
            for strand in [ '+', '-' ]:
                # add one to make array 1-based
                self[ (chrm, strand) ] = numpy.zeros( size+1 )
        """
        
        if strands != ['infer_from_fname', ]:
            if len( fps ) != len( strands ):
                raise ValueError, "If the strand is specified, the number " \
                   + "of strand entries needs to be equal to the number of fps."
        else:
            strands = ["infer_from_fname"]*len( fps )
        
        for fp, strand in zip( fps, strands ):
            self.load_data_from_fp( fp, strand )
        
        self.intervals = None
        self.zero_intervals = None
        
        chrm_sizes_fp.seek( 0 )
        
        return

def main():
    #raise NotImplementedError
    chrm_sizes_fp = open( sys.argv[1] )
    fps = [ open( fname ) for fname in sys.argv[2:] ]
    wiggle = Wiggle( chrm_sizes_fp, fps )
    print wiggle
    #wiggle.calc_all_intervals()
    #print wiggle.write_wiggles( "tmp.plus.gff", "tmp.minus.gff" )

if __name__ == "__main__":    
    main()

