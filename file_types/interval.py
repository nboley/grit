import sys
import math
import operator
import bisect
import re
import copy
from random import random
from itertools import izip
from operator import itemgetter

#from statistical_functions import *

verbose = False
output = sys.stdout

class ParseError(Exception):
    """Indicates that a file format is not to expectation.

    """
    def __init__(self, value="Parsing Error"):
        self.parameter = value
    def __str__(self):
        return self.parameter

def prefix_filename(filename, prefix):
    """Add a prefix to a filename. 

    i.e. filename = ./test_data/anns.bed, prefix = complement
    returns ./test_data/complement_anns.bed
    """
    ## TODO I should be using the sys path handling machinery here
    ## XXX This is probably broken on windows
    import re
    filename_portions = re.split("/", filename)
    filename_portions[-1] = prefix + filename_portions[-1]
    return '/'.join( filename_portions )

##################### Python2.3 Compatability Code #############################
try: set
except NameError:
    from sets import Set as set
##################### END Python2.3 Compatability Code ##########################

class interval(tuple):
    """Store a closed interval of ints. 
    
    The interval is closed at both ends. ie, interval(1,3) is {1,2,3}

    This is designed to be used to indicate genomic regions, thus the int. To
    this end, there are several functions that facilitate comparison.

    intersection: returns the intersection of self and another interval
    union: returns the union of self and another interval
    
    overlap: returns the size ( in bp's ) of an interval and a second interval
    does_overlap: returns whether or not an overlap overlaps

    """
    # this prevents each object from creating a __dict__ ( just a small optimization )
    __slots__ = []

    # creates the class from the object
    def __new__(self, start, end):
        if end < 0 or start < 0 or  not isinstance(start, (int, long)) != int or not isinstance(end, (int, long)):
            raise ValueError, 'All values must be non-negative (integers|longs), not (%i, %i), (%s, %s)' % ( start, end, type(start), type(end) )
        if end < start:
            raise ValueError, 'The end value (%d) must be greater than the start value (%d) ' % ( end, start )
        self = tuple.__new__(interval, (start, end))
        return self

    start = property(lambda self: self[0])
    end = property(lambda self: self[1])
    size = property(lambda self: self[1] - self[0]+1 )

    def intersection( self, otherInterval ):
        """Return the intersection of self and otherInt as a tuple of intervals.
           If the intersecton is empty, returns an empty tuple.
        """
        assert isinstance(otherInterval, interval)
        start = max(self.start, otherInterval.start)
        end = min(self.end, otherInterval.end)
        if start <= end:
            return (interval(start, end),)
        else:
            return ()

    def union( self, otherInterval ):
        """Return the union of self and otherInterval as a tuple of intervals.
        """
        assert isinstance(otherInterval, interval)
        if self.does_overlap \
           or otherInterval.start == self.end+1 \
           or self.start == otherInterval.end+1:
            start = min(self.start, otherInterval.start)
            end = max(self.end, otherInterval.end)
            return (interval(start, end),)
        else:
            return ( self, otherInterval )

    # calculate the length of the overlap
    def overlap(self, otherInterval):
        assert isinstance(otherInterval, interval)
        return max(min(self[1], otherInterval[1]) - max(self[0], otherInterval[0]) + 1, 0)

    def does_overlap(self, other):
        assert isinstance(other, interval)
        return self.overlap(other) > 0

    ############### Comparators ########################################################################################
    def __lt__(self, other):
        if isinstance(other, (int, long, float)):
            return self.end < other
        elif type(other) == interval:
            return self.end < other.start
        else: raise TypeError, "Can't compare an object of type %s to a slice object" % type(other)

    def __le__(self, other):
        if isinstance(other, (int, long, float)):
            return self.start <= other
        elif type(other) == interval:
            return self.end <= other.end
        else: raise TypeError, "Can't compare an object of type %s to a slice object" % type(other)

    def __eq__(self, other):
        # only allow the comparison of intervals or 'numbers'
        if not isinstance(other, (int, long, float, interval)): return False
        return self <= other and self >= other

    def __ne__(self, other):
        return not self == other

    def __ge__(self, other):
        if isinstance(other, (int, long, float)):
            return self.end >= other
        elif type(other) == interval:
            return self.start >= other.start
        else: raise TypeError, "Can't compare an object of type %s to a slice object" % type(other)

    def __gt__(self, other):
        if isinstance(other, (int, long, float)):
            return self.start > other 
        elif type(other) == interval:
            return self.start > other.end
        else: raise TypeError, "Can't compare an object of type %s to a slice object" % type(other)

class region(list):
    def __init__( self, intervals_iter=(), values_iter=(), name=None, length=None ):
        # create a list to store the values
        self.values = []
        list.__init__(self.values, values_iter)
        list.__init__(self, intervals_iter)
        self.name = name
        self._length = length

    def _set_length(self, new_length):
        assert isinstance(new_length, int)
        if len(self) >0 and new_length <= self[-1].end:
            raise ValueError, 'All bp indices must be less than the region length.'
        self._length = new_length

    def _get_length(self):
        return self._length

    length = property(_get_length, _set_length)

    def _merge_intervals( self, start, end, new, new_value ):
        """Merge overlapping intervals. End is *not* inclusive ( like range )

        start, end refers to the start, end indexes of the overlaping intervals 
        ( in self ). new, new_value are the new interval and it's value, 
        respectively.

        We sum the overlaps. For instance, [0,2] v 2 with [1,3] v 1 will yield
        [0,0]v2, [1,2]v3, [3,3]v2. Then, we delete the old intervals and insert
        the new.

        This is not so tough if we note that regardless of the size of end-start
        we can only have 2 interval splits. Everything else is just
        additions or value changes. 
        """
        # store the new intervals and their values
        new_intervals = []
        new_values = []

        start_bp = max(new.start, self[start].start)
        end_bp = min(new.end, self[end-1].end)

        # special case identical intervals
        if end - start == 1 and self[start].start == new.start and \
                self[start].end == new.end:
            self.values[start] += new_value
            return

        old_intervals = list.__getitem__(self, slice(start,end))
        old_intervals.sort(key=itemgetter(0))
        old_interval_starts = [ old_interval.start for old_interval in old_intervals ]
        old_interval_ends = [ old_interval.end for old_interval in old_intervals ]
        for index in xrange(len(old_intervals)-1):
            try: 
                # add the overlapping interval
                new_intervals.append( interval(old_interval_starts[index], old_interval_ends[index]) )
                new_values.append( self.values[start+index] + new_value )
                # add the 'new' interval
                new_intervals.append( interval(old_interval_ends[index]+1, old_interval_starts[index+1]-1) )
                new_values.append( new_value )
            except:
                # FIXME
                pass

        # if appropriate, add the pre-interval
        if new.start < old_interval_starts[0]:
            new_intervals.insert(0, interval(new.start, old_interval_starts[0]-1))
            new_values.insert(0, new_value)
        elif old_interval_starts[0] < new.start:
            new_intervals.insert(0, interval(old_interval_starts[0], new.start-1))
            new_values.insert(0, self.values[start])

        # if appropriate, add the post-interval
        try:
            if new.end > old_interval_ends[-1]:
                new_intervals.append(interval(old_interval_ends[-1]+1, new.end))
                new_values.append(new_value)
            elif old_interval_ends[-1] > new.end:
                new_intervals.append( interval(new.end+1, old_interval_ends[-1]), )
                new_values.append(self.values[end-1])
        except:
            import pdb
            pdb.set_trace()

        # now, remove the old intervals
        del self[start:end]
        del self.values[start:end]
    
        # add the new intervals
        for item, value in reversed(zip(new_intervals, new_values)):
            self.insert( start, item )
            self.values.insert( start, value )

    def add(self, coverage_region, value=1 ):
        # if the item belongs at the end, insert it
        # this optimizes the common append case while ensuring the data
        # remains in the correct order
        if self.length != None and coverage_region.end >= self.length:
            raise ValueError, 'All bps must be less than the region length'
        elif len(self) == 0 or ( len(self) > 0 and coverage_region > self[-1] ):
            list.append(self, coverage_region)
            self.values.append(value)
        # if the item doesn't belong at the end make sure it's unique and
        # insert it into the appropriate location
        else:
            # first, check if the intervals intersects
            l_loc = bisect.bisect_left(self, coverage_region)            
            # if the new element intersects a previous element
            if self[l_loc].does_overlap(coverage_region):
                r_loc = bisect.bisect_right(self, coverage_region)
                self._merge_intervals( l_loc, r_loc, coverage_region, value )
            else: 
                list.insert(self, l_loc, coverage_region)
                list.insert(self.values, l_loc, value)

    def _numRegions(self):
        return len(self)
    numRegions = property( _numRegions )


    def get_subregion(self, start_bp, stop_bp, shift_to_zero=False):
        """Return a copy of the genomic subregion. 
        """
        assert start_bp >= 0 and stop_bp >= 0
        
        if shift_to_zero: shift_value = start_bp
        else: shift_value = 0

        # find the start interval via the bisection algorithm
        start_index = bisect.bisect_left(self, start_bp)
        stop_index = bisect.bisect_right(self, stop_bp)

        intervals = self[slice(start_index, stop_index)]        
        
        return type(self)(      intervals_iter=( interval(
                                    (max(item.start,start_bp)-shift_value), 
                                    (min(item.end,stop_bp)-shift_value)
                                ) for item in intervals.iter_feature_regions() ),
                                values_iter=iter(self.values[start_index:stop_index]), 
                                name=self.name,
                                length=stop_bp-start_bp
        )

    def get_overlapping_intervals( self, inter ):
        """Return a copy of the genomic subregion that overlaps interval. 
        """
        
        # find the start interval via the bisection algorithm
        start_index = bisect.bisect_left(self, inter.start)
        stop_index = bisect.bisect_right(self, inter.end)

        inter = self[slice(start_index, stop_index)]        

        if start_index == stop_index: 
            length = 0
        else:
            length = self[stop_index-1].end  - self[start_index].start + 1

        return type(self)(      intervals_iter=(  interval(item.start, item.end)
                                  for item in inter.iter_feature_regions() ),
                                values_iter=iter(self.values[start_index:stop_index]), 
                                name=self.name,
                                length=length
        )
    
    def bp_score_array(self):
        """Return a numpy array of all scores. 
        
        If numpy is not available, return a list of scores. 
        """
        # FIXME optimize this
        try: 
            import numpy
            return numpy.array(list(self.iter_features()))
        except ImportError:
            return list(self.iter_features())
        
    def mean(self):
        """Return the mean of all values in self"""
        rv = 0
        for item, value in self.iter_intervals_and_values():
            rv += item.size*value
        return float(rv)/self.length

    def var(self):
        """Return the centered variance of all values in self"""
        mean = self.mean()
        ret_val = 0
        total_size = 0
        for item, value in self.iter_intervals_and_values():
            total_size += item.size
            ret_val += item.size*(value-mean)**2
        ret_val += (self.length - total_size)*(mean**2)
        return float(ret_val)/(self.length-1)

    def regionCorrelation(self, otherRegion):
        """Return the pearson correlation of self and other, by basepair"""
        # special case for empty intervals        
        if len(self) == 0 or len(otherRegion) == 0:
            print ValueError, 'Can\'t take the correlation of a 0 length vector'
            return 0.0
            raise ValueError, 'Can\'t take the correlation of a 0 length vector'

        # calculate the sum(x_i*y_i) term
        overlapProduct = 0
        other_index = 0

        """         
        for iv, value in self.iter_intervals_and_values():
            print iv, value
            # remove the other intervals that can never overlap (given the sort)
            while other_index < len(otherRegion) and \
                  otherRegion[other_index].end < iv.start:
                other_index += 1
            # if other index is longer than the length, there can be no matches
            if other_index >= len(otherRegion):
                break
            
            # calculate the length of the overlap
            overlap = iv.overlap( otherRegion[other_index] )
            print iv, value, otherRegion[other_index], otherRegion.values[other_index], overlap

            # if the two regions overlap, they contribute to the sum
            if overlap > 0:
                overlapProduct += overlap*( value*otherRegion.values[other_index] )
                # if so, move onto the next region in self
                continue
            # if not, we know that other is *not* << self and there is not
            # any overlap, so we continue tot he next self anyways.
        """
        d1 = self.bp_score_array()
        d2 = otherRegion.bp_score_array()
        try:
            import numpy
            overlapProduct = float((d1*d2).sum())
        except ImportError:
            overlapProduct = sum( [float(i1)*i2 for i1, i2 in izip(d1,d2)] )
        p_corr = overlapProduct - self.length*self.mean()*otherRegion.mean()
        p_corr = p_corr/((self.length-1)*math.sqrt(self.var())*math.sqrt(otherRegion.var()))            

        return p_corr

    def iter_intervals_and_values(self):
        for index in xrange(len(self)):
            yield self[index], self.values[index]

    def iter_features(self):
        """Iterates through every feature in the interval.
        
        For instance, if the region was of length 10 and had a single
        feature interval [2,8] this would return 0 0 1 1 1 1 1 1 1 0.
        """

        # deal with the start points
        for loop in xrange(self[0][0]):
            yield 0

        # deal with the internal features
        for loop in xrange(len(self)-1):
            value = self.values[loop]
            for iloop in xrange(self[loop].size):
                yield value
            for iloop in xrange(self[loop+1][0] - self[loop][1] - 1):
                yield 0

        # deal with the last feature
        value = self.values[-1]
        for loop in xrange(self[-1].size):
            yield value

        # deal with the non-features until the end
        for loop in xrange(self.length - self[-1][1]):
            yield 0

        return

    def featuresLength(self):
        # BAD python 2.3 compat change
        return sum( [ feature.size for feature in self.iter_feature_regions() ] )

    def iter_feature_regions(self):
        """Iterates through every feature region in the interval.
        
        For instance, if the region was of length 10 and had a single
        feature interval [2,8] this would return interval(2, 8).
        """

        for feature_interval in list.__iter__(self):
            yield feature_interval

        return

    def __iter__(self):
        """Make sure that we arent using iter incorrectly.

        If you want the list iterator call list.__iter__(self).
        """
        raise NotImplementedError, \
            "No defined iterator - call the particular fn directly ( ie iter_feature_regions, etc. )."

    def split(self, split_points):
        """Splits this region into a regions object.
        
        """
        rv = regions()
        
        # deal with the trivial case of no splits
        if len(split_points) == 0:
            rv[self.name] = self
            return rv

        # Determine the prefix of the name
        if self.name == None:
            name = ""
        else:
            name = self.name + "_"

        # extend the split points list to include the boundaries
        # we make a copy because, if split points *is* a list, the list 
        # init will just return a shallow copy.
        split_points = list(copy.copy(split_points))
        split_points.insert(0, 0)
        split_points.append(self.length)

        for index in xrange(0, len(split_points)-1):
            # get the subregion, shifting it to zero
            rv[name + 'split_%i' % (index+1)] = \
                self.get_subregion(split_points[index], split_points[index+1], True)
                
            # name the subregion appropriately
            rv[name + 'split_%i' % (index+1)].name = name + 'split_%i' % (index+1)

        # check to make sure nothing funny happens ( this catches an old bug
        # and can probably be removed, but -O is safer )
        assert self.length == sum( item.length for item in rv.values() )

        return rv

    def __getslice__(self, start, stop):
        return self.__getitem__(slice(start, stop, None))

    def __getitem__(self, key):
        # if the key is an integer, treat it like a normal list
        if type(key) == int:
            return list.__getitem__(self, key)        
        
        # FIXME optimize this
        # otherwise, the key better be a slice
        assert type(key) == slice
        return type(self)( iter(list.__getitem__(self, key)), 
                           iter(list.__getitem__(self.values, key)), 
                           self.name )
    
    def build_cdf(self):
        """Return a copy of this region in cumulative space. 
        
        1,1,1,2,2,2 becomes 1,2,3,5,7,9
        """
        new_region = cum_region( (), (), self.name, self.length )
        curr_value = 0
        for interval, value in self.iter_intervals_and_values():
            new_region.add( interval, (curr_value, curr_value+value*interval.size) )
            curr_value += value*interval.size
        return new_region

    def plot(self, split_points=[]):
        import pylab
        points  = [ entry for entry in self.iter_split_points(0, self.length) ]
        points.append(self.length)
        values = [ self.value(entry) for entry in points ]
        
        pylab.plot(points, values)
        for entry in split_points:
            pylab.axvline(entry)
        pylab.show()

    def iter_split_points(self, minValue, maxValue):
        iterator = list.__iter__(self)
        interval = iterator.next()
        while interval.end < minValue:  
            interval = iterator.next()
        while interval.start-1 < maxValue:  
            yield interval.start-1
            if interval.size > 1 and interval.end < maxValue:
                yield interval.end
            interval = iterator.next()       
        return

    def value(self, bp_index):
        """Return the value at a specific basepair.

        """
        # make sure that we're only trying to index by basepair
        assert isinstance(bp_index, (int, long))
        r_index = bisect.bisect_right(self, bp_index)

        ## if the l_index is 0, we havnt seen any non 0 values yet
        if r_index == 0:
            return 0
        ## else, find the previous interval, noting l_index > 0
        else: 
            prev_interval_start = self[r_index-1].end
        
        ## if we are in a no region area, return 0
        if prev_interval_start <= bp_index:
            return 0
        ## if we are in a region area, return the value at that point
        else:
            return self.values[r_index-1]

class cum_region(region):
    """Store region data as a frozen cumulative sum.

    The 'value' at a bp is the sum of the values *up to and including* that bp. 
    Also, the object is frozen, making way for future optimization.
    """

    def append( self, *args ):
        raise AttributeError, "'cum_region' object has no attribute 'append'"

    def __setitem__( self, *args ):
        raise TypeError, "'cum_region' object does not support item assignment"    

    def value(self, bp_index):
        """Return the value at a specific basepair.

        """
        # make sure that we're only trying to index by basepair
        assert isinstance(bp_index, (int, long))
        r_index = bisect.bisect_right(self, bp_index)

        ## if the l_index is 0, we havnt seen any non 0 values yet
        if r_index == 0:
            return 0
        ## else, find the previous interval, noting l_index > 0
        else: 
            prev_interval_start = self[r_index-1].end
        
        ## if we are in a no region area, return the cumsum at the next point
        if prev_interval_start <= bp_index:
            return self.values[r_index-1][1]
        ## if we are in a region area, return the cumsum at the next point
        ## plus 1 for every bp inside the region we are in. Note that we add one
        ## to this because the intervals are inclusive
        else:
            per_bp_factor = (self.values[r_index-1][1] - self.values[r_index-1][0])/self[r_index-1].size
            return self.values[r_index-1][0] + per_bp_factor*(bp_index - self[r_index-1].start + 1)


class binary_region(region):
    """A region object optimized to work with binary regions.

    """
    def add(self, interval):
        region.add( self, interval, 1 )

    def __init__(self, intervals_iter=(), name=None, length=None, values_iter=()):      
        region.__init__(self, intervals_iter, (), name, length)

    def _merge_intervals( self, start, end, new, new_value ):
        """Merge overlapping intervals down.

        This tells the region class to correctly merge intervals - rather than 
        add them just merge the intervals down.
        """
        
        # TODO fix the renaming hack - this is from a previous version
        l_loc, r_loc, newItem = start, end, new
        
        new_interval = interval( min(self[l_loc].start, newItem.start), max(self[r_loc-1].end, newItem.end) )
        if False:
            if verbose:
                print "WARNING! WARNING! WARNING!\n \tFeature interval %s intersects %i other FI(s) " \
                    % ( newItem, r_loc-l_loc )
                print "\tThe intersecting intervals will be merged."
                print "\t%s is being changed to %s" % ( self[l_loc], new_interval )
            for index in xrange(r_loc-1, l_loc, -1):
                if verbose:
                    print "\tThe interval %s is being removed." % str(self[index])
                del self[index]
            if verbose: print "WARNING! WARNING! WARNING!\n"
        self[l_loc] = new_interval

    # calculates the overlap between this and another set of intervals track
    def overlap(self, otherAnnTrack):
        # special case for empty intervals        
        if len(self) == 0 or len(otherAnnTrack) == 0: return 0

        totalOverlap = 0

        currentMatches = []
        
        thisIter = self.iter_feature_regions()
        nextMatch = thisIter.next()
        
        for iv in otherAnnTrack.iter_feature_regions():
            # first, get rid of any current matches that dont match
            # because of the ordering, these should start not matching at the begining
            for item in currentMatches:
                if item.end >= iv.start: break
                else: del currentMatches[0] 

            # next, add any new items to currentMatches
            while nextMatch != None and nextMatch.start <= iv.end:
                currentMatches.append(nextMatch)
                try: nextMatch = thisIter.next()
                except StopIteration: nextMatch = None

            # finally, calculate the overlap of every item in currentMatches
            for item in currentMatches:
                totalOverlap += item.overlap(iv)

        return totalOverlap

    def regionCorrelation(self, otherRegion):
        """ Calculate the correlation of this region with another region

        """
        
        if len(self) == 0 or len(otherRegion) == 0: return 0
        
        if not(self.length == otherRegion.length):  
            raise ValueError, "Can't calculate correlation of regions of unequal length: (%i, %i)" % (self.length,otherRegion.length)
        
        n = self.length
        overlap = self.overlap(otherRegion)
        bp1 = self.featuresLength()
        bp2 = otherRegion.featuresLength()
        corr = ((n*overlap)-(bp1*bp2))/(math.sqrt(n*bp1-math.pow(bp1,2)) * math.sqrt(n*bp2-math.pow(bp2,2)))

        return corr

    # calculates the region overlap between this and another set of intervals track
    def regionOverlap(self, otherRegion, min_coverage_percent=0.0):
        """Calculate the number of regions that overlap self with the given region

        We take self to be the covered region. ie, for 
        self    ===============================
        other   ---   ----   -----------    -         
        the region overlap is 1, whereas for  
        self    ===========  ==================
        other   ---   ----   -----------    -                 
        it is 2.

        The algorithm is nearly identical to they overlap fn - it's a merge join
        """

        # special case for empty intervals        
        if len(self) == 0: return 0

        totalOverlap = 0
        self_index = 0
        other_index = 0
         
        for iv in list.__iter__(self):
            # remove the other intervals that can never overlap (given the sort)
            while other_index < len(otherRegion) and \
                  otherRegion[other_index].end < iv.start:
                other_index += 1
            # if other index is longer than the length, there can be no matches
            if other_index >= len(otherRegion):
                break

            # check to see if there is an overlap
            overlap = iv.overlap( otherRegion[other_index] )
            if overlap >= 1 and overlap >= min_coverage_percent*iv.size:
                otherRegion[other_index], otherRegion[other_index].size, \
                overlap, min_coverage_percent*iv.size
                totalOverlap += 1
                # if so, move onto the next region in self
                continue
            # if not, we know that other is *not* << self and there is not
            # any overlap, so we continue tot he next self anyways.

        #if self.brute_regionOverlap( otherRegion ) != totalOverlap:
        #    print "WARNING!!! NOT EQUAL!!!!", self.brute_regionOverlap( otherRegion ), totalOverlap

        return totalOverlap

    def bruteOverlap(self, otherRegion):
        # mostly for testing my fancier overlap algorithm
        # this compares every possible combination of intervals
        totalOverlap = 0
        for other_iv in otherRegion:
            for self_iv in self:
                totalOverlap += other_iv.overlap(self_iv)
        return totalOverlap

    def brute_regionOverlap(self, otherRegion):
        # mostly for testing my fancier overlap algorithm
        # this compares every possible combination of regions
        totalOverlap = 0
        for self_iv in list.__iter__(self):
            for other_iv in list.__iter__(otherRegion):
                if other_iv.does_overlap(self_iv):
                    totalOverlap += 1
                    break
        return totalOverlap

    def iter_intervals_and_values(self):
        for index in xrange(len(self)):
            yield self[index], 1

    def iter_features(self):
        """Iterates through every feature in the interval.
        
        Overrides the iter_features for region class and applies to binary regions
        
        For instance, if the region was of length 10 and had a single
        feature interval [2,8] this would return 0 0 1 1 1 1 1 1 1 0.
        """

        # deal with the start points
        for loop in xrange(self[0][0]):
            yield 0

        # deal with the internal features
        for loop in xrange(len(self)-1):
            for iloop in xrange(self[loop][1] - self[loop][0]):
                yield 1
            for iloop in xrange(self[loop+1][0] - self[loop][1]):
                yield 0

        # deal with the last feature
        for loop in xrange(self[-1][1] - self[-1][0]):
            yield 1

        # deal with the non-features until the end
        for loop in xrange(self.length - self[-1][1]):
            yield 0

        return

class genomic_coverage_region( region ):
    def _length(self):
        try: return self._cached_length
        except AttributeError:
            # BAD python 2.3 compat change
            self._cached_length = sum( [ feature.size for feature in self.iter_feature_regions() ] )
        return self._cached_length
    length = property( _length )  

    # FIXME fix the cached length in set_attr

    def append(self, newItem):
        # clear the cached length 
        if hasattr(self, "_cached_length"):
            del self._cached_length
        # call region's append method 
        region.append(self, newItem)

    def intersection(self, inter):
        """Return the intersection of inter and self, *in genomic coordinates*
        
        inter needs to be an interval.
        """

        # first, find the lower and upper bounds of intersection
        l_loc = bisect.bisect_left(self, inter.start)
        r_loc = bisect.bisect_right(self, inter.end)        
        # if there is no overlap, return None
        if l_loc == r_loc: return None

        # make sure that self has *some* regions
        assert len(self) > 0
        # make sure the interval is before the first basepair
        assert self[0].start <= inter.end
        # makre sure the interval is after the last basepair
        assert self[-1].end >= inter.start

        # otherwise, calculate the overlap
        # BAD python 2.3 compat change
        bp_overlap = sum( [item.overlap(inter) for item in list.__iter__(self)] )
        # next, note the interval start is the count of bp's <= self[l_loc] 
        if l_loc > 0:
            new_start = (self[0:l_loc]).length
        else:
            new_start = 0
        # account for the potential overlap
        if inter.start > self[l_loc].start:
            new_start += ( inter.start - self[l_loc].start )
        
        # make sure that the new start is positive - ( explicitly catch ian's bug )
        assert new_start >= 0

        return interval( new_start, new_start+bp_overlap-1 )

class regions(dict):
    """A region container object. 
    
    This is just a dictionary of region's whereby each key is the 
    region's track name. Also contains a method to store an overlap stat.
    """

    #!!! BUG - FIX THE LENGTH CACHING MECHANISM
    # It currently assumes that I never use __setitem__

    def _totalLength(self):
        """Calculate the total length and cache it for later use.

        """
        try: return self._cached_totalLength
        except AttributeError:
            self._cached_totalLength = float(sum([ value.length for value in self.values() ]))
            return self._cached_totalLength

    totalLength = property( _totalLength )

    def regionFraction(self):
        """Returns a dict of the relative region lengths.

        For instance, if there are two regions of names R1 and R2 and they are both
        the same length then this will return { 'R1': 0.5, 'R2': 0.5 }

        """
        # BAD python2.3 compat change
        return dict([ ( key, float(value.length)/self.totalLength ) for key, value in self.iteritems() ])

    def extend(self, other):
        for key in other.keys():
            # make sure that we wont overwrite an existing region
            if self.has_key(key):
                raise ValueError, 'Can not extend %r with %r: they both contain a region named %s' % (self, other, key)
            self[key] = other[key]

        # delete the cached length
        try: 
            del self._cached_totalLength
        except AttributeError:
            pass

    def writeBedFile(self, bed_f, lengths_f):
        keys = self.keys()
        keys.sort()

        for key in keys:
            for (start, end), value in self[key].iter_intervals_and_values():
                bed_f.write("%s\t%i\t%i\n" % (key, start, end))

            lengths_f.write("%s\t0\t%i\n" % (key, self[key].length) )

    def writeGff3File(self, gff3_f):
        keys = self.keys()
        keys.sort()

        for key in keys:
            for i, ((start, end), value) in enumerate(self[key].iter_intervals_and_values()):
                if isinstance(self[key], binary_region):
                    gff3_f.write("\t".join((str(key), "merged_feature", "merged_feature", str(start), str(end), ".", ".", ".", "ID=" + str(i))) + "\n")
                else:  # 
                    gff3_f.write("\t".join((str(key), "merged_feature", "merged_feature", str(start), str(end), str(value), ".", ".", "ID=" + str(i))) + "\n")


class lazy_wiggle_file(object):
    """A memory safe wiggle parser.
    
    We region subselects that are done on disk to increase access efficiency
    when the entire file is too big to fit in memory.    
    """

    def md5file(self, fh):
        """Return the hex digest of a file without loading it all into memory
        
        This caches the data structure and checks that it hasnt changed."""
        import md5
        digest = md5.new()
        while 1:
            buf = fh.read(4096)
            if buf == "":
                break
            digest.update(buf)
        self.f.seek(0)
        return digest.hexdigest()

    def _build_annotation_points(self):
        """Scan through the file and find and store all the declaration lines.

        """
        def parse_dec_line( line, f_pos ):
            values = dict(re.findall('(\S+)=(\S+)', line))
            try: 
                values['start'] = int(values['start'])
            except KeyError:
                values['start'] = 0

            try: values['span'] = int(values['span'])
            except KeyError: values['span'] = 1
            
            try: values['step'] = int(values['step'])
            except KeyError: pass
            
            if line.startswith("variableStep"):
                values['step_type'] = 'Variable'
            elif line.startswith("fixedStep"):
                values['step_type'] = 'Fixed'
            else:
                pass
                assert False

            values['f_pos'] = f_pos

            # store internal points and their file positions. For long 
            # regions, this can significantly speed up the subregion seelcts
            # Each item in the list is a tuple ( start_bp, file position )
            values['internal_indexes'] = []
            
            return values

        def add_dec_line( values ):
            # store the declaration line information
            if not self.declarations.has_key(values['chrom']):
                self.declarations[values['chrom']] = []
            insert_index = bisect.bisect_left( self.declarations[values['chrom']], values['start'])
            self.declarations[values['chrom']].insert( insert_index, values )

        self.f.seek(0)
        # we use the readline method to preserve compatibility with tell
        # we might want to optimize this with a bigger read
        values = None
        f_pos = 0
        line_cntr = 0
        line = self.f.readline()
        while line != '':
            # ignore track lines ( and their multi line spans )
            if line.startswith("track"):
                f_pos = self.f.tell()
                line = self.f.readline()
                while line.endswith("\\") and line != "":
                    f_pos = self.f.tell()
                    line = self.f.readline()
                continue

            # parse declaration lines
            if line.startswith('variableStep') \
               or line.startswith('fixedStep'):
                # reset the line counter
                line_cntr = 0
                # determine the end of the previous set
                if values != None:
                    add_dec_line( values )
                line_data = parse_dec_line( line, f_pos )
                values = line_data
                if values['start'] != None:
                    values['end'] = values['start'] + values['span']
            # look for the end
            elif values != None:
                line_cntr += 1
                if values['step_type'] == 'Variable':
                    curr_bp = int(re.split("\s+", line)[0])
                    if values['start'] == None:
                        values['start'] = curr_bp
                    values['end'] = curr_bp + values['span']
                else:
                    assert values['step_type'] == 'Fixed'
                    values['end'] += values['step']
                # if this is a long region, add an index point in
                if line_cntr >= self.max_lines:
                    values['internal_indexes'].append( (values['end'], f_pos) )
                    line_cntr = 0
            f_pos = self.f.tell()
            line = self.f.readline()
        if values != None:
            add_dec_line( values )

    
    def get_subregion(self, chr_name, start_bp, stop_bp):
        """Get a subregion.
        
        This is a memory efficient way to get a region without having to parse 
        the entire file.
        """
        step_type = None
        
        all_values = [ values for values in self.declarations[chr_name] \
                       if values['start'] <= stop_bp and values['end'] >= start_bp ]
        sub_region_length = max(values['end'] for values in all_values)
        sub_region = region((), (), chr_name, max(values['end'] for values in all_values) )
        for values in all_values:
            self.f.seek(values['f_pos'])
            # this should be a declaration line
            line = self.f.readline()
            assert line.startswith('variableStep') or line.startswith('fixedStep')

            # now, look at the entires in values['internal_indexes'] and fast 
            # forward to a sensible location. Namely, we want the first bp loc
            # st the bp is strictly less than start_bp. Then, it cannot be 
            # included in the correct interval.
            #TODO use a bisection algorithm
            search_loc = None
            for bp_index, f_loc in values['internal_indexes']:
                if bp_index < stop_bp:
                    search_loc = f_loc
                else:
                    break
            if search_loc != None:
                self.f.seek(search_loc)
            
            if values['step_type'] == 'Fixed':
                pos = int(values['start']) - values['step']
            
            # start parsing the values
            line = self.f.readline()
            while line != "" \
                and not line.startswith('variableStep') \
                and not line.startswith('fixedStep'):
                line = line.strip()
                if values['step_type'] == 'Variable':
                    tmp_pos, value = re.split("\s+", line)
                    value = float(value)
                    pos = int(tmp_pos)
                elif values['step_type'] == 'Fixed':
                    pos += values['step']
                    value = float(line.strip())
                
                sub_region.add( interval(int(pos), int(pos)+values['span']-1), value )
                
                # if we've passed the region, there is nothing left to do
                if pos > stop_bp:
                    break

                line = self.f.readline()
        return sub_region.get_subregion(start_bp, stop_bp)
    
    def __init__(self, wig_f):
        import cPickle, os
        self.max_lines = 10000
        self.declarations = sc_dict()
        # open and store the file reference
        if isinstance(wig_f, (str, unicode)):
            self.opened_f = True
            self.f = open(wig_f)
        else:
            self.opened_f = False
            self.f = wig_f
        
        # find the file's checksum
        self.md5sum = self.md5file( self.f )
        # see if we have a cached version of declarations lying around
        # if so, unpickle it and do nothing else

        # try to open the corrrectly named file
        fname = ".cached_annotation_pnts" + self.f.name      
        if os.path.exists(fname) and os.path.isfile(fname):
            # open the file
            f = open( fname, 'rb' )
            # check to make sure that the checksum matches
            self.declarations = cPickle.load(f)
            if self.md5sum != self.declarations.md5sum:
                self.declarations = sc_dict()
            f.close()
        
        if len(self.declarations) == 0:
            self._build_annotation_points()
            self.declarations.md5sum = self.md5sum
            f = open( fname, 'wb' )
            cPickle.dump( self.declarations, f )
            f.close()
    
    def __str__( self ):
        return str(self.declarations)

def parse_wiggle_file(wig_f, lengths_file, BINARY=False):
    """Parse a wiggle file into a regions object.

    """
    lengths_file.seek(0)

    # first, build the lengths file region
    lengths = regions()
    for line_num, line in enumerate(lengths_file):
        values = re.split('\s+', line.strip())
        if len(values) < 3:
            raise ParseError("Error parsing line #%i \"%s\" in the file \"%s\"" % (line_num+1, line.strip(), lengths_file.name)) 
        chName, start, end = values[0:3]

        feature = interval(int(start), int(end))
        # if the lengths object doesnt have a region of this name, create it
        if not lengths.has_key(chName):
            lengths[chName] = genomic_coverage_region(name=chName)
        lengths[chName].append(feature)

    # next, initialize the chromosomes data structure
    chromosomes = regions()
    for chName in lengths.keys():
        length = lengths[chName].length
        if BINARY:
            chromosomes[chName] = binary_region((), chName, length)
        else:
            chromosomes[chName] = region((), (), chName, length)
   
    # loop through each line
    curr_chr = None
    curr_span = None
    last_pos = 1
    line = wig_f.readline()
    while line != '':
        # skip comment lines
        if line.startswith('#'):
            line = wig_f.readline()
            continue
        elif line.startswith('track'):
            line = line.strip()
            while line.endswith('\\'):
                line = line[:-1] + wig_f.readline().strip()
            line = wig_f.readline()
            continue
        elif line.startswith('fixedStep'):
            raise NotImplementedError, \
                'Only VariableStep wiggle annotations are permitted.'
        # we are at a declaration line
        elif line.startswith('variableStep'):
            line = line.strip()
            while line.endswith('\\'):
                line = line[:-1] + wig_f.readline().strip()
                
            values = dict(re.findall('(\S+)=(\S+)', line))
            curr_chr = values['chrom']

            # reset the previous pos counter
            try:
                last_pos = chromosomes[curr_chr][-1].end+1
            # if the list is empty, se tthe counter to 0
            except IndexError:
                last_pos = 0
            # if chromosomes doesnt have the key, curr_chr ignore it
            except KeyError:
                pass

            try: curr_span = int(values['span'])
            except KeyError:
                curr_span = None
        # else, it's a data line
        else:
            line = line.strip()
            pos, value = map(float, re.split("\s+", line))
            # the region objects view empty spaces as 0's
            if value != 0:
                # if the region is binary, non-zero regions have a value of 1
                if BINARY:
                    value = 1
                # XXX? http://genome.cse.ucsc.edu/goldenPath/help/wiggle.html
                # The standard is not so clear on what to do if no span is defined.
                # However, from the example, it seems to default to 1
                if curr_span == None:
                    chromosomes[curr_chr].add( interval(int(pos), int(pos)), value )
                    #chromosomes[curr_chr].add( interval(last_pos, int(pos)), value )
                else:
                    chromosomes[curr_chr].add( interval(int(pos), int(pos)+curr_span-1), value )
                last_pos = int(pos) + 1
        line = wig_f.readline()
    return chromosomes

def parse_bed_line_iterator( bed_f, lengths_file, group=None, is_binary=False ):
    import re
    
    # first, build the lengths file region
    lengths = regions()
    for line_num, line in enumerate(lengths_file):
        values = re.split('\s+', line.strip())
        if len(values) < 3:
            raise ParseError("Error parsing line #%i '%s' in the file '%s'"
                % (line_num+1, line.strip(), lengths_file.name)) 
        chName, start, end = values[0:3]

        feature = interval(int(start), int(end))
        # if the lengths object doesnt have a region of this name, create it
        if not lengths.has_key(chName):
            lengths[chName] = genomic_coverage_region(name=chName)
        lengths[chName].append(feature)
    
    # next, initialize the chromosomes data structure
    chromosomes = regions()
    for chName in lengths.keys():
        length = lengths[chName].length
        if is_binary:
            chromosomes[chName] = binary_region((), chName, length)
        else:
            chromosomes[chName] = region((), (), chName, length)
    for line in bed_f:
        if line.strip() == "": continue
        if line.startswith("#"): continue
        # parse out the region name, and interval and parse it
        values = re.split('\s+', line.strip())
        chName, start, end = values[0:3]
        feature = interval(int(start), int(end))
        
        # if we are filtering by group
        if group is not None:
            try:
                name = values[3]
                if name != group:
                    continue  # Feature not in group: skip it
            except IndexError:
                raise ParseError("Error parsing groups in \"%s\": no name field found" % bed_f.name)
        
        if is_binary:
            # FIXME add warning here for missing region?

            # see what parts - if any - of the feature intersects the genomic region
            # take only the intersecting pieces and shift the coordinates to be
            # measured in the genomic region space. If the key doesnt exist in
            # the lengths file, then ignore the region.
            if not lengths.has_key( chName ):
                continue
            shifted_feature = lengths[chName].intersection(feature) 
            # if there is no intersection, there is nothing else to do
            if shifted_feature == None: continue
            else: 
                try:
                    chromosomes[chName].add(shifted_feature)
                # This assumes that, if the region name is not present in the lengths file,
                # then we want to ignore it
                except KeyError:
                    pass
        else:
            value = float(values[4])  # Score field of BED file

            # see what parts - if any - of the feature intersects the genomic region
            # take only the intersecting pieces and shift the coordinates to be
            # measured in the genomic region space. If the key doesnt exist in
            # the lengths file, then ignore the region.
            if not lengths.has_key( chName ):
                continue
            shifted_feature = lengths[chName].intersection(feature) 
            # if there is no intersection, there is nothing else to do
            if shifted_feature == None: continue
            else: 
                try:
                    chromosomes[chName].add(shifted_feature, value)
                # This assumes that, if the region name is not present in the lengths file,
                # then we want to ignore it
                except KeyError:
                    pass
        
        continue    
    
    lengths_file.seek(0)
    return chromosomes

def parse_gff3_file( gff3_f, lengths_file, group=None ):
    def iter_gffslines_as_bedlines():
        for line in gff3_f:
            if line.startswith( "#" ): continue
            this_line = "\t".join( item for loop, item in \
                                   enumerate(re.split('\s+',line.strip())) \
                                   if loop in (0,3,4) ) \
                      + "\n"
            yield this_line
        
    # BUG TODO FIXME Allow gff3 to parse non-binary formats
    projected_bed = parse_bed_line_iterator(\
        iter_gffslines_as_bedlines(), lengths_file, group=group, is_binary=True)

    return projected_bed  

def parse_bed_file(bed_f, lengths_file, group=None, force_binary=False):
    """Parse a bed and lengths file into a regions object.

       If group is specified, only BED lines whose "name" line matches
       group are added to the regions object
    """

    import re

    # make sure that we are at the beginning of the files
    bed_f.seek(0)
    lengths_file.seek(0)

    # check to see if the features are binary or not
    values = re.split('\s+', bed_f.readline().strip())
    if len(values) < 3:
        print >> sys.stderr, "Offending values: '%s'\n" % '\t'.join( values )
        raise ParseError("Error parsing \"%s\": this is not a valid BED file" \
            % bed_f.name) 
    elif len(values) == 3 or force_binary: 
        BINARY = True
    else: 
        BINARY = False
    bed_f.seek(0)
    
    chromosomes = parse_bed_line_iterator( bed_f, lengths_file, group, BINARY )
    
    return chromosomes
