# Copyright (c) 2011-2012 Nathan Boley

from collections import namedtuple
GenomicInterval = namedtuple('GenomicInterval', ['chr', 'strand', 'start', 'stop'])

class GenomicIntervals( dict ):
    """Store a set of genomic intervals parsed for quick retrieval.
    
    This should probably always be subclassed. 

    """
    def __init__( self ):
        dict.__init__( self )
        self.is_frozen = False
        
        return
    
    def add( self, chrm, strnd, start, stop ):
        if self.is_frozen:
            raise TypeError( "The intervals list is frozen." )
        if chrm.startswith( "chr" ):
            chrm = chrm[3:]
        try:
            self[ (chrm, strnd ) ].append( (start, stop) )
        except KeyError:
            self[ (chrm, strnd ) ] = [ (start, stop), ]
        
        return
    
    def freeze( self ):
        """Convert the intervals into a sorted numpy array, for faster access.

        """
        for key in self:
            self[ key ] = numpy.array( self[ key ] )
            sorted_indices = self[ key ][:,0].argsort()
            self[ key ] = self[ key ][ sorted_indices ]
        
        self.is_frozen = True
        
        return
    
    def iter_indices_overlapping_a_region( self, chrm, strnd, start=0, stop=-1 ):
        """Iterate through the indicies of the regions which overlap a region.
        This can be optimized with both a start and stop sorted numpy.array.
        """
        if ( chrm, strnd ) not in self:
            return
        
        for index, (int_start, int_stop) in enumerate(self[ ( chrm, strnd ) ]):
            # if the internal region overlaps the provided region
            if int_stop >= start and ( stop == -1 or int_start <= stop ):
                yield index
        
        return
    
    def iter_intervals_overlapping_a_region( self, chrm, strnd, start=0, stop=-1 ):
        """ Iterate through the internal regions which 
        
        """
        for index in self.iter_indices_overlapping_a_region( chrm, strand, start, stop ):
            yield self[(chrm, strand)][index].tolist()
        
        return

    def get_overlapping_intervals( self, chrm, strand, start, stop ):
        return set( self.iter_intervals_overlapping_a_region( \
                chrm, strand, start, stop ) )
    
    def iter_joined_overlapping_intervals( self, chrm, strand, intervals, \
                                               interval_type_callback=None ):
        def does_overlap( ext_interval, int_interval ):
            if ext_interval[0] > int_interval[1] or ext_interval[1] < int_interval[0]:
                return False
            return True
        
        # sort the intervals by genomic start location
        intervals.sort()
        
        # materialize internal intervals
        internal_intervals = self[ ( chrm, strand ) ]
        assert isinstance( internal_intervals, numpy.ndarray )
        
        # sort internal_intervals
        internal_intervals = internal_intervals[ internal_intervals[:,0].argsort() ]
        
        # find allowable start pointer, and initialize to 1
        start_index = 0
        
        for ext_interval in intervals:
            overlapping_intervals = []
            for internal_interval in internal_intervals[start_index:]:
                if does_overlap( ext_interval, internal_interval ):
                    if interval_type_callback == None:
                        overlapping_intervals.append( internal_interval )
                    else:
                        type_of_overlap = interval_type_callback( \
                            ext_interval, internal_interval )
                        overlapping_intervals.append( \
                            (internal_interval, type_of_overlap) )
                
                # there can be no intervals before the current internal_interval that 
                # overlap because all previous starts will be less than the current 
                # start, due to sort order
                if internal_interval[1] < ext_interval[0]:
                    start_index += 1
                    
                # there can be no other intervals that overlap in this condition
                # because all subsequenct starts will be greater, due to sort order
                if internal_interval[0] > ext_interval[1]: 
                    break
                
            yield ext_interval, overlapping_intervals
        
        return
