# Copyright (c) 2011-2012 Nathan Boley

import re
import numpy

from itertools import product

from collections import namedtuple
GenomicInterval = namedtuple('GenomicInterval', ['chr', 'strand', 'start', 'stop'])

class GenomicIntervals( dict ):
    """Store a set of genomic intervals parsed for quick retrieval.
    
    This should probably always be subclassed. 

    """
    def __init__( self ):
        dict.__init__( self )
        self.is_frozen = False
    
    def add( self, chrm, strnd, start, stop ):
        if self.is_frozen:
            raise TypeError( "The intervals list is frozen. " )
        if chrm.startswith( "chr" ):
            chrm = chrm[3:]
        try:
            self[ (chrm, strnd ) ].append( (start, stop) )
        except KeyError:
            self[ (chrm, strnd ) ] = [ (start, stop), ]
    
    def freeze( self ):
        """Convert the intervals into a sorted numpy array, for faster access.

        """
        for key in self:
            self[ key ] = numpy.array( self[ key ] )
            sorted_indices = self[ key ][:,0].argsort()
            self[ key ] = self[ key ][ sorted_indices ]
        
        self.is_frozen = True

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
    
class Junctions( GenomicIntervals ):
    def iter_connected_exons( self, chrm, strand, start, stop, exon_bndrys, \
                                  include_scores=False ):
        """Find the exon connections that this junction list implies.
        if include_scores then also return connected_exons scores
        """
        jn_indices = self.iter_indices_overlapping_a_region( \
            chrm, strand, start, stop )
        paired_exons = []
        paired_exon_scores = []
        for index in jn_indices:
            jn_start, jn_stop = self[(chrm, strand)][index]
            jn_score = self._scores[(chrm, strand)][index]
            # we use +/- 1, because the junctions are intron start/ends, inclusive
            start_exons = [ i for i, bndries in enumerate( exon_bndrys ) 
                            if bndries != None and jn_start - 1 == bndries[1] ]
            stop_exons = [ i for i, bndries in enumerate( exon_bndrys ) 
                            if bndries != None and jn_stop + 1 == bndries[0] ]
            models = [ (i,j) for (i,j) in  product( start_exons, stop_exons ) if i <= j ]
            paired_exons.extend( models )
            paired_exon_scores.extend( [jn_score]*len(models) )
        
        if include_scores:
            return paired_exons, paired_exon_scores
        return paired_exons
    
    def __init__( self ):
        GenomicIntervals.__init__( self )
        self._scores = {}
        
    def add( self, chrm, strnd, start, stop, score=0 ):
        GenomicIntervals.add( self, chrm, strnd, start, stop )
        try:
            self._scores[ (chrm, strnd) ].append( score )
        except KeyError:
            self._scores[ (chrm, strnd) ] = [ score, ]
    
    def freeze( self ):
        """Convert the intervals into a sorted numpy array, for faster access.
        Also order scores to maintain indices between intervals
        """
        for key in self:
            # order start, stops by start coordinate
            self[ key ] = numpy.array( self[ key ] )
            sorted_indices = self[ key ][:,0].argsort()
            self[ key ] = self[ key ][ sorted_indices ]
            
            # use same sorted_indices to sort _scores in the same order
            self._scores[ key ] = numpy.array( self._scores[ key ] )
            self._scores[ key ] = self._scores[ key ] [sorted_indices ]
        
        self.is_frozen = True

def parse_junctions_file( junctions_fp ):
    jns = Junctions()
    for line in junctions_fp:
        gene_name, trans_name, feature_type, region, score = parse_gff_line( line, True )            
        # if this is a valid gff/gtf line (i.e. not a header)
        if feature_type != None:
            jns.add( region.chr, region.strand, region.start, region.stop, score )            
    
    jns.freeze()
    return jns

class TSSs( GenomicIntervals ):    
    def get_start_exons( self, gene_bndry, exon_bndrys ):
        tss_s = sorted( self.get_overlapping_intervals( gene_bndry ) )
        start_exons = set()
        for p_start, p_stop in tss_s:
            for i, (e_start, e_stop) in enumerate(exon_bndrys):
                # they overlap if the tss is bot not strictly before
                # the exon and not stricly after.
                if not( p_stop < e_start or p_start > e_stop ):
                    start_exons.add( i )
        
        return sorted( start_exons )

def parse_tss_file( tss_fp ):
    tss_s = TSSs()
    for line in tss_fp:
        gene_name, trans_name, feature_type, region = parse_gff_line( line )
        # if this is a valid gff/gtf line (i.e. not a header)
        if feature_type != None:
            tss_s.add( region.chr, region.strand, region.start, region.stop  )
        
    tss_s.freeze()
    return tss_s

parse_tes_file = parse_tss_file

def parse_gff_line( line, include_score=False ):
    data = line.split()
    
    if len( data ) < 8:
        return None, None, None, None
    # the type of element - ie exon, transcript, etc.
    feature_type = data[2]
    
    # check that this is a valid gff/gtf line else return None
    # check that start and stop are valid integers
    try:
        start = int(data[3])
        stop = int(data[4])
    except ValueError:
        return None, None, None, None
    # check for valid strand
    strand = data[6]
    if strand not in [ '+', '-', '.' ]:
        return None, None, None, None
    
    if include_score:
        try: score = int( data[5] )
        except ValueError: score = 0
    
    # parse the meta data, and grab the gene name
    meta_data = dict( zip( data[8::2], data[9::2] ) )
    # remove ';' from value fields if nessecary
    for key, val in meta_data.iteritems():
        if val.endswith( ';' ):
            meta_data[ key ] = val[:-1]

    # get gene and transcript name if parsing a gtf line 
    # else it is a gff line and does not have gene or trans names
    gene_name = None
    trans_name = None
    if ('gene_id' in meta_data) and ('transcript_id' in meta_data):
        gene_name = meta_data[ 'gene_id' ]
        trans_name = meta_data[ 'transcript_id' ]
    
        if gene_name.startswith('"'):
            gene_name = gene_name[1:-1]
        if trans_name.startswith('"'):
            trans_name = trans_name[1:-1]

    # fix the chromosome if necessary
    chrm = data[0]
    if chrm.startswith("chr"):
        chrm = data[0][3:]
        
    if include_score:
        return gene_name, trans_name, feature_type, GenomicInterval( \
            chrm, strand, start, stop ), score
    return gene_name, trans_name, feature_type, GenomicInterval( \
        chrm, strand, start, stop )

class GeneBoundaries( dict ):
    def _validate_exons( self ):
        """Test all of the genes for consistency. 
        If a gene is bad, move it to _invalid_genes
        
        """
        self._invalid_genes = {}
        # fix exons that overlap
        for gene_name in self.keys():
            # get the set of exons for this gene
            exons = list( self[ gene_name ] )
            
            try:
                # validate the exon inputs
                if not all( exon.chr == exons[0].chr for exon in exons ):
                    raise ValueError, \
                        "There appears to be a gene with exons from different chromosomes."
                if not all( exon.strand == exons[0].strand for exon in exons ):
                    raise ValueError, \
                        "There appears to be a gene with exons from different strands."
                strand = exons[0].strand
                chr = exons[0].chr
            except ValueError, inst:
                print( "Error in %s: %s" % ( gene_name, str(inst) ) )
                self._invalid_genes[ gene_name ] = self[ gene_name ]
                del self[ gene_name ]
        
    
    def __init__( self, fp, fasta=None ):
        """Parse the gtf file at fp and return the genes and their unique exons.
        
        """
        self.filename = fp.name
        self.fasta = fasta

        # iterate through each line in the gtf, and organize exons
        # based upon the genes that they came from
        for line in fp:
            gene_name, trans_name, feature_type, region = parse_gff_line( line )
            # skip non-exon lines
            if feature_type != 'exon': continue
            
            try:
                self[gene_name].add( region )
            except KeyError:
                self[gene_name] = set()
                self[gene_name].add( region )
        
        # make sure that all of the exons make sense. ie, from the 
        # same chromosome, etc.
        self._validate_exons()
        
        # since we know that all of the remaining genes are on the same strand
        # and in the same chromosome, change the format so that the genes are
        # a chr, strand, and array of starts and stops
        for gene_name in self.keys():
            # get the set of exons for this gene
            exons = list( self[ gene_name ] )
            # build an array of the exon starts and stops to faciliate easier
            # inclusion tests
            self[ gene_name ] = \
                GeneModel( gene_name, exons[0].strand, exons[0].chr, \
                               ( ( exon.start, exon.stop ) for exon in exons ) )
        
        return

class GeneModel( object ):
    def exons_do_overlap( self, bndry1, bndry2 ):
        if isinstance( bndry1, ( list, tuple ) ):
            return( bndry1[1] >= bndry2[0] )
        else:
            bndry1 = self.bp_bndry( bndry1 )
            bndry2 = self.bp_bndry( bndry2 )
            return( bndry1[1] >= bndry2[0] )
    
    def bp_bndry( self, index ):
        """Get the basepair bounds for exon 'index'.
        """
        return self.exon_bndrys[index]
    
    def exon_index( self, start, stop ):
        """Get the exon index for ( start, stop ).
        """
        return self.exon_bndrys.index( (min(start, stop), max(start, stop)) )
    
    def _get_boundaries( self ):
        min_val = min( i[0] for i in self.exon_bndrys if i != None )
        max_val = max( i[1] for i in self.exon_bndrys if i != None )
        
        return GenomicInterval( \
            self.chromosome, self.strand, min_val, max_val )
    
    def _get_nonoverlapping_exon_boundaries( self ):
        boundaries = set()
        for start, stop in self.exon_bndrys:
            boundaries.add( start )
            # we add one so that the boundaries are always
            # region starts, noting that gtf's are closed-closed
            boundaries.add( stop+1 )
        
        return numpy.array( sorted( boundaries ) )
    
    def __init__( self, name, strand, chromosome, exon_bndrys ):
        self.name = name
        self.strand = strand
        self.chromosome = chromosome
        
        # if this is a dictioanry, preserve the exon indicies.
        if isinstance( exon_bndrys, dict ):
            self.exon_bndrys = [None]*( 1+ max( exon_bndrys.keys() ) )
            for key, val in exon_bndrys.iteritems():
                self.exon_bndrys[ key ] = val
            for index, val in enumerate( self.exon_bndrys[1:] ):
                if val == None:
                    self.exon_bndrys[ index + 1 ] = self.exon_bndrys[ index ] 
        else:
            self.exon_bndrys = sorted(exon_bndrys)
        
        self.exon_starts = \
            numpy.array(sorted(set( item[0] for item in self.exon_bndrys )))
        self.exon_stops = \
            numpy.array(sorted(set( item[1] for item in self.exon_bndrys )))
        self.exon_lens = \
            [ stop - start for i, (start, stop) in enumerate(self.exon_bndrys) ]
        
        self.boundaries = self._get_boundaries()
        self.nonoverlapping_exon_boundaries = \
            self._get_nonoverlapping_exon_boundaries()
        self.nonoverlapping_exon_lens = \
            [ stop - start for start, stop in \
                  zip(self.nonoverlapping_exon_boundaries[:-1], \
                          self.nonoverlapping_exon_boundaries[1:]) ]

        
        self._possible_junctions = None
        
        return
    
    def find_nonoverlapping_exons( self, exon ):
        start, stop = self.exon_bndrys[ exon ]
        start_nonoverlapping_exon = self.nonoverlapping_exon_boundaries.searchsorted( \
            start, side='right' ) - 1
        stop_nonoverlapping_exon = self.nonoverlapping_exon_boundaries.searchsorted( \
            stop + 1, side='right' ) - 2
        
        return range( start_nonoverlapping_exon, stop_nonoverlapping_exon + 1 )
    
    def __str__( self ):
        ret_strs = ["Gene: " + self.name]
        for i, bndry in enumerate(self.exon_bndrys):
            ret_strs.append( "%i\t%s" % ( i, bndry ) )
        
        return "\n".join( ret_strs )
    
    def _find_possible_junctions( self ):
        """
        """
        # I moved this method from an external function, and I don't want to
        # put in all the 'self.'s, so I'm just creating a local variable.
        exon_bndrys = self.exon_bndrys
        
        def junction_is_valid( bndry1, bndry2 ):
            """Stub for future junction filtering.
            
            """
            return True
        
        def iter_junctions( ):
            for i in xrange( len( exon_bndrys ) ):
                # we always include the bin with reads fully within that bin
                for j in xrange( i+1, len( exon_bndrys ) ):
                    # check whether this junction is valid
                    dont_overlap = ( not self.exons_do_overlap( \
                            exon_bndrys[i], exon_bndrys[j] ) )
                    # Stub for a potential junction is valid function
                    valid = junction_is_valid( exon_bndrys[i], exon_bndrys[j]  )
                    if dont_overlap and valid:
                        yield (i,j)
        
        self._possible_junctions = list( iter_junctions() )
        
        return
    
    def iter_possible_junctions( self ):
        if self._possible_junctions == None:
            self._find_possible_junctions()
        
        return iter( self._possible_junctions )
