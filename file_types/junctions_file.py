import numpy

from itertools import product, izip

from genomic_intervals import GenomicIntervals, GenomicInterval
from gtf_file import parse_gff_line, create_gff_line
from pysam import Fastafile
import itertools
from collections import defaultdict

CONSENSUS_PLUS = 'GTAG'
CONSENSUS_MINUS = 'CTAC'

def get_jn_type( chrm, upstrm_intron_pos, dnstrm_intron_pos, fasta, jn_strand ):
    # get first 2 bases from 5' and 3' ends of intron to determine 
    # splice junction strand by comparing to consensus seqs
    # subtract one from start since fasta is 0-based closed-open
    intron_seq = fasta.fetch( 'chr'+chrm , upstrm_intron_pos - 1, upstrm_intron_pos + 1) + \
        fasta.fetch( 'chr'+chrm , dnstrm_intron_pos - 2 , dnstrm_intron_pos )
    
    # if junction matchs consensus set strand
    # else return non-canonical
    if intron_seq.upper() == CONSENSUS_PLUS:
        canonical_strand = '+'
    elif intron_seq.upper() == CONSENSUS_MINUS:
        canonical_strand = '-'
    else:
        return 'non-canonical'
    
    if jn_strand == canonical_strand:
        return 'canonical'
    return 'canonical_wrong_strand'

def build_jn_line( region, group_id, count=0, fasta_obj=None ):        
    group_id_str = 'group_id "{0}";'.format( str(group_id) )
    if fasta_obj != None:
        intron_type = get_jn_type( region.chr, region.start, region.stop, \
                                       fasta_obj, region.strand )
        group_id_str += ' type "{0}";'.format( intron_type )
    
    jn_line = create_gff_line( region, group_id_str, score=count )
    return jn_line

def write_junctions( junctions, out_fp, scores=None, groups=None, \
                         fasta_fn=None, track_name="discovered_jns" ):
    """Output junctions with more than MIN_NONOVERLAPPING_READS in gff 
       format if apply_filter
    
    """
    if isinstance( junctions, Junctions ):
        assert scores == None and groups == None
        jns_iter = junctions.iter_jns_and_cnts_and_grps()
    else:
        if scores == None: scores = itertools.repeat( 0 )
        if groups == None: groups = itertools.repeat( '.' )
        jns_iter = izip( junctions, scores, groups )
        
    fasta_obj = None if fasta_fn == None else Fastafile( fasta_fn )
    
    out_fp.write( "track name={0}\n".format(track_name) )
    for region, count, grp in jns_iter:
        jn_line = build_jn_line( region, count, grp, fasta_obj )
        out_fp.write( jn_line + '\n' )
    
    return

def get_intron_bndry_sets( jns_fp ):
    """Get a dict of sets containing jns for each chrm, strand pair.
    
    The first returned value, upstream_jns, refers to the side of the
    intron that is *nearest* the promoter. 
    
    The eecond returned value, downstream_jns, refers to side
    nearest the transcript's 3' end. 
    """
    upstream_jns = defaultdict( set )
    downstream_jns = defaultdict( set )
    jns = defaultdict( set )
    for line in jns_fp:
        data = parse_gff_line( line )
        if data == None: continue
        pos = data[0]
        if pos.strand == '+':
            upstream_jns[ (pos.chr, pos.strand) ].add( pos.start  )
            downstream_jns[ (pos.chr, pos.strand) ].add( pos.stop  )
        else:
            assert pos.strand == '-'
            upstream_jns[ (pos.chr, pos.strand) ].add( pos.stop  )
            downstream_jns[ (pos.chr, pos.strand) ].add( pos.start  )

    return upstream_jns, downstream_jns


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
        self._scores = defaultdict( list )
        self._groups = defaultdict( list )
        
        return
    
    def add( self, chrm, strnd, start, stop, score=0, group='.' ):
        GenomicIntervals.add( self, chrm, strnd, start, stop )
        self._scores[ (chrm, strnd) ].append( score )
        self._groups[ (chrm, strnd) ].append( group )
        
        return
    
    def iter_jns_and_cnts_and_grps( self ):
        for (chrm, strand), chrm_jns in self.iteritems():
            for (start, stop), score, grp in izip( \
                    chrm_jns, self._scores[(chrm, strand)], \
                        self._groups[(chrm, strand)] ):
                yield ( GenomicInterval( chrm, strand, int(start), int(stop) ), \
                            int(score), grp )
        
        return
        
    def freeze( self ):
        """Convert the intervals into a sorted numpy array, for faster access.
        Also order scores to maintain indices between intervals
        """
        new_grps = defaultdict( list )
        for key in self:
            # order start, stops by start coordinate
            self[ key ] = numpy.array( self[ key ] )
            sorted_indices = self[ key ][:,0].argsort()
            self[ key ] = self[ key ][ sorted_indices ]
            
            # use same sorted_indices to sort _scores in the same order
            self._scores[ key ] = numpy.array( self._scores[ key ] )
            self._scores[ key ] = self._scores[ key ] [sorted_indices ]

            for i in sorted_indices:
                new_grps[key].append( self._groups[key][ i ] )
        
        self._groups = new_grps
        self.is_frozen = True
        
        return

def add_jns_from_gff_file( jns, junctions_fp ):
    for line in junctions_fp:
        gff = parse_gff_line( line )
        # skip invlaid lines
        if gff == None: continue

        jns.add( gff.region.chr, gff.region.strand, gff.region.start, \
                     gff.region.stop, gff.score, gff.group )
    
    return

def parse_junctions_file_dont_freeze( junctions_fp ):
    jns = Junctions()
    add_jns_from_gff_file( jns, junctions_fp )
    return jns

def parse_junctions_file( junctions_fp ):
    jns = parse_junctions_file_dont_freeze( junctions_fp )
    jns.freeze()
    return jns

def parse_junctions_files( junctions_fps ):
    jns = Junctions()
    for junctions_fp in junctions_fps:
        add_jns_from_gff_file( jns, junctions_fp )
    
    jns.freeze()
    return jns
