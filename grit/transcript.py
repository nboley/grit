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

import os, sys

import tempfile
import cPickle as pickle

from itertools import chain, izip
from collections import namedtuple, defaultdict

import files.gtf
import config

GenomicInterval = namedtuple('GenomicInterval', 
                             ['chr', 'strand', 'start', 'stop'])

CategorizedExons = namedtuple( "CategorizedExons", ["TSS", "TES", "internal", 
                                                    "full_transcript"] )

VERBOSE = False

def partition_coding_and_utr_segments( exons, cds_start, cds_stop ):
    """Split the exons into UTR and CDS exons.

    """
    # find the exon index the the start codon intersects
    cds_start_i = [ i for i, (start, stop) in enumerate(exons) 
                    if cds_start >= start and cds_start <= stop ]
    assert len( cds_start_i ) == 1
    cds_start_i = cds_start_i[0]
    
    # we start at the cds_start exon because we know this index must be >= 
    assert cds_stop >= cds_start
    cds_stop_i = [ i for i, (start, stop) in enumerate(exons[cds_start_i:]) 
                   if cds_stop >= start and cds_stop <= stop ]
    assert len( cds_stop_i ) == 1
    cds_stop_i = cds_stop_i[0] + cds_start_i
    
    def mod_external_bndrys( exons, lower_bnd, upper_bnd ):
        """If necessary, shrink the external boundaries"""
        # if it's an empty set, there is nothing to be done
        if len( exons ) == 0: return exons
        exons[0] = ( max(exons[0][0], lower_bnd), exons[0][1] )
        exons[-1] = ( exons[-1][0], min(upper_bnd, exons[-1][1] ) )
        return exons
    
    cds_exons = mod_external_bndrys( 
        list(exons[cds_start_i:cds_stop_i+1]), cds_start, cds_stop )
    
    us_utr_stop_i = cds_start_i if exons[cds_start_i][0] < cds_start \
        else cds_start_i - 1
    us_utr_exons = mod_external_bndrys(
        list(exons[:us_utr_stop_i+1]), 1, cds_start-1)
    
    ds_utr_stop_i = cds_stop_i if exons[cds_stop_i][1] > cds_stop \
        else cds_stop_i + 1
    ds_utr_exons = mod_external_bndrys(
        list(exons[ds_utr_stop_i:]), cds_stop+1, 1e100)
    
    return us_utr_exons, cds_exons, ds_utr_exons

class Gene( object ):
    def __init__(self, id, name, 
                 chrm, strand, start, stop, 
                 transcripts, meta_data={}):
        self.id = id
        self.name = name
        
        self.chrm = chrm
        self.strand = strand
        self.start = start
        self.stop = stop
        self.transcripts = transcripts
        self.meta_data = meta_data
        return    

    def find_nonoverlapping_boundaries( self ):
        boundaries = set()
        for transcript in self.transcripts:
            for exon in transcript.exons:
                boundaries.add( exon[0] )
                boundaries.add( exon[1]+1 )
        
        return sorted( boundaries )
    
    def write_to_file(self, ofname=None):
        if ofname == None:
            opdir = tempfile.mkdtemp()
            ofname = os.path.join(opdir, self.id + ".gene")
        with open(ofname, "w") as ofp:
            pickle.dump(self, ofp)
        return ofname
    
    def find_transcribed_regions( self ):
        exons = set()
        for transcript in self.transcripts:
            exons.update( transcript.exons )
        
        return files.gtf.flatten( sorted( exons ) )
    
    def calc_bpkm(self, read_cov):
        base_cov, length = 0.0, 0
        for start, stop in self.find_transcribed_regions():
            length += stop - start + 1
            base_cov += read_cov[(self.chrm, self.strand)][start:stop+1].sum()
        
        return base_cov/length
        
    def extract_elements(self):
        """Extract the different element types.
        
        """
        elements = {'gene': set(),
                    'internal_exon': set(),
                    'intron': set(),
                    'polya': set(),
                    'promoter': set(),
                    'single_exon_gene': set(),
                    'tes_exon': set(),
                    'tss_exon': set(),
                    'exon': set() }
        elements['gene'].add( (self.start, self.stop) )
        for t in self.transcripts:
            if len( t.exons ) == 1:
                elements['single_exon_gene'].add(
                    (t.exons[0][0], t.exons[0][1]))
            else:
                elements['exon'].update(t.exons)
                tss_exon = t.exons[0] if t.strand == '+' else t.exons[-1]
                elements['tss_exon'].add(tss_exon)
                elements['internal_exon'].update( t.exons[1:-1] )
                elements['intron'].update( t.introns )
                tes_exon = t.exons[-1] if t.strand == '+' else t.exons[0] 
                elements['tes_exon'].add(tes_exon)

            promoter = t.find_promoter()
            elements['promoter'].add(promoter)

            polya = t.find_polya_region()
            elements['polya'].add(polya)

        return elements

class Transcript( object ):
    def __init__(self, trans_id, chrm, strand, exons, cds_region,
                 gene_id, score=None, fpkm=None, fpk=None, 
                 promoter=None, polya_region=None, coding_sequence=None,
                 conf_lo=None, conf_hi=None, frac=None,
                 gene_name=None, name=None):
        self.gene_id = gene_id
        self.id = trans_id

        self.gene_name = gene_name
        self.name = name

        self.ref_gene = None
        self.ref_trans = None
        self.ref_match_class_code = None
        
        self.chrm = chrm
        self.strand = strand
        
        self.score = score
        self.fpkm = fpkm
        self.fpk = fpk
        self.conf_lo = conf_lo
        self.conf_hi = conf_hi
        self.frac = frac
        
        exon_bnds = list( chain( *exons ) )
        self.exon_bnds = exon_bnds
        self.start = exon_bnds[0]
        self.stop = exon_bnds[-1]
        assert self.start <= self.stop

        self.exons = tuple(zip(exon_bnds[:-1:2], exon_bnds[1::2]))
        self.introns = tuple([ (x+1, y-1) for x, y in 
                               izip(exon_bnds[1:-2:2], 
                                    exon_bnds[2:-1:2]) ])
        
        self.is_protein_coding = ( cds_region != None )
        self.coding_sequence = None
        
        self.cds_region = cds_region
        self.start_codon = None
        self.stop_codon = None
        
        self.cds_exons = None
        self.fp_utr_exons = None
        self.tp_utr_exons = None
        self.us_exons = None
        self.ds_exons = None
        
        self.promoter = promoter
        self.polya_region = polya_region
        
        self._seq = None
        
        if cds_region != None:
            self.add_cds_region( cds_region, coding_sequence )
    
    def add_cds_region( self, cds_region, coding_sequence=None ):
        self.is_protein_coding = True
        self.coding_sequence = coding_sequence
        self.cds_region = cds_region
        self.us_exons, self.cds_exons, self.ds_exons = \
            partition_coding_and_utr_segments( 
            self.exons, self.cds_region[0], self.cds_region[1] )
    
        us_codon, ds_codon = self.cds_region[0], self.cds_region[1]
        # if this is a reverse strand transcript, rev 5' and 3' ends
        if self.strand == '+':
            self.fp_utr_exons, self.tp_utr_exons \
                = self.us_exons, self.ds_exons
            self.start_codon, self.stop_codon = us_codon, ds_codon
        else:
            self.fp_utr_exons, self.tp_utr_exons \
                = self.ds_exons, self.us_exons
            self.stop_codon, self.start_codon = us_codon, ds_codon
        
    def __hash__( self ):
        if self.cds_region != None:
            return hash(( self.chrm, self.strand, 
                          self.exons, tuple(self.cds_region) ))
        else:
            return hash( (self.chrm, self.strand, self.exons, None) )
    
    def IB_key( self ):
        """Return a key for matching transcripts on their external bounds.
        
        """
        if len( self.exon_bnds ) == 2:
            return ( self.chrm, self.strand, "SE_GENE", self.cds_region )
        else:
            return (self.chrm, self.strand, 
                    tuple(self.exon_bnds[1:-1]), self.cds_region)
    
    def relative_pos( self, genome_coord ):
        """Convert a genme coordinate into a transcript coordinate.
        
        """
        return genome_coord - self.start - sum( 
            i_stop - i_start + 1 
            for i_start, i_stop in self.introns 
            if i_stop < genome_coord )

    def genome_pos( self, trans_coord ):
        """Convert a transcript coordinate into a genome coordinate.
        
        To do this, we simply need to count the total intron length
        before this point in the transcript.
        """
        # store the total length of exons before the current point
        # in the below loop
        tot_exons_len = 0
        # store the total intron length, before the points in the 
        # below loop
        insert_len = 0
        for (i_start, i_stop), (e_start, e_stop) \
                in izip( self.introns, self.exons ):
            e_len = e_stop - e_start + 1
            # if adding this exon to the exon lengths would move us 
            # past the point in transcription coordinates, we are done
            if tot_exons_len + e_len > trans_coord:
                break
            
            # otherwise, update the current locations
            tot_exons_len += e_len
            
            i_len = i_stop - i_start + 1
            insert_len += i_len
        
        # the location in genome coordinates is simply the location 
        # of the transcript, plus the number of spliced out bases( introns ) 
        # before this position in the transcripts, plus the number of bases we 
        # are into the transcript
        return trans_coord + self.start + insert_len
    
    def build_gtf_lines( self, meta_data, source='.'):
        assert self.gene_id != None
        assert self.id != None
        if self.gene_name != None and self.gene_name != self.gene_id:
            meta_data['gene_name'] = self.gene_name
        if self.name != None and self.name != self.id:
            meta_data['transcript_name'] = self.name
        ret_lines = []
        def build_lines_for_feature( exons, feature, is_CDS=False ):
            current_frame = 0
            score = str(self.score) if self.score != None else '.'
            for start, stop in exons:
                region = GenomicInterval( self.chrm, self.strand, start, stop )
                frame = current_frame if is_CDS else '.'
                yield files.gtf.create_gtf_line( 
                    region, self.gene_id, self.id, meta_data,
                    score, feature=feature, frame=str(frame),
                    source=source)
                current_frame = ( current_frame + stop - start + 1 )%3
            return
        
        if self.promoter != None:
            ret_lines.extend( build_lines_for_feature( 
                    [self.promoter,], 'promoter', False ) )

        ret_lines.extend( build_lines_for_feature( 
                self.exons, 'exon', False ) )
    
        if self.polya_region != None:
            ret_lines.extend( build_lines_for_feature( 
                    [self.polya_region,], 'polya', False ) )
        
        if self.cds_region != None:
            us_exons, ds_exons = self.fp_utr_exons, self.tp_utr_exons
            us_label, ds_label = 'five_prime_UTR', 'three_prime_UTR'
            if self.strand == '-': 
                us_exons, ds_exons = ds_exons, us_exons
                us_label, ds_label = ds_label, us_label
            
            ret_lines.extend( build_lines_for_feature( 
                    us_exons, us_label, False ) )
                    
            ret_lines.extend( build_lines_for_feature( 
                    self.cds_exons, 'CDS', True ) )
            
            ret_lines.extend( build_lines_for_feature( 
                    ds_exons, ds_label, False ) )

        fp_element_order = {'promoter': 1, 'five_prime_UTR': 2, 
                            'exon': 3, 'CDS': 4, 
                            'three_prime_UTR': 5, 'polya': 6}
        tp_element_order = {'promoter': -1, 'five_prime_UTR': -2, 
                            'exon': -3, 'CDS': -4, 
                            'three_prime_UTR': -5, 'polya': -6}
        
        def order(line):
            data = line.split()
            if self.strand == '-':
                return tp_element_order[data[2]], -int(data[3])
            else:
                return fp_element_order[data[2]], -int(data[3])

        ret_lines.sort(key=order)
        
        return "\n".join( ret_lines )
    
    def calc_length(self):
        return sum( x[1]-x[0]+1 for x in self.exons )

    def build_tracking_line(self):
        contig_name = self.chrm
        if config.FIX_CHRM_NAMES_FOR_UCSC:
            contig_name = fix_chrm_name_for_ucsc(contig_name)
        return [
            t.id, # tracking ID
            '-' if t.class_code == None else t.class_code,  # class code
            '-' if t.ref_trans == None else t.ref_trans, # nearest ref id
            t.gene_id, # gene unique id
            '-' if t.ref_gene == None else t.ref_gene, # nearest reference gene
            '-', # TSS ID
            "%s:%s:%i-%i" % (contig_name, t.strand, t.start, t.stop), 
            str(t.calc_length()) # transcript length
            ]
    
    def find_promoter(self, inferred_promoter_length=50):
        # if there is an annotated promoter, then return it
        if self.promoter != None:
            return self.promoter
        
        # otherwise, return *up to* the first inferred_promoter_length bases
        # at the beggining of the transcript
        if self.strand == '+':
            return ( self.exons[0][0],
                     min(self.exons[0][0]+inferred_promoter_length-1, 
                         self.exons[0][1]) )
        else:
            return ( max( self.exons[-1][1]-inferred_promoter_length+1, 
                          self.exons[-1][0] ),
                     self.exons[-1][1] )
        
        assert False

    def find_polya_region(self, inferred_polya_region_length=20):
        # if there is an annotated polya region, then return it
        if self.polya_region != None:
            return self.polya_region
        
        # otherwise, return *up to* the first inferred_promoter_length bases
        # at the beggining of the transcript
        if self.strand == '+':
            return ( max( self.exons[-1][1]-inferred_polya_region_length+1, 
                          self.exons[-1][0] ),
                     self.exons[-1][1] )
        else:
            return ( self.exons[0][0],
                     min(self.exons[0][0]+inferred_polya_region_length-1, 
                         self.exons[0][1]) )
        
        assert False
