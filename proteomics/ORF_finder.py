#!/usr/bin/python

# import python mods
import os, sys

import multiprocessing
import Queue
import time
import numpy
import re

from collections import defaultdict, namedtuple
from bx.intervals.intersection import Intersecter, Interval
from itertools import izip
from operator import itemgetter

# import biological mods
from pysam import Fastafile

# declare constants
MIN_AAS_PER_ORF = 120

GENCODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

COMP_BASES = { 'A':'T', 'T':'A', 'C':'G', 'G':'C' }

""" ORF Classes:
## ORF Coding Classes ##
KI - Known Isoform - Same start and stop codons and all splices 
AS = Alternative Splice - Same start and stop codons, but some different splice sites
AA = Alternative stArt - Same stop codon, but different start codon
AO = Alternative stOp - Same start codon, but different stop codon
NI = Novel Isoform - Novel Isoform which shares splice site, but not start or stop codon
SN = alternative Splice No shared splice - Same start and stop codons, but all different spli ce sites
AN = alternative stArt No shared spice - Same stop codon, but different start codon and no shared splice sites
ON = alternative stOp No shared splice - Same start codon, but different stop codon and no shared splice sites
CG = Coding General overlap - Overlaps known ORF, but does not share start or stop codon or a single splice site
NO = NOvel - Does not overlap any known ORF region on same strand

## Non-coding Classes ##
NC = Non-Coding - Does not overlap any known ORF region
NS = Non-Sense - Overlaps known ORF start and stop codon
NG = Non-coding General overlap - Overlaps ORF, but not a known start and stop codon
SO = Short Orf - Overlaps known ORF that is shorter than MIN_AAS_PER_ORF
"""

# Variables effecting .annotation.gtf output
PRINT_NON_ORF_TRANS = True
PRINT_ORF_EXONS = False
PRINT_CODONS = False
PRINT_UTRS = True

AAS_PER_FASTA_LINE = 80
OUTPUT_PROTEINS = False
VERBOSE = False
MIN_VERBOSE = False

DO_PROFILE = False

# add parent(slide) directory to sys.path and import SLIDE mods
sys.path.append( os.path.join( os.path.dirname(__file__), "..", "sparsify" ) )
import sparsify_transcripts
from sparsify_transcripts import build_gene_transcripts
from gene_models import GeneBoundaries
from transcripts import Transcripts, Transcript, TranscriptsFile

sys.path.append( os.path.join(os.path.dirname(__file__), "..", "file_types") )
from gtf_file import parse_gtf_line

class multi_safe_file( file ):
    def __init__( self, *args ):
        args = list( args )
        args.insert( 0, self )
        file.__init__( *args )
        self.lock = multiprocessing.Lock()
        return
    
    def write( self, line ):
        self.lock.acquire()
        file.write( self, line + "\n" )
        self.flush()
        self.lock.release()
        return

class OrfData( dict ):
    def __init__( self ):
        dict.__init__( self )
        
        self.splice_sites = {}
        
        self.ORF = namedtuple( 
            "ORF", [ "int_struct", "is_short", "trans_name" ] )
        return
    
    def add( self, key, orf_struct, short_orf, trans_name ):
        if key in self:
            self[key].append( self.ORF(orf_struct, short_orf, trans_name) )
        else:
            self[key] = [ self.ORF(orf_struct, short_orf, trans_name) ]
        
        return
    
    def are_all_short( self, intervals, chrm, strand ):
        if all( item.is_short for intrvl in intervals for item in self[
                (chrm, strand, intrvl.start, intrvl.end) ] ):
            return True
        return False
    
    def get_trans_name( self, interval, chrm, strand ):
        # return the first transcript name with the specified orf interval
        trans_name = \
            self[(chrm, strand, interval.start, interval.end)][0].trans_name
        return trans_name
    
    def shares_splice( self, interval, trans_orf_struct, chrm, strand ):
        """determine if the trans_orf_struct shares any splice sites in interval
        """
        key = (chrm, strand, interval.start, interval.end)
        if key not in self.splice_sites:
            # get the known orf splice sites for all orfs in interval
            known_splice_sites = set( site for orf in self[key] for 
                                      site in orf.int_struct )
            # cache the set of splice sites for later transcript checks
            self.splice_sites[key] = known_splice_sites
        else:
            known_splice_sites = self.splice_sites[key]
        
        # if there are shared splice site(s) with the known orfs
        if len( known_splice_sites.intersection( trans_orf_struct ) ) > 0:
            return True
        
        return False
    
    def has_exact_match( self, interval, trans_orf_struct, chrm, strand ):
        # loop through the known orf internal structures
        for orf in self[ (chrm, strand, interval.start, interval.end) ]:
            # if an exact match is found
            if trans_orf_struct == orf.int_struct:
                return True
        # if no exact match is found return false
        return False
    
    def get_exact_trans_name( self, interval, trans_orf_struct, chrm, strand ):
        for orf in self[ (chrm, strand, interval.start, interval.end) ]:
            if orf.int_struct == trans_orf_struct:
                return orf.trans_name
        
        return

def write_class_specific_orfs( out_prefix ):
    ann_fname = out_prefix + '.annotated.gtf'
    ki_fname = out_prefix + '.KI_class.gtf'
    as_fname = out_prefix + '.AS_class.gtf'
    
    # get only exact protein matches and same start stop class
    ki_cmd = "awk '$16 ~ /^\"KI/' {0} > {1}".format( ann_fname, ki_fname )
    os.system( ki_cmd )
    as_cmd = "awk '$16 ~ /^\"AS/' {0} > {1}".format( ann_fname, as_fname )
    os.system( as_cmd )
    
    return

def write_orfs_gtf( all_orfs, non_coding_genes, out_prefix ):
    if MIN_VERBOSE:
        print >> sys.stderr, 'Writing annotated CDS files.'
    with TranscriptsFile( out_prefix + '.annotated.gtf', 'w' ) as trans_fp:
        trans_fp.add_transcripts_with_orfs(
            all_orfs, print_non_orfs=PRINT_NON_ORF_TRANS,
            print_orf_exons=PRINT_ORF_EXONS, print_codons=PRINT_CODONS,
            print_UTRs=PRINT_UTRS )
    
    with TranscriptsFile( out_prefix + '.non_coding_genes.gtf', 'w' ) as non_coding_fp:
        non_coding_fp.add_transcripts_with_orfs( 
            non_coding_genes, print_non_orfs=True )
    
    return

def get_trans_class( genomic_region, known_orfs_data, gene, trans ):
    """Classify trans that do not contain an orf as non-coding or non-sense
    """
    # if known orfs are not provided return None
    if known_orfs_data == None: return None
    
    chrm = gene.chromosome
    strand = gene.strand
    
    # unpack known orf structures and use them to find orf match classes
    known_orfs_tree, known_orfs_details = known_orfs_data
    known_orfs = known_orfs_tree[(chrm, strand)].find(
        genomic_region[0], genomic_region[1] )
    
    # if this transcript does not encode an orf and there are no known orfs 
    # overlapping this region then it is a non-coding region
    if len( known_orfs ) == 0: return ('NC', None)
    
    def contains_coord( exon, coord ):
        if exon[0] <= coord and \
                exon[1] >= coord:
            return True
        return False
    
    def contains_start_and_stop_codon( orf_interval ):
        """Determine if a non-coding trans contains both a start and stop codon
        """
        # set initial values to false
        contains_start_codon = False
        contains_stop_codon = False
        
        # loop through each exon to see if it contains a start or stop codon
        for exon in [gene.bp_bndry(i) for i in trans]:
            if contains_coord( exon, orf_interval.start ):
                contains_start_codon = True
            if contains_coord( exon, orf_interval.end ):
                contains_stop_codon = True
            # short circuit the loop once both start and stop overlaps are found
            if contains_start_codon and contains_stop_codon:
                return True
        
        return False
    
    
    # initialize found short orf overlap trans name var
    short_orf_trans_name = None
    # if this is transcript does not have an orf, but overlaps a known orf
    # check whether the transcript overlaps the actual start and stop codons
    for orf_interval in known_orfs:
        if contains_start_and_stop_codon( orf_interval ):
            # if the overlapping orf is short then save the short orf trans name
            if known_orfs_details.are_all_short( 
                    [ orf_interval, ], chrm, strand ):
                short_orf_trans_name = known_orfs_details.get_trans_name( 
                    orf_interval, chrm, strand )
            # if there is a known orf that is not short then this class takes 
            # precident so 'NS' is returned immediately
            else:
                orf_trans_name = known_orfs_details.get_trans_name( 
                    orf_interval, chrm, strand )
                return ( 'NS', orf_trans_name )
    
    # if a short overlap is found, but a 'NS' class was not return it
    if short_orf_trans_name != None:
        return ( 'SO', short_orf_trans_name  )
    
    # get a trans_name for general overlap type
    orf_trans_name = known_orfs_details.get_trans_name( 
        known_orfs[0], chrm, strand )
    # if the transcript does not overlap the start and stop codon then it is 
    # a general overlap of an ORF, but not a non-sense
    return ( 'NG', orf_trans_name )

def get_orf_class( genomic_region, known_orfs_data, gene, trans ):
    """Get orf classes described above compared to known orfs
    """
    # if known orfs are not provided return None
    if known_orfs_data == None: return None
    
    chrm = gene.chromosome
    strand = gene.strand
    
    # unpack known orf structures and use them to find orf match classes
    known_orfs_tree, known_orfs_details = known_orfs_data
    known_orfs = known_orfs_tree[(chrm, strand)].find(
        genomic_region[0], genomic_region[1] )
    
    # if there are no overlapping known orfs then this is a novel orf
    if len( known_orfs ) == 0: return ( 'NO', None )
    
    def get_trans_orf_struct():
        # find the set of boundaries internal to the orf in this trans
        trans_orf_struct = []
        for (start, stop) in [ gene.bp_bndry(i) for i in trans ]:
            # add start and/or stop only if they fall within the trans_orf
            if start > genomic_region[0] and start < genomic_region[1]:
                trans_orf_struct.append( start )
            if stop > genomic_region[0] and stop < genomic_region[1]:
                trans_orf_struct.append( stop )
        # ensure sorted for comparison to known_orf_struct
        trans_orf_struct = tuple(sorted( trans_orf_struct ))
        
        return trans_orf_struct
    
    def shares_start_and_stop( interval ):
        if interval.start == genomic_region[0] and \
                interval.end == genomic_region[1]:
            return True
        return False
    
    def shares_start( interval ):
        if (strand == '+' and interval.start == genomic_region[0]) or \
                (strand == '-' and interval.end == genomic_region[1]):
            return True
        return False
    
    def shares_stop( interval ):
        if (strand == '+' and interval.end == genomic_region[1]) or \
                (strand == '-'and interval.start == genomic_region[0]):
            return True
        return False
    
    # orf classes enumeration
    ( exact_match, same_start_stop, same_stop, same_start, shared_splice,
      start_stop_no_splice, stop_no_splice, start_no_splice,
      general_overlap ) = range(9)
    orf_classes = { exact_match: 'KI', same_start_stop:'AS', same_stop:'AO',
                    same_start:'AA', shared_splice:'NI',
                    start_stop_no_splice:'SN', stop_no_splice:'ON',
                    start_no_splice:'AN', general_overlap:'CG' }
    
    def update_orf_class( interval, orf_class ):
        # if the transcript shares a splice splice with the known orf
        if known_orfs_details.shares_splice(
                interval, trans_orf_struct, chrm, strand ):
            if shares_start_and_stop( interval ):
                if known_orfs_details.has_exact_match(
                        interval, trans_orf_struct, chrm, strand ):
                    return exact_match
                return same_start_stop
            elif orf_class > same_start and shares_start( interval ):
                return same_start
            elif orf_class > same_stop and shares_stop( interval ):
                return same_stop
            elif orf_class > shared_splice:
                return shared_splice
            return orf_class
        
        # this start stop does not share a splice site with this trans
        else:
            if orf_class <= shared_splice:
                return orf_class
            
            if shares_start_and_stop( interval ):
                return start_stop_no_splice
            elif orf_class > start_no_splice and shares_start( interval ):
                return start_no_splice
            elif orf_class > stop_no_splice and shares_stop( interval ):
                return stop_no_splice
        
        return orf_class
    
    
    # initialize orf class to no_match (possible novel isoform)
    orf_class = general_overlap
    # get the trans orf structure from gene and trans object
    trans_orf_struct = get_trans_orf_struct()
    
    # loop through intervals and analyze each one
    for interval in known_orfs:
        orf_class = update_orf_class( interval, orf_class )
        if orf_class == exact_match:
            trans_name = known_orfs_details.get_exact_trans_name( 
                interval, trans_orf_struct, chrm, strand )
            return ( orf_classes[orf_class], trans_name )
    
    trans_name = known_orfs_details.get_trans_name(known_orfs[0], chrm, strand)
    return (orf_classes[orf_class], trans_name)

def translate( seq ):
    """Emulate Biopythons translate method, but faster
    """
    prot_seq = ''
    for codon in [ seq[i*3:(i+1)*3] for i in xrange(int(len(seq)/3)) ]:
        if codon not in GENCODE:
            print >> sys.stderr, 'Invalid codon found: ' + codon
            raise ValueError
        prot_seq += GENCODE[ codon ]
    
    return prot_seq

def format_fasta( gene, genomic_orf, orf_seq, trans_names, orf_class, orf_id ):
    """ output fasta format protein records for input orfs
    """
    try:
        # get protein sequence using translate
        prot = translate( orf_seq )
    except ValueError:
        prot = 'X'
        print >> sys.stderr, 'WARNING: ' + gene.name + ' contains a valid open reading ' + \
            'frame, but contains an invalid base in the provided fasta file.'
    
    trans_names = ';'.join( trans_names )
    # format the orf class with class type and matching transcript
    if orf_class == None:
        class_string = ''
    else:
        class_string = 'class:' + orf_class[0]
        if orf_class[1] != None:
            class_string += ':' + orf_class[1]
    
    # produce fasta line header
    header = '|'.join( (
            '>' + gene.name, 'ORF_id:' + str(orf_id), 'length:'+str(len(prot)),
            'chr' + gene.chromosome + ':' + str(genomic_orf[0]) + '-' +
            str(genomic_orf[1]) + ':' + gene.strand, 
            class_string, trans_names ) )
    fasta_output = header + '\n'
    
    # format protein sequence lines
    pos = 0
    while len( prot ) - pos > AAS_PER_FASTA_LINE:
        fasta_output += str(prot[ pos : pos + AAS_PER_FASTA_LINE ]) + '\n'
        pos += AAS_PER_FASTA_LINE
    fasta_output += prot[pos:]
    
    return fasta_output

def convert_to_genomic( orfs, gene, trans ):
    """ Convert open reading frame coordinates from transcript to genomic space
    """
    # store tuples of genomic exon start positions
    exon_starts = [ gene.bp_bndry(index)[0] for index in trans ]
    # store the length in the transcript up to the exon start
    cumsum_exon_lengths = numpy.cumsum( [
            (gene.bp_bndry(index)[1] - gene.bp_bndry(index)[0] + 1)
            for index in trans] )
    # get trans_len for readability
    trans_len = cumsum_exon_lengths[-1]
    
    def calc_genomic_offset( bps_into_trans ):
        # get index in exon_lengths of exon containing pos
        exon_index = numpy.searchsorted( cumsum_exon_lengths, bps_into_trans )
        # get the genomic start position of that exon
        exon_genomic_start = exon_starts[ exon_index ]
        # get the exon start position along the transcript
        trans_bps_before_exon_start = 0 if exon_index == 0 else \
            cumsum_exon_lengths[ exon_index - 1 ]
        
        return exon_genomic_start + \
            ( bps_into_trans - trans_bps_before_exon_start )
    
    genomic_orfs = []
    for orf_start, orf_stop in orfs:
        # if neg stranded gene convert start and stop to pos strand equvalents
        # when printing this will be resolved using gene.strand
        if gene.strand == '-':
            tmp_start = orf_start
            orf_start = trans_len - orf_stop - 1
            orf_stop = trans_len - tmp_start - 1
        
        genomic_orf_start = calc_genomic_offset( orf_start )
        genomic_orf_stop = calc_genomic_offset( orf_stop )
        
        genomic_orfs.append( (genomic_orf_start, genomic_orf_stop) )
    
    return genomic_orfs

def find_all( sequence, codon ):
    """ Returns a list of positions within sequence that are the start of codon
    """
    return [ x.start() for x in re.finditer( codon, sequence ) ]

def find_orfs( sequence ):
    """ Finds all valid open reading frames in the string 'sequence', and
    returns them as tuple of start and stop coordinates
    """    
    def grp_by_frame( locs ):
        locs_by_frame = { 0: [], 1:[], 2:[] }
        for loc in locs:
            locs_by_frame[ loc%3 ].append( loc )
        for frame in locs_by_frame.keys():
            locs_by_frame[ frame ].reverse()
        return locs_by_frame

    def find_orfs_in_frame( starts, stops ):
        prev_stop = -1
        while len( starts ) > 0 and len( stops ) > 0:
            start = starts.pop()
            if start < prev_stop: continue
            stop = stops.pop()
            while stop < start and len( stops ) > 0:
                stop = stops.pop()
            if start > stop:
                assert len( stops ) == 0
                break
            if (stop - start) >= (MIN_AAS_PER_ORF * 3):
                yield ( start, stop-1 )
            prev_stop = stop

        return

    # find all start and stop codon positions along sequence
    starts = find_all( sequence, 'ATG' )
    stop_amber = find_all( sequence, 'TAG' )
    stop_ochre = find_all( sequence, 'TAA' )
    stop_umber = find_all( sequence, 'TGA' )
    stops = stop_amber + stop_ochre + stop_umber
    stops.sort()
    
    orfs = []    
    # group the starts and stops by their frame
    
    starts_by_frame = grp_by_frame( starts )
    stops_by_frame = grp_by_frame( stops )
    
    
    for frame in (0, 1, 2):
        starts = starts_by_frame[ frame ]
        stops = stops_by_frame[ frame ]
        orfs.extend( find_orfs_in_frame(starts, stops) )
    
    return orfs

def get_trans_seq( trans, gene, gene_seq ):
    """ get the mRNA sequence of the transcript from the gene seq
    """
    trans_seq = ''
    for exon_index in trans:
        # get the exon coords from the gene object
        start, stop = gene.bp_bndry( exon_index )
        
        # convert the coords from genomic to gene-relative
        gene_start = start - gene.boundaries.start
        gene_stop = stop - gene.boundaries.start
        
        if gene.strand == '+':
            # get the portion of the gene sequence for the current exon
            # add 1 to stop since string slice is closed-open
            trans_seq += gene_seq[ gene_start : gene_stop + 1 ]
        else:
            # if the gene is neg strand reverse coords as gene_seq is rev_comp
            tmp_start = gene_start
            gene_start = len(gene_seq) - gene_stop - 1
            gene_stop = len(gene_seq) - tmp_start - 1
            # add the new sequence at the beginning since seq is rev_comp
            trans_seq = gene_seq[ gene_start : gene_stop + 1 ] + trans_seq
    
    return trans_seq

def reverse_complement( seq ):
    """Emulate Biopython reverse_complement method, but faster
    """
    rev_comp_seq = ''
    # loop through sequence in reverse sequence
    for base in seq[::-1]:
        if base in COMP_BASES:
            # else add the compelemntary base
            rev_comp_seq += COMP_BASES[ base ]
        else:
            # if the base is invalid just add it to the rev_comp
            rev_comp_seq += base
    
    return rev_comp_seq

def get_gene_seq( gene, fasta ):
    chrm = gene.chromosome
    if not chrm.startswith( 'chr' ):
        chrm = 'chr' + chrm
    
    # get the raw sequence from the gene object and the fasta file
    # subtract one from start since fasta is 0-based closed-open
    gene_seq = fasta.fetch(
        chrm, gene.boundaries.start - 1, gene.boundaries.stop )
    
    # convert the sequence to upper case
    gene_seq = gene_seq.upper()
    
    if gene.strand == '-':
        gene_seq = reverse_complement( gene_seq )
    
    return gene_seq

def add_trans_orfs( trans, gene, trans_seq, orfs, genomic_orfs,
                    gene_orfs, gene_orf_id, num_trans_w_orf, known_orfs_data ):
    """Add the orfs in this trans to the gene_orfs
    """
    # initiailize the classes and orf_ids lists for trans orfs
    orf_classes = []
    orf_ids = []
    # if there are orfs to be recorded from this trans
    if len( orfs ) > 0:
        num_trans_w_orf += 1
        # add unique orfs to gene_orfs dict
        for i, genomic_region in enumerate(genomic_orfs):
            # create minimal unique orf idenifier consisting 
            # of genomic start stop and orf seq
            key = ( genomic_region, trans_seq[orfs[i][0]:orfs[i][1]] )
            # if this is a trans that has an orf that has already been
            # identified in another transcript just add the tran_name
            if key in gene_orfs:
                gene_orfs[ key ][2].append( trans.name )
            # if this orf has not been identified yet then find its orf
            # class and give it a unique orf_id
            else:
                # determine orf class for this orf
                orf_class = get_orf_class(
                    genomic_region, known_orfs_data, gene, trans )
                # initialize the unique orf entry
                gene_orfs[ key ] = \
                    [ orf_class, gene_orf_id, [trans.name,] ]
                gene_orf_id += 1
            
            # get orf_class and orf_id which were either just calculated
            # or previously stored for this orf
            orf_classes.append( gene_orfs[key][0] )
            orf_ids.append( gene_orfs[key][1] )
    else:
        # get the transcript boundaries to compare to known orfs
        genomic_region = ( gene.bp_bndry( trans.up_strm_exons[0] )[0],
                           gene.bp_bndry( trans.dn_strm_exons[-1] )[1] )
        orf_classes.append( get_trans_class(
                genomic_region, known_orfs_data, gene, trans ) )
    
    return gene_orfs, orf_classes, orf_ids, gene_orf_id, num_trans_w_orf

def apply_orf_rank( gene_orfs, gene_ann_trans ):
    """ convert the orf_ids to a rank of orf lengths
    """
    def apply_orf_rank_to_annotation( gene_ann_trans, orf_ids_map ):
        # loop through transcripts and convert the orf_ids for each
        for i, (trans, gene, genomic_orfs, orf_classes, orf_ids) in \
                enumerate(gene_ann_trans):
            new_orf_ids = []
            for old_orf_id in orf_ids:
                new_orf_ids.append( orf_ids_map[ old_orf_id ] )
            
            # replace old orf ids with size rank orf ids
            gene_ann_trans[i][4] = new_orf_ids
        
        return gene_ann_trans
    
    # get a list of orf lengths for this locus
    orf_lengths = []
    for genomic_region, orf_seq in gene_orfs.iterkeys():
        orf_lengths.append( ( len( orf_seq ),
                              (genomic_region, orf_seq) ) )
    
    # create a map from old_orf_id to orf_length_rank and 
    # convert gene_orf's orf_ids
    orf_ids_map = {}
    for new_orf_id, ( orf_len, key ) in enumerate( 
        sorted( orf_lengths, key=itemgetter(0), reverse=True ) ):
        orf_class, old_orf_id, orf_trans_names = gene_orfs[key]
        gene_orfs[key] = [ orf_class, new_orf_id, orf_trans_names ]
        orf_ids_map[old_orf_id] = new_orf_id
    
    # change the longest orf_ids in gene_ann_trans to 0 as well
    gene_ann_trans = apply_orf_rank_to_annotation( gene_ann_trans, orf_ids_map )
    
    return gene_orfs, gene_ann_trans

def find_gene_orfs( gene, transcripts, fasta, output_queue, known_orfs_data, 
                    protein_fp ):
    """Find all of the unique open reading fromes in a gene
    """
    # get the gene sequence from the fasta file
    # note that the sequence is reverse_complement if neg strand gene
    gene_seq = get_gene_seq( gene, fasta )
    
    # a unique set of ORFs accessed by genomic ORF coords and orf seq
    # this is used to create the .proteins.fa records for this gene
    gene_orfs = {}
    gene_ann_trans = []
    gene_orf_id = 1
    num_trans_w_orf = 0
    for trans in transcripts:
        # get sequence of transcript
        trans_seq = get_trans_seq( trans, gene, gene_seq )
        # find open reading frames
        orfs = find_orfs( trans_seq )
        # convert orf coord from trans-relative to genomic-relative
        genomic_orfs = convert_to_genomic( orfs, gene, trans )
        
        gene_orfs, orf_classes, orf_ids, gene_orf_id, num_trans_w_orf = \
            add_trans_orfs( trans, gene, trans_seq, orfs, genomic_orfs, 
                            gene_orfs, gene_orf_id, num_trans_w_orf, 
                            known_orfs_data )
        
        # add orf info for each transcript even if genomic_orfs is empty
        gene_ann_trans.append(
            [trans, gene, genomic_orfs, orf_classes, orf_ids] )
    
    has_orf = gene_orf_id > 1
    # if any of the transcripts in this gene contain an orf
    if gene_orf_id > 1:
        # change the orf_id to be ordered by length in both data structs
        gene_orfs, gene_ann_trans = apply_orf_rank( \
            gene_orfs, gene_ann_trans )
    
    output_queue.put( ( has_orf, gene_ann_trans, num_trans_w_orf, len(gene_orfs) ) )
    
    if OUTPUT_PROTEINS:
        fasta_lines = []
        # output unique proteins for current gene
        for (genomic_orf, orf_seq), (orf_class, orf_id, trans_names) in \
                gene_orfs.iteritems():
            # store protein sequences
            fasta_output = format_fasta(
                gene, genomic_orf, orf_seq, trans_names, orf_class, orf_id )
            fasta_lines.append( fasta_output )
        
        if len(fasta_lines) > 0:
            # write protein sequences out to multi_safe_file
            protein_fp.write( '\n'.join( fasta_lines ) )
    
    return

def find_gene_orfs_worker( input_queue, output_queue, fasta_fn, 
                           protein_fp, known_orfs_data ):
    # open fasta file in each thread separately
    fasta = Fastafile( fasta_fn )
    
    # process genes for orfs until input queue is empty
    while not input_queue.empty():
        try:
            ( gene, transcripts ) = input_queue.get(block=False)
        except Queue.Empty:
            break
        
        if VERBOSE: print >> sys.stderr, '\tProcessing ' + gene.name
        find_gene_orfs( gene, transcripts, fasta, output_queue, 
                        known_orfs_data, protein_fp )
        if VERBOSE: print >> sys.stderr, '\tFinished ' + gene.name
    
    return

def find_all_orfs( genes, all_trans, fasta_fn, known_orfs_data, 
                   out_prefix, threads ):
    # create queues to store input and output data
    manager = multiprocessing.Manager()
    input_queue = manager.Queue()
    output_queue = manager.Queue()
    
    if MIN_VERBOSE: print >> sys.stderr, 'Processing all transcripts for ORFs.'
    
    # open protein output file
    protein_fp = multi_safe_file( out_prefix + '.proteins.fa', 'w' ) \
        if OUTPUT_PROTEINS else None
    
    total_genes = 0
    total_trans = 0
    # populate input_queue and calculate summary stats
    for gene_name, transcripts in all_trans.iteritems():
        total_genes += 1
        total_trans += len( transcripts )
        
        gene = genes[gene_name]
        
        input_queue.put( ( gene, transcripts ) )
    
    #find_gene_orfs_worker( input_queue, output_queue, fasta_fn, protein_fp, known_orfs_data )
    #sys.exit( 0 )
    #return
    
    args = ( input_queue, output_queue, fasta_fn, protein_fp, known_orfs_data )

    # spawn threads to estimate genes expression
    processes = []
    for thread_id in xrange( threads ):
        p = multiprocessing.Process( target=find_gene_orfs_worker, args=args )
        p.start()
        processes.append( p )
    
    # initialize all_orfs return list 
    # this is used for writing the .annotation.gtf file
    all_orfs = []
    non_coding_genes = []
    
    # initialize summary variables
    genes_w_orf = 0
    trans_w_orf = 0
    num_orfs = 0
    
    def process_queue():
        trans_w_orf, num_orfs, num_genes_w_orf = (0,0,0)
        while not output_queue.empty():
            try:
                ( has_orf, gene_ann_trans, gene_trans_w_orf, num_gene_orfs ) = \
                    output_queue.get(block=False)
            except Queue.Empty:
                break
            
            trans_w_orf += gene_trans_w_orf
            num_orfs += num_gene_orfs
            
            # store gene_ann_trans in appropriate list
            if has_orf:
                all_orfs.extend( gene_ann_trans )
                num_genes_w_orf += 1
            else:
                non_coding_genes.extend( gene_ann_trans )
        
        return trans_w_orf, num_orfs, num_genes_w_orf
    
    # process output queue
    while any( p.is_alive() for p in processes ):
        gene_trans_w_orf, num_gene_orfs, num_genes_w_orf = process_queue()
        genes_w_orf += num_genes_w_orf
        trans_w_orf += gene_trans_w_orf
        num_orfs += num_gene_orfs
        time.sleep(1)
    
    # process any remaining lines
    gene_trans_w_orf, num_gene_orfs, num_genes_w_orf = process_queue()
    genes_w_orf += num_genes_w_orf
    trans_w_orf += gene_trans_w_orf
    num_orfs += num_gene_orfs
    
    if OUTPUT_PROTEINS:
        protein_fp.close()
    
    if MIN_VERBOSE:
        print >> sys.stderr, '\t' + str(genes_w_orf) + ' genes out of ' + \
            str(total_genes) + ' contained an open reading frame.'
        print >> sys.stderr, '\t' + str(trans_w_orf) + ' transcripts out of ' + \
            str(total_trans) + ' contained an open reading frame.'
        print >> sys.stderr, '\t' + 'A total of ' + str(num_orfs) + \
            ' valid unique open reading frames were identified.'
    
    return all_orfs, non_coding_genes

def parse_known_orfs( gtf_fp ):
    """Parse gtf file with known orfs
    orfs are stored as an interval searchable tree and a dict indexed by the
    tree search results
    """
    def get_name_from_field( name ):
        if name.startswith('"'):
            name = name[1:]
        if name.endswith(';'):
            name = name[:-1]
        if name.endswith('"'):
            name = name[:-1]
        return name
    
    
    if gtf_fp == None:
        return None
    
    # store all orf cds regions grouped by their chrm, strand and trans_id
    all_known_orfs = defaultdict( list )
    for line in gtf_fp:
        gtfl = parse_gtf_line( line )
        if gtfl == None or gtfl.feature != 'CDS': continue
        # get the transcript name if it is provided as this is likely the common
        # name. The common name may not be unique though so also store trans_id
        trans_name = gtfl.trans_id
        #trans_name = gtfl.meta_data[ 'transcript_name' ] if \
        #    'transcript_name' in gtfl.meta_data else gtfl.trans_id
        trans_name = get_name_from_field( trans_name )
        # store just the boundary positions here
        # assume that all input orf CDS regions are valid
        all_known_orfs[
            ( gtfl.region.chr, gtfl.region.strand, 
              gtfl.trans_id, trans_name ) ].extend(
            ( gtfl.region.start, gtfl.region.stop ) )
    
    # convert structure to a searchable bx-python tree and 
    # matching internal structure dictionary
    known_orfs_tree = defaultdict( Intersecter )
    known_orfs_details = OrfData()
    for (chrm, strand, trans_id, trans_name), cdss in \
            all_known_orfs.iteritems():
        # determine if this known orf is shorter than MIN_AAS_PER_ORF
        orf_len = 0
        for start, stop in izip(cdss[::2], cdss[1::2]):
            orf_len += stop - start + 1
        short_orf = (orf_len / 3 < MIN_AAS_PER_ORF)
        
        # sort the bndry positions so that first and last are the orf boundaries
        cdss = sorted( cdss )
        # store the orf boundaries in an interval searchable tree
        known_orfs_tree[ (chrm, strand) ].add_interval(
            Interval( cdss[0], cdss[-1] ) )
        # store the rest of the internal orf stucture in a dict indexed with the
        # results from the results from a tree search from known_orfs_tree
        known_orfs_details.add( (chrm, strand, cdss[0], cdss[-1]),
                                tuple( cdss[1:-1] ), short_orf, trans_name )
    
    return known_orfs_tree, known_orfs_details

def build_objects( gtf_fp, ann_fp ):
    if MIN_VERBOSE:
        print >> sys.stderr, 'Parsing input files.'
    
    # load transcripts and gene objects
    genes = GeneBoundaries( gtf_fp )
    gtf_fp.seek(0)
    all_trans = build_gene_transcripts( gtf_fp, genes )
    gtf_fp.close()
    
    known_orf_data = parse_known_orfs( ann_fp )
    if ann_fp != None: ann_fp.close()
    
    return genes, all_trans, known_orf_data

def parse_arguments():
    global MIN_AAS_PER_ORF
    global OUTPUT_PROTEINS
    global VERBOSE
    import argparse
    
    parser = argparse.ArgumentParser(
        description = 'Find open reading frames(ORF) in the input GTF file '
        'and output gtf file with annotated reading frames.' )
    parser.add_argument(
        'gtf', type=file,
        help='GTF file to search for ORFs.' )
    parser.add_argument(
        'fasta',
        help='Fasta file from which to search for open reading frames.' )
    parser.add_argument(
        '--protein', default=False, action='store_true',
        help='Flag to output protein fasta file to ' +
        '<out_prefix>.proteins.fa (default: False).' )
    parser.add_argument(
        '--annotation', '-a', type=file,
        help='GTF file of annotated CDS regions. '
        'This will be used to catagorize ORFs.' )
    parser.add_argument(
        '--min_aas', '-m', type=int,
        help='Number of amino acids to require for an open read frame. ' +
        '(default: {0:d})'.format( MIN_AAS_PER_ORF ))
    
    parser.add_argument(
        '--threads', '-t', type=int, default=1,
        help='Number of threads with which to find ORFs. ' +
        '(default: %(default)d)')
    
    parser.add_argument(
        '--out_prefix', '-o',
        help='Prefix of output file. (default: gtf)')
    parser.add_argument( 
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    if args.min_aas != None: MIN_AAS_PER_ORF = args.min_aas
    
    # create default if no prefix provided or if same as gtf filename
    out_prefix = args.out_prefix if args.out_prefix != None else \
        os.path.basename( args.gtf.name )
    
    # set flag args
    OUTPUT_PROTEINS = args.protein
    VERBOSE = args.verbose
    sparsify_transcripts.VERBOSE = VERBOSE
    
    return args.gtf, args.fasta, args.annotation, args.threads, out_prefix

def main():
    gtf_fp, fasta_fn, ann_fp, threads, out_prefix = parse_arguments()
    genes, all_trans, known_orfs_data = build_objects( gtf_fp, ann_fp )
    
    all_orfs, non_coding_genes = find_all_orfs(
        genes, all_trans, fasta_fn, known_orfs_data, out_prefix, threads )
    
    
    write_orfs_gtf( all_orfs, non_coding_genes, out_prefix )
    
    #write_class_specific_orfs( out_prefix )
    
    return

if __name__ == '__main__':
    if DO_PROFILE:
        import cProfile
        cProfile.run( 'main()' )
    else:
        main()
