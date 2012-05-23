"""Methods for building transcript models. 

Outputs tracking and combined gtf file information
"""
VERBOSE = True

USE_ALL_POSSIBLE_STARTS = False
USE_ALL_POSSIBLE_STOPS = False
FILTER_SPLICED_UTR_EXONS = True
USE_DISCOVERED_JNS = True
USE_PROVIDED_JNS = True
ALWAYS_INCLUDE_FULL_MODEL = False
FIX_GTF_FOR_UCSC_BROWSER = True

import copy
from operator import attrgetter, itemgetter
from itertools import product
from collections import defaultdict, namedtuple
import sys
import os
from scipy import optimize
import multiprocessing

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), "../" ) )
import sparsify_transcripts
FREQ_FILTER = sparsify_transcripts.FREQ_FILTER

class Transcript( object ):
    """Store a list of exons, and their associated boundaries.

    """
    def __init__( self, up_strm_exons, dn_strm_exons, internal_exons=(), name=None ):
        self._hash = None
        self.internal_exons = tuple( sorted( internal_exons ) )
        self.up_strm_exons = tuple( sorted( up_strm_exons ) )
        self.dn_strm_exons = tuple( sorted( dn_strm_exons ) )
        self.name = name
    
    def iter_possible_fragment_types( self ):
        """Iterate through all possible fragment classes.
        
        A 'fragment type' groups fragments into the exons that they cover. So,
        all fragments that start in exon 1 and end in exon 2 are fragment type
        ( 1, 2 ).
        """
        # we need to materialise the exons because of the UTR mess.
        exons = list( self )
        for index_1 in xrange(len(exons)):
            for index_2 in xrange(index_1, len(exons)):
                yield tuple(exons[index_1:index_2+1])
        
        return
    
    def build_nonoverlapping_transcript( self, gene ):
        """Build the equivalent transcript identified by non-overlapping exons.
        
        """
        nonoverlapping_exons = []
        for exon in self:
            nonoverlapping_exons.extend( gene.find_nonoverlapping_exons( exon ) )
        
        return Transcript( (nonoverlapping_exons[0],),  \
                           (nonoverlapping_exons[-1],), \
                            nonoverlapping_exons[1:-1] )
            
    def __iter__( self ):
        """Iterate through each exon of the longest possible transcript.
        
        The longest trancript is always the first fp end, and the last 3p 
        end due to sorting, unless it is a single exon transcript.
        """
        yield self.up_strm_exons[0]
        for e in self.internal_exons:
            yield e
        
        # check if the 5' and 3' utrs are the same
        if self.up_strm_exons[0] != self.dn_strm_exons[-1]:
            yield self.dn_strm_exons[-1]
        
        return
    
    def __hash__( self ):
        if self._hash == None:
            self._hash = hash( tuple( self ) )
        
        return self._hash
    
    def __cmp__( self, other ):
        return hash(self) - hash(other)
    
    def __eq__( self, other ):
        return True if hash(self) == hash(other) else False
    
    def __len__( self ):
        if self.up_strm_exons[0] == self.dn_strm_exons[-1]:
            return 1
        return 2 + len( self.internal_exons )

    def __str__( self ):
        rv = "[{0}]{1}[{2}]".format( ",".join(map( str, self.up_strm_exons )), \
                                         ",".join(map( str, self.internal_exons )), \
                                         ",".join(map( str, self.dn_strm_exons )) )
        return rv
    
    def get_str_identifier( self ):
        # get all of the exons in this transcript
        exons = set( self )
        assert len( exons ) > 0
        
        # return a list of '0''s and '1's, where a one correpsonds to the exon
        # in that index existing within this trancript.
        identifier = [ str(int( i in exons )) for i in xrange( max( exons ) + 1 ) ]
        assert len( identifier ) > 0
        
        return str( int( '1' + "".join( identifier ), 2 ) )
    
    def __repr__( self ):
        return str( self )
    
class Transcripts( list ):
    """Store a collection of transcripts.

    """
    MetaDataTuple = namedtuple('MetaDataTuple', \
                               ['sourcefile', 'lasso_lambda', 'freq', \
                                'read_coverage', 'improbability'])
    
    def __init__( self, gene, iterator=() ):
        """Initialize a transcripts object.
        
        Iterator should return either: transcript objects
        or tuples of the form ( transcript, (meta_data...) )
        where meta_data is a tuple ( sourcefile, lasso_lambda, freq )
        
        """
        self.gene = gene
        
        self.sourcefiles = []
        self.lasso_lambdas = []
        self.freqs = []
        self.read_coverages = []
        self.improbabilities = []
        
        self._frozen = False
        
        for item in iterator:
            if isinstance( item, Transcript ):
                self.add_transcript( item )
            else:
                assert len( item ) == 2
                self.add_transcript( item[0], item[1][0], \
                                     item[1][1], item[1][2], \
                                     item[1][3], item[1][4] )
        
        return

    def append( self, item ):
        raise NotImplementedError( "Use the add transcript method." )
    
    def add_transcript( self, transcript, sourcefile=None, \
                        lasso_lambda=0.00, freq=0.00, read_coverage=0.0, log_prb=0.0 ):
        assert not self._frozen
        if not isinstance( transcript, Transcript ):
            raise ValueError, "Transcripts can only store items of type Transcript"
        list.append( self, transcript )
        self.sourcefiles.append( sourcefile )
        self.lasso_lambdas.append( lasso_lambda )
        self.freqs.append( freq )
        self.read_coverages.append( read_coverage )
        self.improbabilities.append( log_prb )
        
    def add_frequencies( self, freqs ):
        assert len(self) == len(freqs)
        self.freqs = freqs

    def add_read_coverages( self, read_coverages ):
        assert len(self) == len(read_coverages)
        self.read_coverages = read_coverages
    
    def iter_transcripts_and_freqs( self ):
        for t, f in zip( self, self.freqs ):
            yield t, f
        return

    def iter_transcripts_and_metadata( self ):
        for t, m in zip( self, zip( self.sourcefiles, self.lasso_lambdas, \
                                    self.freqs, self.read_coverages, \
                                    self.improbabilities ) ):
            yield t, self.MetaDataTuple( *m )
        return

    def iter_connected_exons( self ):
        paired_exons = set()
        # iterate through each transcript in self
        for trans in self:
            # add every consecutive exon pair to connected exons set
            prev_exon_i = trans.up_strm_exons[0]
            for i, exon_i in enumerate( trans ):
                # add only the pairs after the first exon
                if i > 0:
                    paired_exons.add( (prev_exon_i, exon_i) )
                prev_exon_i = exon_i
        return sorted( paired_exons )
    
    def clear(self):
        del self[:]
        self.sourcefiles = []
        self.lasso_lambdas = []
        self.freqs = []
        self.read_coverages = []
        self.improbabilities = []
        self._frozen = False
        
    def sort( self, order='lambda' ):
        data = list( self.iter_transcripts_and_metadata() )
        
        if order == 'lambda':
            data.sort( key = lambda x:x[1].freq, reverse=True )
            data.sort( key = lambda x:x[1].lasso_lambda, reverse=True )
        # otherwise, by frequency
        elif order == 'freq':
            data.sort( key = lambda x:x[1].lasso_lambda, reverse=True )
            data.sort( key = lambda x:x[1].freq, reverse=True )
        else:
            raise ValueError, "Unrecognized sort order: {0} " + \
                "( should be lambda or freq )".format( order )
        
        self.clear()

        for item in data:
            self.add_transcript( item[0], item[1][0], item[1][1], \
                                 item[1][2], item[1][3], item[1][4] )
        
    def pretty_str( self, filter_0_freq=True ):
        if len( self ) == 0:
            return "\n"
        
        # find the maximum length of the transcript string
        max_len = max( len( str(tr) ) for tr in self )

        strs = []
        for t, freq in self.iter_transcripts_and_freqs():
            # we check for None, as the frequency may not have
            # been estimated
            if not filter_0_freq or round( freq, 8 ) > 0:
                strs.append( str(t).ljust( max_len + 5 ) \
                             + str(round( freq, 3 )) )
        
        return "\n".join( strs )
    
    def __hash__( self ):
        self._frozen = True
        return hash( tuple( self ) )


class TranscriptsFile( file ):
    """A multiprocess safe version of a file.
    """
    def __init__( self, *args ):
        args = list( args )
        args.insert( 0, self )
        file.__init__( *args )
        self.write_lock = multiprocessing.Lock()
        self.write(\
            "track name=\"SLIDE\" description=\"SLIDE\" visibility=full\n")
        self.flush()
    
    @staticmethod
    def _build_exon_gtf_line( exon_gene_index, exon_trans_index, \
                                  gene, transcript_id, feature_type='exon', \
                                  region=None, score='.', frame='.', \
                                  lasso_lambda=None, freq=None, \
                                  source_file=None, readcnt=None, orf_id=None,\
                                  orf_classes=None):
        """
        seqname - The name of the sequence. Must be a chromosome or scaffold.
        source  - The program that generated this feature.
        feature - The name of this type of feature. Some examples of 
                  standard feature types are "CDS", "start_codon", 
                  "stop_codon", and "exon".
        start   - The starting position of the feature in the sequence. 
                  The first base is numbered 1.
        end     - The ending position of the feature (inclusive).
        score   - A score between 0 and 1000. If the track line useScore 
                  attribute is set to 1 for this annotation data set, the score 
                  value will determine the level of gray in which this feature 
                  is displayed (higher numbers = darker gray). If there is no 
                  score value, enter ".".
        strand  - Valid entries include '+', '-', or 
                  '.' (for don't know/don't care).
        frame   - If the feature is a coding exon, frame should be a number 
                  between 0-2 that represents the reading frame of the first 
                  base. If the feature is not a coding exon, 
                  the value should be '.'.
        group   - All lines with the same group are linked together into a 
                  single item. 
        """
        gtf_line = list()
        if FIX_GTF_FOR_UCSC_BROWSER:
            chr_name = gene.chromosome
            if not chr_name.startswith( "chr" ):
                chr_name = "chr" + chr_name
            gtf_line.append( chr_name )
        else:
            gtf_line.append( gene.chromosome )
        
        gtf_line.append( "slide" )
        gtf_line.append( feature_type )
        if region != None:
            start, stop = region
        else:
            start, stop = gene.bp_bndry( exon_gene_index )
        gtf_line.append( str(start))
        assert start > 0
        gtf_line.append( str(stop))
        gtf_line.append( str(score) )
        gtf_line.append( gene.strand )
        gtf_line.append( str(frame) )
        
        meta_info_pat = \
            'gene_id "{0:s}"; transcript_id "{1:s}"; exon_number "{2:d}"'
        meta_info = meta_info_pat.format( \
            gene.name, transcript_id, exon_trans_index )
        if lasso_lambda != None:
            meta_info += '; lasso_lambda "{0:.7f}"'.format( lasso_lambda )
        if freq != None:
            meta_info += '; freq "{0:.7f}"'.format( freq )
        if source_file != None:
            meta_info += '; source "{0:s}"'.format( source_file )
        if readcnt != None:
            meta_info += '; readcnt "{0:d}"'.format( readcnt )
        if orf_id != None:
            meta_info += '; orf_id "{0:d}"'.format( orf_id )
        if orf_classes != None:
            orf_classes = '|'.join( ( ':'.join( map( str, item) ) for \
                                          item in orf_classes if item != None ))
            meta_info += '; orf_class "{0}"'.format( orf_classes )
        gtf_line.append( meta_info )
        
        return '\t'.join( gtf_line )
    
    def add_transcripts( self, transcripts, gene, include_meta_data=True ):
        """write out each transcript.
        
        """
        ## FIXME BUG XXX
        readcnt = 1
        
        all_transcript_lines = []
        for trans_index, (trans, md) \
            in enumerate(transcripts.iter_transcripts_and_metadata()):

            trans_id = gene.name + "_" + trans.get_str_identifier()
            if len( trans_id ) > 200:
                # make sure the transcript ids aren't too long
                trans_id = str( hash(trans_id) )
            # a hash should never be longer than a word, which ( at least for the 
            # forseeable future ) should be shorter than 254 bytes 
            assert len( trans_id ) < 254
            
            # need to enumerate so that transcript has sequentially numbered exons
            # not exons numbered (arbitrarily) by gene exon number
            for exon_trans_index, exon_gene_index in enumerate( trans ):
                # for the merged form, we want to include meta data
                if include_meta_data:
                    all_transcript_lines.append( \
                        self._build_exon_gtf_line( \
                            exon_gene_index, exon_trans_index+1, transcripts.gene, \
                                trans_id, lasso_lambda=md.lasso_lambda, freq=md.freq, \
                                source_file=md.sourcefile, readcnt=readcnt ) )
                else:
                    all_transcript_lines.append( 
                        self._build_exon_gtf_line( \
                            exon_gene_index, exon_trans_index+1, gene, trans_id ) )

        self.write_lock.acquire()
        for line in all_transcript_lines:
            self.write( line + "\n" )
        self.flush()
        self.write_lock.release()
        
        return
    
    def add_unique_transcripts( self, all_transcripts, \
                                    filter_by_nonzero_freq=False ):
        """write each transcript present in any sample from a single gene

        """
        gene = None
        
        unique_transcripts = set()
        for transcripts in all_transcripts:
            # find the gene name, and make sure that it is the same for all 
            # of the transcripts objects
            if gene == None:
                gene = transcripts.gene
            else:
                assert gene.name == transcripts.gene.name
            
            for transcript, md in transcripts.iter_transcripts_and_metadata():
                if filter_by_nonzero_freq and md.freq < FREQ_FILTER:
                    continue
                
                unique_transcripts.add( transcript )
        
        unique_transcripts = Transcripts( \
            gene, sorted(unique_transcripts, key=lambda x:str(x)) )
        
        self.add_transcripts( unique_transcripts, gene, False )
        
        return
    
    @staticmethod
    def _get_orf_regions( gene, trans, orfs, orf_classes, orf_ids, \
                              print_UTRs, print_codons ):
        """Get all regions corresponding to each open reading frame
        """
        # initialize list of orf regions
        trans_orf_regions = []
        
        reverse_strand = True if gene.strand == '-' else False
        # set upstream and downstream feature type values 
        # for gtf lines based on strandedness
        feature_strings = \
            { 'up':('3UTR','stop_codon'), 'dn':('5UTR','start_codon') } \
            if reverse_strand else \
            { 'up':('5UTR','start_codon'), 'dn':('3UTR','stop_codon') }
        
        # get exon boundaries from transcript and gene objects
        exon_bndrys = \
            [ (trans_index+1, gene.bp_bndry( exon_index )) for \
                  trans_index, exon_index in enumerate(trans) ]
        
        def add_cds_regions( orf_start, orf_stop, orf_index ):
            cds_regions = [ \
                (index, (max( exon_start, orf_start ), \
                             min( exon_stop, orf_stop ))) for \
                    index, (exon_start, exon_stop) in exon_bndrys \
                    if exon_stop >= orf_start and exon_start <= orf_stop ]
            
            # if strand is reverse iterate through cds regions in 
            # reverse in order to calculate frame correctly
            if reverse_strand:
                cds_regions.reverse()
            
            # first frame is always 0
            frame = 0
            for exon_trans_index, region in cds_regions:
                trans_orf_regions.append( \
                    ( region, None, exon_trans_index, gene, trans.name,'CDS', \
                          frame, orf_ids[orf_index], [orf_classes[orf_index],] ))
                
                # calculate the frame of the next cds region
                frame = ((( -(( region[1]-region[0]+1 )%3 ))%3 ) + frame )%3
            
            # return the exon number for the start and stop exons
            if reverse_strand:
                return cds_regions[-1][0], cds_regions[0][0]
            return cds_regions[0][0], cds_regions[-1][0]
        
        def add_up_strm_UTRs( orf_start, orf_stop, orf_index ):
            # get utr exons or portions of exons upstream of the orf
            up_strm_regions = [ \
                (index, (exon_start, min( exon_stop, orf_start - 1 ))) \
                    for index, (exon_start, exon_stop) in exon_bndrys \
                    if exon_start < orf_start ]
            
            # store utr regions upstream of cds
            for exon_trans_index, region in up_strm_regions:
                trans_orf_regions.append( \
                    ( region, None, exon_trans_index, gene, trans.name,\
                          feature_strings['up'][0], '.', orf_ids[orf_index], \
                          [orf_classes[orf_index],] ))
            
            return
        
        def add_up_strm_codon( up_strm_codon_start, exon_num, orf_index ):
            # get start or stop codon region depending on strand
            up_strm_codon_region = \
                (up_strm_codon_start, up_strm_codon_start+2) \
                if not reverse_strand \
                else (up_strm_codon_start-3, up_strm_codon_start-1)
            
            # add the codon region
            trans_orf_regions.append( \
                ( up_strm_codon_region, None, exon_num, gene, \
                      trans.name, feature_strings['up'][1], '0', \
                      orf_ids[orf_index], [orf_classes[orf_index],]))
            
            return
        
        def add_dn_strm_UTRs( orf_start, orf_stop, orf_index ):
            # get utr regions down stream of the orf
            dn_strm_regions = [ \
                (index, (max( exon_start, orf_stop+1 ), exon_stop)) \
                    for index, (exon_start, exon_stop) in exon_bndrys \
                    if exon_stop > orf_stop ]
            
            # store utr regions downstream of cds
            for exon_trans_index, region in dn_strm_regions:
                trans_orf_regions.append( \
                    ( region, None, exon_trans_index, gene, trans.name,\
                          feature_strings['dn'][0], '.', orf_ids[orf_index], \
                          [orf_classes[orf_index],] ))
            
            return
        
        def add_dn_strm_codon( dn_strm_codon_stop, exon_num, orf_index ):
            # get start or stop codon region depending on strand
            dn_strm_codon_region = \
                (dn_strm_codon_stop+1, dn_strm_codon_stop+3) \
                if not reverse_strand \
                else (dn_strm_codon_stop-2, dn_strm_codon_stop)
            
            # add the codon region
            trans_orf_regions.append( \
                ( dn_strm_codon_region, None, exon_num, gene,\
                      trans.name, feature_strings['dn'][1], '0', \
                      orf_ids[orf_index], [orf_classes[orf_index],]))
            
            return
        
        
        # loop through distinct ORFs and get all coorosponding regions
        for orf_index, (orf_start, orf_stop) in enumerate( orfs ):
            up_strm_exon_num, dn_strm_exon_num = \
                add_cds_regions( orf_start, orf_stop, orf_index )
            
            if print_UTRs:
                add_up_strm_UTRs( orf_start, orf_stop, orf_index )
                add_dn_strm_UTRs( orf_start, orf_stop, orf_index )
            if print_codons:
                add_up_strm_codon( orf_start, up_strm_exon_num, orf_index )
                add_dn_strm_codon( orf_stop, dn_strm_exon_num, orf_index )
        
        return trans_orf_regions
    
    def add_transcripts_with_orfs( self, all_orfs, print_non_orfs=True, \
                                       print_orf_exons=True, \
                                       print_codons=True, \
                                       print_UTRs=False ):
        """Write regions associated with orfs within genes and trans
        regions include: exon, CDS, 5UTR, 3UTR, start_codon and stop_codon
        Note CDS regions will always be printed all other regions have switches
        """
        for trans, gene, orfs, orf_classes, orf_ids in all_orfs:
            # initialize the list of regions to be added
            all_regions = []
            
            # if there are no orfs then optionally add the exons of the trans
            if len( orfs ) == 0:
                if print_non_orfs:
                    # store transcript exons if requested
                    for exon_trans_index, exon_gene_index in enumerate( trans ):
                        all_regions.append( \
                            (None, exon_gene_index, exon_trans_index+1, gene, \
                                 trans.name, 'exon', '.', None, orf_classes))
            
            # this transcript has an orf to be parsed
            else:
                if print_orf_exons:
                    # store transcript exons if requested
                    for exon_trans_index, exon_gene_index in enumerate( trans ):
                        all_regions.append( \
                            (None, exon_gene_index, exon_trans_index+1, gene, \
                                 trans.name, 'exon', '.', None, orf_classes))
                
                all_regions.extend( self._get_orf_regions( \
                        gene, trans, orfs, orf_classes, orf_ids, \
                            print_UTRs, print_codons ) )
            
            # sort regions by their position
            # Note that exons will not be sorted as their region is None
            # and will be determined later from gene and trans
            all_regions.sort( key = lambda x:x[0] )
            
            def iter_trans_lines():
                for region, gene_index, trans_index, gene, trans_name, \
                        feature_type, frame, orf_id, orf_classes \
                        in all_regions:
                    yield self._build_exon_gtf_line( \
                        gene_index, trans_index, gene, trans_name, \
                            region=region, feature_type=feature_type, \
                            frame=frame, orf_id=orf_id, \
                            orf_classes=orf_classes )
                
                return
            
            # write gtf lines for current transcript to file
            self.write_lock.acquire()
            for line in iter_trans_lines():
                self.write( line + "\n" )
            self.flush()
            self.write_lock.release()
        
        return

class GeneTracking( object ):
    """Store information from single gene for comparing expression across samples
    
    Emulates information provided by the cuffcompare tracking file
    Using gene_name for Cufflinks transfrag id, cufflinks locus id, and reference gene
    Would have to include other id values in initial parsing of exon lines
    """
    
    def __init__( self, gene, bam_fn ):
        self.gene = gene
        self.transcripts = {}
        
        return
    
    def _get_freq_max( self, transcripts ):
        """Get highest freq and assert that freqs total to 1
        """
        freq_max = max( transcripts.freqs )
        if not freq_max: freq_max = 1.0
        
        return freq_max
    
    def _get_trans_len_and_cov( self, trans, exon_cnts ):
        """Calculate the length and fragment coverage over a list of exons

        """
        length = 0
        cov = 0
        exon_ids = [ trans.up_strm_exons[0], ]
        exon_ids.extend( trans.internal_exons )
        exon_ids.append( trans.dn_strm_exons[-1] )
        
        for exon_index in exon_ids:
            start, stop = self.gene.exon_bndrys[ exon_index ]
            length += ( stop - start )
            cov += exon_cnts[ exon_index ]
                
        return length, cov
    
    def add_transcripts_from_single_experiment( self, transcripts, binned_reads ):
        """Add transcripts to GeneTracking object for one source

        """
        # get max of freq to determine FMI (fraction of major isoform)
        freq_max = self._get_freq_max( transcripts )

        exon_cnts = defaultdict(int)
        """
        # find the exon coverage for each exon this experiment, using the binned reads
        for read_grp_cnts in binned_reads.binned_reads.values():
            for bin, cnt in read_grp_cnts.iteritems():
                # take the first exon in each read
                exon_cnts[ bin[0] ] += 1
                exon_cnts[ bin[2] ] += 1
        """
        
        # add relevant information to a dictionary so that identical transcripts
        # from different bam files are added to the same list
        for trans_index, (trans, md) \
            in enumerate(transcripts.iter_transcripts_and_metadata()):
            length, cov = self._get_trans_len_and_cov( trans, exon_cnts )
            trans_id = self.gene.name + trans.get_str_identifier()
            
            # store a transcript entry for each source-transcript combination
            # fill entry only if transcript exists in that source
            try:
                self.transcripts[ trans_id ][ self.source_index[ md.sourcefile ] ] = \
                    ( md.sourcefile, trans_index, int(100*md.freq/freq_max), cov, length, trans )
            except KeyError:
                self.transcripts[ trans_id ] = \
                    [ None for i in xrange( self.num_sources ) ]
                self.transcripts[ trans_id ][ self.source_index[ md.sourcefile ] ] = \
                    ( md.sourcefile, trans_index, int(100*md.freq/freq_max), cov, length, trans )
        return

    def iter_ids_and_trans( self ):
        """iter transcript id and transcript objs in transcripts

        """
        for trans_id, source_list in self.transcripts.iteritems():
            for source_transcript in source_list:
                if source_transcript:
                    yield trans_id, source_transcript[5]
                    break
        return

class TrackingFile( file ):
    """Output tracking file object with function to write out each gene

    """
    def __init__( self, *args ):
        args = list( args )
        args.insert( 0, self )
        file.__init__( *args )

    def add_gene_tracking( self, gene_tracking ):
        """Write tracking information for a single gene to the tracking file

        """
        assert isinstance( gene_tracking, GeneTracking )
        
        gene_name = gene_tracking.gene.name
        transcript_lines = []
        for trans_id, sources_list in gene_tracking.transcripts.iteritems():
            # initialize tracking line with stub for:
            # transfrag_id, locus_id, reference_gene, trans_id, class_code
            transcript_line = "\t".join( ( gene_name, gene_name, gene_name, \
                                               trans_id, "-" ) )
            # fill tracking line with one field for each source ("-" for none)
            for source_index, source_transcript in enumerate( sources_list ) :
                if source_transcript:
                    source_index += 1
                    [ source, trans_index, FMI, cov, length, trans ] = source_transcript
                    source_field = "|".join( ( "\tq" + str( source_index ) + ":" + \
                                                   source, source + "." + \
                                                   str(trans_index), str(FMI), \
                                                   "0.000000",  "0.000000", \
                                                   "0.000000", str(cov), str(length) ) )
                else:
                    source_field = "\t-"
                transcript_line += source_field
            transcript_lines.append( transcript_line )

        if len( transcript_lines ):
            self.write( "\n".join( transcript_lines ) + "\n" )
            self.flush()

def iter_transcripts( exon_bndrys, junctions, start_exons, stop_exons ):
    """Iterate all transcripts.
    
    THIS ALWAYS WORKS IN THE DIRECTION OF GENOME COORDINATES.
    
    
    junctions - an iterable over tuples of connected exon indices
    start_exons - an iterable over exon indices that transcripts 
          are allowed to start in. ( ie, a tss list )
    
    Transcripts are the set of all exons 
    such that no 2 exons overlap. The algorithm 
    is recursive and as follows:
    
    Initialization:
        1) sort the exon list by starts
        2) add each exon to the list ( any exon could
           start the transcript )
    
    Recursive Step:
        We want to add any exons that could immediately 
        follow the last exon in the transcript. That is, 
        we add the first non-overlapping exon, and then
        any exons that overlap that exon.
    
    Final Step:
        There are no more exons to be added.
    
    """
    stop_exons = set( stop_exons )
    
    exons = exon_bndrys

    junctions_set = set( jn for jn in junctions if jn[0] != jn[1] )

    if VERBOSE:
        print >> sys.stderr, "{0:d} Exons and {1:d} Junctions".format( len(exons), len(junctions)  )
    
    Exon = namedtuple('Exon', ['start', 'stop'])

    exons = [ Exon(exon[0], exon[1]) for exon in  exons ]
    exons.sort( key=attrgetter('start') )
    
    if len( exons ) == 0:
        return
    
    def get_next_index( curr_index ):
        """Return the index of the next non-overlapping exon.
        
        """
        next_index = curr_index + 1
        while next_index < len( exons ) \
                and exons[next_index].start <= exons[curr_index].stop:
            next_index += 1
        return next_index
    
    def get_overlapping_exons( index ):
        """Get all exons that overlap exon 'index', including itself
        
        """
        next_index = index
        while next_index < len( exons ) \
                and exons[index].start <= exons[index].stop:
            yield next_index
            next_index += 1
        return
    
    # get the list of first exons
    # filter out start exons that appear in the equivalent first exon map
    # ( due to them being very close on the 5' end ) and sort
    new_start_exons = set()
    for exon in start_exons:
        new_start_exons.add( exon )
        
    start_exons = sorted(new_start_exons)
    transcripts = [ list((index,)) for index in start_exons ]
    
    while len( transcripts ) > 0:
        # get the next transcript off of the stack
        transcript = transcripts.pop()
        
        # make sure that this doesnt have an equivalent 3' end
        last_exon = transcript[-1]
        if last_exon in stop_exons:
            yield transcript
        
        # get the index of the next non-overlapping index
        next_nonoverlapping_index = get_next_index( transcript[-1] )
        if next_nonoverlapping_index > len( exons ):
            continue
        
        for next_exon_index in get_overlapping_exons( next_nonoverlapping_index ):
            # make sure the junction exists
            # if it does, then add this to the stack
            if ( transcript[-1], next_exon_index ) in junctions_set:
                transcript_copy = copy.copy( transcript )
                transcript_copy.append( next_exon_index )
                transcripts.append( transcript_copy )
    
    return        


def find_UTRs( exon_bndrys, connected_exon_pairs, \
                   prior_start_exons=(), prior_stop_exons=() ):
    """Find the ( potential ) UTR regions. 
    
    These are always in order of genomic coordinates ( ie + stranded genes )
    """
    if VERBOSE: print "Prior Start Exons:", prior_start_exons
    if VERBOSE: print "Prior Stop Exons:", prior_stop_exons
    
    # if the caller doesnt provide a list of start exons, then we assume that 
    # any transcript can start anywhere
    jn_starts = set( s for f, s in connected_exon_pairs )
    jn_ends = set( f for f, s in connected_exon_pairs )


    # initialize the start exon sets
    start_exon_indices = set( prior_start_exons )
    stop_exon_indices = set( prior_stop_exons )
    
    if USE_ALL_POSSIBLE_STARTS:
        start_exon_indices.update( xrange(len(exon_bndrys)) )
    elif FILTER_SPLICED_UTR_EXONS:
        # find all of the exons that don't have a junction that *ends* with them.
        # Such exons must be start locations, if they exist at all.
        for exon_i in xrange(len(exon_bndrys)):
            if exon_i not in jn_starts:
                start_exon_indices.add( exon_i )
            if exon_i not in jn_ends:
                stop_exon_indices.add( exon_i )
    # otherwise, just use the passed list of start exons
    else:
        assert len( prior_start_exons ) > 0
        pass
    
    if VERBOSE: 
        print "upstream UTRs:", start_exon_indices
        print "downstream UTRs :", stop_exon_indices

    return start_exon_indices, stop_exon_indices
    
def find_connected_exons( binned_reads, prior_connected_exons=() ):
    ## Find all of the valid transcripts
    if USE_PROVIDED_JNS: connected_exons = sorted( prior_connected_exons )
    else: connected_exons = list()
    
    if VERBOSE: print "Prior Connected Exons:", connected_exons

    if USE_DISCOVERED_JNS:
        connected_exons.extend( copy.deepcopy( binned_reads.connected_exon_pairs ) )    
        if VERBOSE: print "Discovered Connected Exons:", sorted(binned_reads.connected_exon_pairs)
    
    connected_exons= sorted( set( connected_exons ) )
    if VERBOSE: print "Connected Exons:", connected_exons
    if VERBOSE: print
    
    return connected_exons

def build_possible_transcripts( gene, binned_reads, prior_connected_exons=[], \
                                    start_exons=None ):
    """Build all the possible transcripts, using several heuristics.
    
    gene: a gene object
    binned_reads: a bined_reads object
    prior_connected_exons: an iterable over tuples of connected exon indices
    start_exons: an iterable of exons that transcripts are allowed to start in
    
    TODO - document the heuristics here.
    """
    # we only consider exons with less than MAX_NUM_EXONS exons
    if len( gene.exon_bndrys ) > sparsify_transcripts.py.MAX_NUM_EXONS:
        if VERBOSE:
            print "TOO MANY EXONS...."
        raise sparsify_transcripts.py.GeneProcessingError( \
            gene, binned_reads.reads, \
                "Too many exons ( {0:d} )".format(len( gene.exon_bndrys ) ) )
    
    # check to make sure we have at least 1 exon 
    if len( gene.exon_bndrys ) == 0:
        if PAUSE_ON_ERROR:
            print( "BUG PROBABLY - can't fit gene with 0 exons." )
            raw_input("Press enter to conitnue...")
        return Transcripts( gene )

    # Make sure that we have enough bins to continue
    if len( binned_reads.binned_reads.keys() ) == 0:
        if VERBOSE:
            print "Cannot estimate transcript frequencies for genes without valid bins."
        raise sparsify_transcripts.py.GeneProcessingError( \
            gene, binned_reads.reads, "No observed reads at locus." )
    
    # Find all of the exons that have been observed to be adjoint in a transcript.
    connected_exons = find_connected_exons( binned_reads, prior_connected_exons )
    
    ## Find the possible transcript starts and ends
    start_exons, stop_exons = find_UTRs( gene.exon_bndrys, connected_exons )

    raw_trans = iter_transcripts( \
        gene.exon_bndrys, connected_exons, start_exons, stop_exons )
    
    # group transcripts that have the same internal structure
    grouped_transcripts = {}
    for entry in raw_trans:
        if VERBOSE:
            print "RAW TRANSCRIPT:", entry
        key = ( gene.bp_bndry( entry[0] )[1], tuple( entry[1:-1] ), \
                    gene.bp_bndry( entry[-1] )[0] )
        if not grouped_transcripts.has_key( key ):
            grouped_transcripts[ key ] = ( set(), set()  )
        grouped_transcripts[ key ][0].add( entry[0] )
        grouped_transcripts[ key ][1].add( entry[-1] )

    transcripts = Transcripts( gene )
    for key, ( starts, stops ) in grouped_transcripts.iteritems():
        transcripts.add_transcript( \
            Transcript( starts, stops, key[1] ),  binned_reads.sourcefile )
    
    if VERBOSE:
        print( "Valid transcripts: " + gene.strand )
        print transcripts.pretty_str(False)
    
    return transcripts
