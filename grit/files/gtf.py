# Copyright (c) 2011-2012 Nathan Boley

import os, sys
from collections import namedtuple, defaultdict
import itertools

GenomicInterval = namedtuple('GenomicInterval', 
                             ['chr', 'strand', 'start', 'stop'])
GffLine = namedtuple( "GffLine", ["region", "feature", "score", 
                                  "source", "frame", "group"]  )
GtfLine = namedtuple( "GtfLine", ["region", "gene_id", "trans_id", "feature", 
                                  "score", "source", "frame", "meta_data"])

VERBOSE = True
DEBUG = False

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


class Transcript( object ):
    def __init__(self, trans_id, chrm, strand, exons, cds_region,
                 gene_id=None, score=None, rpkm=None, rpk=None, promoter=None ):
        self.gene_id = gene_id
        self.id = trans_id
        self.chrm = chrm
        self.strand = strand
        
        self.score = score
        self.rpkm = rpkm
        self.rpk = rpk
        
        exon_bnds = list( itertools.chain( *exons ) )
        self.exon_bnds = exon_bnds
        self.start = exon_bnds[0]
        self.stop = exon_bnds[-1]
        assert self.start <= self.stop

        self.exons = tuple(zip(exon_bnds[:-1:2], exon_bnds[1::2]))
        self.introns = tuple([ (x+1, y-1) for x, y in 
                               itertools.izip(exon_bnds[1:-2:2], 
                                              exon_bnds[2:-1:2]) ])
        
        self.is_protein_coding = ( cds_region != None )
        
        self.cds_region = cds_region
        self.start_codon = None
        self.stop_codon = None
        
        self.cds_exons = None
        self.fp_utr_exons = None
        self.tp_utr_exons = None
        self.us_exons = None
        self.ds_exons = None
        
        self.promoter = promoter
        if cds_region != None:
            self.add_cds_region( cds_region )
            
    def add_cds_region( self, cds_region ):
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
                in itertools.izip( self.introns, self.exons ):
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
    
    def build_gtf_lines( self, gene_id, meta_data, source='.'):
        ret_lines = []
        def build_lines_for_feature( exons, feature, is_CDS=False ):
            current_frame = 0
            score = str(self.score) if self.score != None else '.'
            for start, stop in exons:
                region = GenomicInterval( self.chrm, self.strand, start, stop )
                frame = current_frame if is_CDS else '.'
                yield create_gtf_line( region, gene_id, self.id, meta_data,
                                       score, feature=feature, frame=str(frame),
                                       source=source)
                current_frame = ( current_frame + stop - start + 1 )%3
            return
        
        ret_lines.extend( build_lines_for_feature( 
                self.exons, 'exon', False ) )
        
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
        
        return "\n".join( ret_lines )

def flatten( regions ):
    regions.sort()
    new_regions = []
    curr_start = regions[0][0]
    curr_end = regions[0][1]
    for i,(start,end) in enumerate(regions):
        if curr_end > end:
            end = curr_end
        if i+1 == len( regions ): 
            if len(new_regions) == 0:
                new_regions = [curr_start, curr_end]
                break
            if new_regions[-1][1] == end:
                break
            else:
                new_regions.append( [ curr_start, curr_end] )
                break
        if regions[i+1][0]-end <= 1:
            curr_end = max( regions[i+1][1], end ) 
        else:
            new_regions.append( [curr_start, curr_end] )
            curr_start = regions[i+1][0]
            curr_end = regions[i+1][1]
    if type(new_regions[0]) == int:
        return [new_regions]
    else:
        return new_regions
    
class Gene( list ):
    def __init__(self, id,chrm, strand, start, stop, transcripts, meta_data={}):
        self.id = id
        self.chrm = chrm
        self.strand = strand
        self.start = start
        self.stop = stop
        self.transcripts = transcripts
        self.meta_data = meta_data
        list.extend( self, ( id, chrm, strand, start, stop, transcripts ) )
        return    

    def find_transcribed_regions( self ):
        exons = set()
        for transcript in self.transcripts:
            exons.update( transcript.exons )
        
        return flatten( sorted( exons ) )
    
    def calc_bpkm(self, read_cov):
        base_cov, length = 0.0, 0
        for start, stop in self.find_transcribed_regions():
            length += stop - start + 1
            base_cov += read_cov[(self.chrm, self.strand)][start:stop+1].sum()
        
        return base_cov/length
        
    def extract_elements(self):
        """Extract the different element types.
        
        """
        for t in gene.transcripts:
            if len( t.exons ) == 1:
                single_exon_genes.add( ( t.chrm, t.strand, t.exons[0][0], t.exons[0][1] ) )
            else:
                tss_exon = t.exons[0] if t.strand == '+' else t.exons[-1]
                tss_exons.add( (t.chrm, t.strand, tss_exon[0], tss_exon[1]) )
                internal_exons.update( (t.chrm, t.strand, x1, x2) for x1, x2 in t.exons[1:-1] )
                introns.update( (t.chrm, t.strand, x1, x2) for x1, x2 in t.introns )
                tes_exon = t.exons[-1] if t.strand == '+' else t.exons[0] 
                tes_exons.add( (t.chrm, t.strand, tes_exon[0], tes_exon[1]) )



def parse_gff_line( line, fix_chrm=True ):
    if line.startswith( '#' ): return None
    
    data = line.split()
    if len( data ) < 9: 
        return None
    
    # fix the chromosome if necessary
    if fix_chrm and data[0].startswith("chr"):
        data[0] = data[0][3:]
    
    # the source is required and always the 2nd entry
    source = data[1]
    
    # the type of element - ie exon, transcript, etc.
    # its required
    feature_type = data[2]
    
    # check that this is a valid gff/gtf line else return None
    # check that start and stop are valid integers, and convert
    # them to integers
    try:
        data[3] = int( data[3] )
        data[4] = int( data[4] )
    except ValueError:
        return None
    
    # check the score
    if data[5] == '.':
        data[5] = '.'
    else:
        try:
            data[5] = float(data[5])
        except ValueError:
            return None
    
    # check for valid strand
    if data[6] not in "+-.":
        return None
    
    # check that the frame is valid
    if data[7] not in ( '0', '1', '2', '.' ):
        return None
    if data[7] in ( '0', '1', '2' ):
        data[7] = int( data[7] )
    
    data[8] = " ".join( data[8:] )
    
    return GffLine( GenomicInterval(data[0], data[6], data[3], data[4]), \
                        data[2], data[5], data[1], data[7], data[8] )

def parse_gtf_line( line, fix_chrm=True, use_name_instead_of_id=False ):
    gffl = parse_gff_line( line, fix_chrm=fix_chrm )
    if gffl == None: return None
    
    # get gene and transcript name if parsing a gtf line 
    # else it is a gff line and does not have gene or trans names
    def get_name_from_field( name ):
        if name.startswith('"'):
            name = name[1:]
        if name.endswith(';'):
            name = name[:-1]
        if name.endswith('"'):
            name = name[:-1]
        return name
    
    meta_data_items = (gffl.group).split()
    # parse the meta data, and grab the gene name
    meta_data = dict( zip( meta_data_items[::2], 
                           ( get_name_from_field(x) 
                             for x in meta_data_items[1::2] ) ) )
    
    if "gene_id" not in meta_data:
        raise ValueError, "GTF lines require a gene_id field."
    if "transcript_id" not in meta_data:
        raise ValueError, "GTF lines require a transcript_id field."
        
    if use_name_instead_of_id and 'gene_name' in meta_data:
        gene_name = meta_data[ 'gene_name' ]
    else:
        gene_name = meta_data[ 'gene_id' ]
    gene_name = get_name_from_field( gene_name )

    if use_name_instead_of_id and 'transcript_name' in meta_data:
        trans_name = meta_data[ 'transcript_name' ]
    else:
        trans_name = meta_data[ 'transcript_id' ]    
    trans_name = get_name_from_field( trans_name )
    
    return GtfLine( gffl[0], gene_name, trans_name,
                    gffl[1], gffl[2], gffl[3], gffl[4], meta_data )

def load_transcript_from_gtf_data(transcript_lines):
    exons = []
    scores = set()
    rpk = None
    rpkm = None
    promoter = None
    CDS_start, CDS_stop = None, None
    for line in transcript_lines:
        scores.add( line.score )
        if rpk == None and 'rpk' in line.meta_data: 
            rpk = line.meta_data['rpk']
        if rpkm == None and 'rpkm' in line.meta_data: 
            rpkm = line.meta_data['rpkm']
        
        if line.feature == 'exon':
            exons.append( (line.region.start, line.region.stop) )
        elif line.feature == 'CDS':
            exons.append( (line.region.start, line.region.stop) )
            CDS_start = line.region.start \
                if CDS_start == None or line.region.start < CDS_start else CDS_start
            CDS_stop = line.region.stop \
                if CDS_stop == None or line.region.stop > CDS_stop else CDS_stop
        elif line.feature == 'promoter':
            promoter = line.region
    
    if len( exons ) == 0:
        return None
    
    CDS_region = None if CDS_start == None or CDS_stop == None \
        else (CDS_start, CDS_stop)
    
    score = next(iter(scores)) if len(scores) == 1 else None
    line = transcript_lines[0]
    return Transcript( line.trans_id, line.region.chr, line.region.strand, 
                       flatten(sorted(exons)), CDS_region,
                       line.gene_id, score, rpkm, rpk, promoter )

def _load_gene_from_gtf_lines( gene_id, gene_lines, transcripts_data ):
    if len( gene_lines ) > 1:
        raise ValueError, "Multiple gene lines for '%s'" % gene_id
    
    transcripts = []
    for trans_id, transcript_lines in transcripts_data.iteritems():
        transcript = load_transcript_from_gtf_data(transcript_lines)
        if transcript == None: continue
        transcripts.append(transcript)

    # if there are no gene lines, then get the info from the transcripts
    if len( gene_lines ) == 0:
        gene_chrm, gene_strand = transcripts[0].chrm, transcripts[0].strand
        gene_start = min(t.start for t in transcripts )
        gene_stop = max(t.stop for t in transcripts )
        gene_meta_data = {}
    else:
        gene_data = gene_lines[0]
        gene_chrm, gene_strand, gene_start, gene_stop = gene_data.region
        gene_meta_data = gene_data.meta_data

    if gene_start != min(t.start for t in transcripts ) \
            or gene_stop != max(t.stop for t in transcripts ):
        if VERBOSE: print >> sys.stderr, "Skipping '%s': gene boundaries dont match the transcript boundaries." % gene_id
        return None
    
    return Gene( gene_id, gene_chrm, gene_strand, 
                 gene_start, gene_stop,
                 transcripts, gene_meta_data )
    

def load_gtf(fname_or_fp, use_name_instead_of_id=False):
    if isinstance( fname_or_fp, str ):
        fp = open( fname_or_fp )
    else:
        assert isinstance( fname_or_fp, file )
        fp = fname_or_fp
    
    gene_lines = defaultdict(lambda: ( defaultdict(list), [] ))
    for line in fp:
        data = parse_gtf_line(line, fix_chrm=True, 
                              use_name_instead_of_id=use_name_instead_of_id)
        if None == data: continue
        # add gene lines directly to the gene object
        if data.feature == 'gene': 
            gene_lines[data.gene_id][1].append( data )
        else:
            gene_lines[data.gene_id][0][data.trans_id].append(data)
    if VERBOSE: print >> sys.stderr, "Finished parsing gene lines."
    
    genes = []
    for gene_id, ( transcripts_data, gene_lines ) in gene_lines.iteritems():
        try:
            gene = _load_gene_from_gtf_lines(gene_id, gene_lines, transcripts_data)
        except Exception, inst:
            print >> sys.stderr, "ERROR : Could not load '%s'" % gene_id
            print >> sys.stderr, "DETAIL: %s" % str( inst )
            if DEBUG: raise
            gene = None
        
        if gene == None: continue
        genes.append( gene )

    if VERBOSE: print >> sys.stderr, "Finished building gene objects."
    
    if isinstance( fname_or_fp, str ):
        fp.close()
    
    return genes

def create_gff_line( region, grp_id, score=0, \
                     feature='.', source='.', frame='.' ):
    r = region
    chrm = region.chr
    if not chrm.startswith( 'chr' ):
        chrm = 'chr' + chrm
    
    gff_line = [ chrm, source, feature, str(r.start), str(r.stop), \
                 str(score), r.strand, frame, str( grp_id ) ]
    return '\t'.join( gff_line )

def create_gff3_line( chrm, strand, start, stop, type,
                      ID, name, parents=[],
                      source='.', frame='.', score='.' ):
    parents_str = "Parent=" + ",".join(parents)+";" if len(parents) > 0 else ""
    meta_data = "ID=%s;%sName=%s;" % ( ID, parents_str, name )
    return "\t".join( map( str, [
                'chr' + chrm, source, type, start, stop, score, 
                strand, stop, frame, meta_data
                ] ) )

def create_gtf_line( region, gene_id, transcript_id, meta_data, score=0, \
                         feature='.', source='.', frame='.' ):
    r = region
    chrm = region.chr
    if not chrm.startswith( 'chr' ):
        chrm = 'chr' + chrm
    
    # if gene and trans id are included in meta data, then remove them so that
    # we dont build them twice
    meta_data_str = 'gene_id "%s"; transcript_id "%s";' \
        % ( gene_id, transcript_id )
    # an optimization for the common case that only gene and transcript id
    # are included
    if len( meta_data ) > 0:
        # remove geneid and transctipt id from the dict if they
        # exist so that we dont double count them
        try: del meta_data['gene_id']
        except KeyError: pass
        try: del meta_data['transcript_id']
        except KeyError: pass
        
        items = [ " " + str(k) + ' "%s";' % str(v) \
                     for k, v in meta_data.iteritems() ]
        meta_data_str += "".join( items )
    
    gtf_line = [ chrm, source, feature, str(r.start), str(r.stop), \
                 str(score), r.strand, frame, meta_data_str ]
    return '\t'.join( gtf_line )

def build_iters( items ):
    iters = []
    for item in items:
        # if it's a string assume we want it to be the same on every line
        if isinstance( item, (str, int) ):
            iters.append( itertools.repeat( str(item) ) )
        # otherwise, assume it's an iterable
        else:
            iters.append( item )
    
    return iters 

def iter_gff_lines( regions_iter, grp_id_iter=None, score='.', \
                         feature='.', source='.', frame='.' ):
    if isinstance( grp_id_iter, str ):
        raise ValueError, "Group ID must be an iterator."
    
    if grp_id_iter == None:
        grp_id_iter = ( str(i) for i in itertools.count(1) )
    
    iters = [ regions_iter, grp_id_iter ] \
        + build_iters(( score, feature, source, frame ))
    
    for data in itertools.izip( *iters ):
        yield create_gff_line( *data )
    
    return

def iter_gtf_lines( regions_iter, gene_id_iter, trans_id_iter,               \
                    grp_id='.', score='.', feature='.', source='.', frame='.', \
                    meta_data_iter=None ):
    # if we didnt provide any additional meta data, then make an empty list iter
    if meta_data_iter == None:
        meta_data_iter = itertools.repeat( [] )
    # build the data iterators, using repeat for defaults
    iters = [ regions_iter, gene_id_iter, trans_id_iter, meta_data_iter ] \
          +  build_iters(( score, feature, source, frame ))
    for data in itertools.izip( *iters ):
        yield create_gtf_line( *data )
    
    return

if __name__ == '__main__':
    with open( sys.argv[1] ) as fp:
        for gene in load_gtf( fp ):
            for t in gene.transcripts:
                print t.build_gtf_lines( gene.id, {} )
