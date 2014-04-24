# Copyright (c) 2011-2012 Nathan Boley

import os, sys
from collections import namedtuple, defaultdict
import itertools

from ..transcript import Gene, Transcript, GenomicInterval
    
from reads import clean_chr_name

GffLine = namedtuple( "GffLine", ["region", "feature", "score", 
                                  "source", "frame", "group"]  )
GtfLine = namedtuple( "GtfLine", ["region", "gene_id", "trans_id", "feature", 
                                  "score", "source", "frame", "meta_data"])

VERBOSE = True
DEBUG = False
    
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
    fpk = None
    fpkm = None
    conf_lo = None
    conf_hi = None
    frac = None
    
    promoter = None
    polya_region = None
    CDS_start, CDS_stop = None, None
    for line in transcript_lines:
        scores.add( line.score )
        if fpk == None and 'fpk' in line.meta_data: 
            fpk = float(line.meta_data['FPK'])
        if fpkm == None and 'FPKM' in line.meta_data: 
            fpkm = float(line.meta_data['FPKM'])
        if fpkm == None and 'fpkm' in line.meta_data: 
            fpkm = float(line.meta_data['fpkm'])
        if conf_lo == None and 'conf_lo' in line.meta_data: 
            conf_lo = float(line.meta_data['conf_lo'])
        if conf_hi == None and 'conf_hi' in line.meta_data: 
            conf_hi = float(line.meta_data['conf_hi'])
        if frac == None and 'frac' in line.meta_data: 
            frac = float(line.meta_data['frac'])
        
        # subtract one to make the coordinates 0 based
        line_start, line_stop = line.region.start-1, line.region.stop-1
        if line.feature == 'exon':
            exons.append( (line_start, line_stop) )
        elif line.feature == 'CDS':
            exons.append( (line_start, line_stop) )
            CDS_start = line_start \
                if CDS_start == None or line_start < CDS_start else CDS_start
            CDS_stop = line_stop \
                if CDS_stop == None or line_stop > CDS_stop else CDS_stop
        elif line.feature == 'promoter':
            promoter = (line_start, line_stop)
        elif line.feature == 'polya':
            polya_region = (line_start, line_stop)
    
    if len( exons ) == 0:
        return None
    
    CDS_region = None if CDS_start == None or CDS_stop == None \
        else (CDS_start, CDS_stop)
    
    score = next(iter(scores)) if len(scores) == 1 else None
    line = transcript_lines[0]
    return Transcript( line.trans_id, line.region.chr, line.region.strand, 
                       flatten(sorted(exons)), CDS_region,
                       gene_id=line.gene_id, score=score, fpkm=fpkm, fpk=fpk, 
                       promoter=promoter, polya_region=polya_region,
                       conf_lo=conf_lo, conf_hi=conf_hi, frac=frac)

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
    elif len(gene_lines) == 1:
        gene_data = gene_lines[0]
        gene_chrm, gene_strand, gene_start, gene_stop = gene_data.region
        # since gff's are 1 based
        gene_start -= 1
        gene_stop -= 1
        gene_meta_data = gene_data.meta_data
    else:
        if VERBOSE: print >> sys.stderr, "Skipping '%s': multiple gene lines." % gene_id
        return None
        
    if gene_start != min(t.start for t in transcripts ) \
            or gene_stop != max(t.stop for t in transcripts ):
        if VERBOSE: print >> sys.stderr, "Skipping '%s': gene boundaries dont match the transcript boundaries." % gene_id
        return None
    
    return Gene( gene_id, gene_chrm, gene_strand, 
                 gene_start, gene_stop,
                 transcripts, gene_meta_data )
    

def load_gtf(fname_or_fp, use_name_instead_of_id=False, 
             contig=None, strand=None):
    if isinstance( fname_or_fp, str ):
        fp = open( fname_or_fp )
    else:
        assert isinstance( fname_or_fp, file )
        fp = fname_or_fp
    
    gene_lines = defaultdict(lambda: ( defaultdict(list), [] ))
    for line in fp:
        data = parse_gtf_line(line, fix_chrm=True, 
                              use_name_instead_of_id=use_name_instead_of_id)
        # skip unparseable lines, or lines without gene ids
        if None == data: continue
        if data.gene_id == "": continue
        if contig != None and data.region.chr != contig: continue
        if strand != None and data.region.strand != strand: continue
        # add gene lines directly to the gene object
        if data.feature == 'gene': 
            gene_lines[data.gene_id][1].append( data )
        else:
            gene_lines[data.gene_id][0][data.trans_id].append(data)
    
    genes = []
    for gene_id, ( transcripts_data, gene_lines ) in gene_lines.iteritems():
        try:
            gene = _load_gene_from_gtf_lines(
                gene_id, gene_lines, transcripts_data)
        except Exception, inst:
            print >> sys.stderr, "ERROR : Could not load '%s'" % gene_id
            print >> sys.stderr, "DETAIL: %s" % str( inst )
            if DEBUG: raise
            gene = None
        
        if gene == None: continue
        genes.append( gene )
    
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
    
    # if gene and trans id are included in meta data, then remove them so that
    # we dont build them twice
    meta_data_str = ( 'gene_id "%s"; transcript_id "%s";' 
        % ( gene_id, transcript_id ) )
    # an optimization for the common case that only gene and transcript id
    # are included
    if len( meta_data ) > 0:
        # remove geneid and transctipt id from the dict if they
        # exist so that we dont double count them
        try: del meta_data['gene_id']
        except KeyError: pass
        try: del meta_data['transcript_id']
        except KeyError: pass
        
        items = [ " " + str(k) + ' "%s";' % str(v) 
                     for k, v in meta_data.iteritems() ]
        meta_data_str += "".join( items )
    
    # add one because gtf lines are 1 based
    gtf_line = [ r.chr, source, feature, str(r.start+1), str(r.stop+1), 
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

def iter_gff_lines( regions_iter, grp_id_iter=None, score='.', 
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
