import sys
from genomic_intervals import GenomicInterval
from collections import namedtuple
import itertools

GffLine = namedtuple( "GffLine", ["region", "feature", "score", 
                                  "source", "frame", "group"]  )
GtfLine = namedtuple( "GtfLine", ["region", "gene_id", "trans_id", "feature", 
                                  "score", "source", "frame", "meta_data"])

def parse_gff_line( line, fix_chrm=True ):
    data = line.split()
    if len( data ) < 9: 
        return None
    
    # fix the chromosome if necessary
    if fix_chrm and data[0].startswith("chr"):
        data[0] = data[0][3:]
    
    # the source is required and always the 2nd entry
    # source = data[1]
    
    # the type of element - ie exon, transcript, etc.
    # its required
    # feature_type = data[2]
    
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
        data[5] = 0
    else:
        try:
            data[5] = int( data[5] )
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

def parse_gtf_line( line, fix_chrm=True ):
    gffl = parse_gff_line( line, fix_chrm=fix_chrm )
    if gffl == None: return None
    
    meta_data_items = (gffl.group).split()
    # parse the meta data, and grab the gene name
    meta_data = dict( zip( meta_data_items[::2], meta_data_items[1::2] ) )
    
    if "gene_id" not in meta_data:
        raise ValueError, "GTF lines require a gene_id field."
    if "transcript_id" not in meta_data:
        raise ValueError, "GTF lines require a transcript_id field."
    
    # get gene and transcript name if parsing a gtf line 
    # else it is a gff line and does not have gene or trans names
    def get_name_from_field( name ):
        if name.startswith('"'):
            name = name[1:]
        if name.endswith('";'):
            name = name[:-2]
        if name.endswith('"'):
            name = name[:-1]
        return name
    
    gene_name = meta_data[ 'gene_id' ]
    gene_name = get_name_from_field( gene_name )
    
    trans_name = meta_data[ 'transcript_id' ]
    trans_name = get_name_from_field( trans_name )
    
    return GtfLine( gffl[0], gene_name, trans_name, \
                        gffl[1], gffl[2], gffl[3], gffl[4], meta_data )

def create_gff_line( region, grp_id, score=0, \
                         feature='.', source='.', frame='.' ):
    r = region
    chrm = region.chr
    if not chrm.startswith( 'chr' ):
        chrm = 'chr' + chrm
    
    gff_line = [ chrm, source, feature, str(r.start), str(r.stop), \
                 str(score), r.strand, frame, str( grp_id ) ]
    return '\t'.join( gff_line )

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
