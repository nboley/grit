import sys, os

from collections import defaultdict, OrderedDict
from itertools import izip, repeat
sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../file_types/fast_gtf_parser/" ) )
from gtf import load_gtf

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../file_types/" ) )
from gtf_file import iter_gtf_lines
from genomic_intervals import GenomicInterval

from compare_annotations import cluster_overlapping_genes

VERBOSE = True


def build_gtf_lines( chrm, strand, \
                     gene_id, trans_id, \
                     exon_bndries, \
                     sources_str=None ):
    feature="exon"
    source="slimerge"
    gene_name = "merged_cluster_%i" % gene_id
    gene_name_iter = repeat( gene_name )
    trans_name_iter = repeat( gene_name + "_%i" % trans_id )
    # build the soruces meta data, if necessary
    meta_data_iter = None
    if sources_str != None:
        meta_data_iter = repeat( \
            OrderedDict((("type", "mrna"), ("sources", sources_str))) )
    
    exons_iter = enumerate( izip( exon_bndries[0::2], exon_bndries[1::2] ) )
    regions_iter = ( GenomicInterval( chrm, strand, start, stop ) \
                         for exon_num, (start, stop) in exons_iter )
    
    return iter_gtf_lines( regions_iter, gene_name_iter, trans_name_iter, \
                               feature=feature, source=source, \
                               meta_data_iter = meta_data_iter )


def cluster_transcripts( ):
    # first, cluster the exons
    pass

def parse_arguments():
    import argparse
    desc = 'Merge transcripts.'
    parser = argparse.ArgumentParser(description=desc)    
    parser.add_argument( "gtfs", type=str, nargs="+", \
        help='GTF file contraining transcripts. (i.e. from build_transcripts)' )
    
    parser.add_argument( '--out-fname', '-o', \
        help='Output file. default: stdout')

    parser.add_argument( '--verbose', '-v', \
                         default=False, action='store_true', \
        help='Whether or not to print status information.')

    parser.add_argument( '--produce-sources', \
                         default=False, action='store_true', \
        help='Whether or not to add source information for each transcript.')
                         
    args = parser.parse_args()
    
    out_fp = open( args.out_fname ) if args.out_fname != None else sys.stdout
    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.gtfs, out_fp, args.produce_sources

def cluster_overlapping_exons(exons):
    exons.sort()
    exon_groups = []
    grp_stop = -1
    current_grp = []
    for start, stop in exons:
        # if start is > the current group stop add the current group and 
        # start a new group
        if start > grp_stop and current_grp:
            exon_groups.append( current_grp )
            current_grp = []
        # if this stop extends the group stop do so
        if stop > grp_stop:
            grp_stop = stop
        # append this exon to the current group
        current_grp.append( (start, stop) )

    # if exons is not empty and not exons were added then add the last group
    if current_grp:
        exon_groups.append( current_grp )        
    
    return exon_groups

def main():
    gtf_fnames, ofp, produce_sources = parse_arguments()
    
    # parse and load all of the gtf files
    transcriptomes = []
    for gtf_fname in gtf_fnames:
        try:
            gtf_data = load_gtf( gtf_fname )
        except IOError:
            print >> sys.stderr, \
                "Warning: couldn't extra any data from ", gtf_fname
            transcriptomes.append( [] )
            continue
        
        transcriptomes.append( gtf_data )
        if VERBOSE: print >> sys.stderr, "Finished parsing", gtf_fname
    
    
    # build the sets of transcripts
    transcript_maps = defaultdict( list )
    for source, transcriptome in izip( gtf_fnames, transcriptomes ):
        for cluster_id, chrm, strand, start, stop, transcripts in transcriptome:
            for trans_id, exons in transcripts:
                key = ( chrm, strand, tuple( exons ) )
                transcript_maps[ key ].append( source )
        """
    # cluster the transcripts:
    grpd_exons = defaultdict( set )
    for trans_id, (chrm, strand, exons) in \
            enumerate( sorted( transcript_maps ) ):
        grpd_exons[(chrm, strand)].update( izip( exons[::2], exons[1::2] ) )
    
    # sort and cluster the exons
    for exons in grpd_exons.itervalues():
        clustered_exons = cluster_overlapping_exons( sorted(exons) )
        for cluster in clustered_exons:
            print cluster
    
    return
    """
    # iterate through unique transcripts, and write them to the output file
    for trans_id, (chrm, strand, exons) in \
            enumerate( sorted( transcript_maps ) ):
        if produce_sources:
            source_str = ",".join( transcript_maps[ (chrm, strand, exons) ] )
        else:
            source_str = None
        
        gtf_lines = build_gtf_lines( \
            chrm, strand, trans_id, 0, exons, source_str )
        
        ofp.write( "\n".join( gtf_lines ) + "\n" )
    
    return
    
if __name__ == "__main__":
    main()
