# Copyright (c) 2011-2012 Nathan Boley

import sys, os
from itertools import izip
from collections import defaultdict, namedtuple

sys.path.append( os.path.join( os.path.dirname( __file__), 
                               "../file_types/fast_gtf_parser/" ) )
from gtf import load_gtfs, Transcript

sys.path.append( os.path.join( os.path.dirname( __file__), 
                               "../elements/transcripts/" ) )
from merge_transcripts import cluster_transcripts

sys.path.append( os.path.join( os.path.dirname( __file__), 
                               "../file_types/" ) )
from interval import binary_region, interval

VERBOSE = False

def match_clustered_transcripts( transcripts, distance_fn ):
    ref_transcripts = [ trans for trans, source in transcripts
                        if source == 'ref' ]
    input_transcripts = [ trans for trans, source in transcripts
                          if source == 'input' ]

    # if there are no reference transcripts in this cluster, then 
    # just yield the original and a None to indicate this
    if len( ref_transcripts ) == 0:
        return zip( input_transcripts, [None]*len(input_transcripts) )
    
    matches = []
    for i_trans in input_transcripts:
        best_match = ref_transcripts[0]
        best_score = distance_fn( i_trans, ref_transcripts[0] )
        for ref_trans in ref_transcripts[1:]:
            new_score = distance_fn( i_trans, ref_trans )
            if new_score < best_score:
                best_match = ref_trans
                best_score = new_score
        
        matches.append( ( i_trans, best_match ) )
    
    return matches

class HexavijezimalCntr( object ):
    def __init__( self, upper_case=False ):
        self.value = 2
        self._upper_case = upper_case
    
    _char_mapping = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.lower()
    
    def __str__( self ):
        return_chars = []
        assert self.value > 0
        tmp_val = self.value-1
        
        while tmp_val > 0:
            remainder = (tmp_val-1)%26
            return_chars.append( self._char_mapping[remainder] )
            tmp_val =  (tmp_val + 1 - remainder)/26
        
        if self._upper_case:
            return "".join( reversed(return_chars) ).upper()
        else:
            return "".join( reversed(return_chars) )
    
    def increment( self ):
        self.value += 1

def base_overlap( t1, t2 ):
    if max(t1.exon_bnds) < min(t2.exon_bnds): return 0
    if min(t1.exon_bnds) > max(t2.exon_bnds): return 0
    t1r = binary_region(interval(*x) for x in t1.exons)
    t2r = binary_region(interval(*x) for x in t2.exons)
    return float(t1r.overlap(t2r))/max(t1r.featuresLength(), t2r.featuresLength())

def exact_match_distance( i_trans, ref_trans ):
    return int( i_trans.IB_key() == ref_trans.IB_key() )

def lossy_distance( i_trans, ref_trans ):
    # check to see if the internal bounaries match exactly
    CDS_match = 1-int( i_trans.IB_key() == ref_trans.IB_key() )
    # check to see if the internal exon boundaries match
    UB_match = 1-int( i_trans.exon_bnds[1:-1] == ref_trans.exon_bnds[1:-1] )
    # next, look at how close the exons are
    return ( CDS_match, UB_match, 1-base_overlap(i_trans, ref_trans) )

def match_transcripts( clustered_transcripts ):
    trans_to_name_mapping = []
    for (chrm, strand), clusters in clustered_transcripts.iteritems():
        if VERBOSE:
            print >> sys.stderr, "Matching clusters for chr %s, strand %s"\
                    % ( chrm, strand )

        for cluster in clusters:
            matches = match_clustered_transcripts(cluster, lossy_distance)
            trans_to_name_mapping.append( matches )
    
    return trans_to_name_mapping

def build_trans_to_name_mapping( matched_transcript_clusters ):
    # for gene clusters with no association in the reference
    NOVEL_GENE_ID_CNTR = 0
    # for transcripts that have a matching gene, but no transcript
    # with the same internal structure
    NOVEL_TRANS_IDS = defaultdict( HexavijezimalCntr )
    # transcripts that have the same internal structure as a reference trans
    # we use upper case to distinguish between potential name collisions
    MODIFIED_TRANS_IDS = defaultdict(lambda: HexavijezimalCntr(upper_case=True))
    
    trans_to_name_mapping = {}
    for cluster in matched_transcript_clusters:
        for trans, match in cluster:
            if match == None:
                new_gene_id = "MGN.%i" % NOVEL_GENE_ID_CNTR
                new_trans_id = new_gene_id
                NOVEL_GENE_ID_CNTR =  NOVEL_GENE_ID_CNTR + 1
            else:
                new_gene_id = match.gene_id
                if trans.IB_key() == match.IB_key():
                    new_trans_id = "%s.%s" % (
                        match.id, MODIFIED_TRANS_IDS[new_gene_id] )
                    MODIFIED_TRANS_IDS[match.id].increment()
                else:
                    new_trans_id = "%s.%s" % (
                        match.gene_id, NOVEL_TRANS_IDS[new_gene_id] )
                    NOVEL_TRANS_IDS[new_gene_id].increment()
            
            trans_to_name_mapping[trans] = ( new_gene_id, new_trans_id )
    
    return trans_to_name_mapping

def load_sources( fp ):
    sources = {}
    for line in fp:
        data = line.strip().split( "\t" )
        sources[(data[0], data[1])] = data[2]
    return sources

def parse_arguments():
    import argparse
    desc = 'Rename transcripts.'
    parser = argparse.ArgumentParser(description=desc)    
    parser.add_argument(
        "reference_gtf", type=file,
        help='GTF file containing transcripts to take names from.' )

    parser.add_argument(
        "input_gtf", type=file,
        help='GTF file containing transcripts to rename.' )
    
    parser.add_argument(
        '--sources-fname', type=file,
        help='A sources file to be renamed.')
    
    parser.add_argument(
        '--output-prefix', '-o', type=str, 
        help='Output filename prefix. Output will be written to PREFIX.gtf and PREFIX.sources. default: stdout')
    
    parser.add_argument(
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    # make sure that the output prefix is set if a sources file is provided
    if args.sources_fname != None and args.output_prefix == None:
        raise ValueError, "If the sources file is provided you must provide an output filename prefix."
    
    # open the output files
    if args.output_prefix != None:
        ofp = open( args.output_prefix + ".gtf", "w" ) \
            if args.output_prefix != None else sys.stdout
        sources_ofp = open( args.output_prefix + ".sources", "w" ) \
            if args.sources_fname != None else None
    
    return args.reference_gtf, args.input_gtf, \
        args.sources_fname, ofp, sources_ofp


def get_non_cds_trans_name( name ):
    """This is a bit of a hack. Since we don't write out an alternate sources 
    file when we run ORF finder, we need to modify the transcript name to get 
    the correct sources. 

    """
    if -1 != name.find( "CDS" ):
        return "_".join( name.split( "_" )[:3] )
    else:
        return name

def main():
    ref_fp, input_fp, sources_fp, ofp, sources_ofp = parse_arguments()
    
    # if necessary, laod the sources file
    if sources_fp != None:
        sources = load_sources( sources_fp )
    
    # load all transcripts file and store corresponding data
    if VERBOSE: print >> sys.stderr, "Loading Transcripts..."
    ref_genes, input_genes = load_gtfs( [ref_fp.name, input_fp.name], 1 )
    
    # cluster transcripts by shared overlapping exons. We use the clusters
    # to reduce the space of mathcing transcripts to search. 
    if VERBOSE: print >> sys.stderr, "Clustering Transcripts..."
    clustered_transcripts=cluster_transcripts( 
        [ref_genes, input_genes], ["ref", "input"] )
    
    # for each cluster, find the best match
    matched_transcript_clusters = match_transcripts( clustered_transcripts )
    
    # build a mapping from the names of the input transcripts to their new names
    trans_to_name_mapping = build_trans_to_name_mapping( 
        matched_transcript_clusters )
    
    # finally, write out the new transcripts
    ofp.write( "track name='merged_renamed'\n" )
    for gene in input_genes:
        for trans in gene.transcripts:
            # get the new names
            new_gene_id, new_trans_id = trans_to_name_mapping[trans]
            
            # if we have a sources file, write out the line corresponding 
            # to this transcript
            if sources_ofp != None:
                new_sources = sources[
                    (trans.gene_id, get_non_cds_trans_name(trans.id))]
                new_data = ( new_gene_id, new_trans_id, new_sources )
                sources_ofp.write( "\t".join(new_data) + "\n" )
            
            # finally, rename the transcript and write ou tthe gtf lines
            trans.id = new_trans_id
            trans.gene_id =  new_gene_id
            ofp.write( trans.build_gtf_lines( new_gene_id, {} ) + "\n" )
    
    # cleanup the ( manually ) opened files
    if ofp != sys.stdout: ofp.close()
    if sources_ofp != None: sources_ofp.close()
    
    return

if __name__ == "__main__":
    main()


