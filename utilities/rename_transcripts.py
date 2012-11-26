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
        '--out-fname', '-o', type=argparse.FileType('w'), default=sys.stdout,
        help='Output file. default: stdout')
    
    parser.add_argument(
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.reference_gtf, args.input_gtf, \
        args.sources_fname, args.out_fname

def match_clustered_transcripts( transcripts, distance_fn ):
    ref_transcripts = [ trans for trans, source in transcripts
                        if source == 'ref' ]
    input_transcripts = [ trans for trans, source in transcripts
                          if source == 'input' ]

    # if there are no reference transcripts in this cluster, then 
    # just yield the original and a None to indicate this
    if len( ref_transcripts ) == 0:
        for t in input_transcripts:
            yield ( t, None )
        return
    
    for i_trans in input_transcripts:
        best_match = ref_transcripts[0]
        best_score = distance_fn( i_trans, ref_transcripts[0] )
        for ref_trans in ref_transcripts[1:]:
            new_score = distance_fn( i_trans, ref_trans )
            if new_score < best_score:
                best_match = ref_trans
                best_score = new_score
            yield( i_trans, ref_trans )
    
    return

class HexavijezimalCntr( object ):
    def __init__( self ):
        self.value = 2
    
    _char_mapping = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.lower()
    
    def __str__( self ):
        return_chars = []
        assert self.value > 0
        tmp_val = self.value-1
        
        while tmp_val > 0:
            remainder = (tmp_val-1)%26
            return_chars.append( self._char_mapping[remainder] )
            tmp_val =  (tmp_val + 1 - remainder)/26
        
        return "".join( reversed(return_chars) )
    
    def increment( self ):
        self.value += 1

def base_overlap( t1, t2 ):
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

def main():
    ref_fp, input_fp, sources_fp, ofp = parse_arguments()
    ofp.write( "track name='merged_renamed'\n" )
    
    # load all transcripts file and store corresponding data
    transcriptomes = load_gtfs( [ref_fp.name, input_fp.name], 1 )
    
    # cluster transcripts by shared overlapping exons. We use the clusters
    # to reduce the space of mathcing transcripts to search. 
    clustered_transcripts=cluster_transcripts(transcriptomes, ["ref", "input"])
    
    # for each cluster, find the best match
    trans_to_name_mapping = []
    for (chrm, strand), clusters in clustered_transcripts.iteritems():
        for cluster in clusters:
            matches = match_clustered_transcripts(cluster, lossy_distance)
            trans_to_name_mapping.append( matches )

    # for gene clusters with no association in the reference
    NOVEL_GENE_ID_CNTR = 0
    MODIFIED_TRANS_IDS = defaultdict( HexavijezimalCntr )
    
    for cluster in trans_to_name_mapping:
        for trans, match in cluster:
            if match == None:
                new_gene_id = "MGN.%i" % NOVEL_GENE_ID_CNTR
                new_trans_id = new_gene_id
                NOVEL_GENE_ID_CNTR =  NOVEL_GENE_ID_CNTR + 1
            else:
                new_gene_id = match.gene_id
                if 0 == exact_match_distance( trans, match ):
                    new_trans_id = match.id
                else:
                    new_trans_id = "%s.%s" % (
                        match.gene_id, MODIFIED_TRANS_IDS[new_gene_id] )
                    MODIFIED_TRANS_IDS[new_gene_id].increment()
            
            trans.gene_id = new_gene_id
            trans.id = new_trans_id
            ofp.write( trans.build_gtf_lines( new_gene_id, {} )+"\n")
    
    return

if __name__ == "__main__":
    main()
