import sys, os
from itertools import izip, combinations
from collections import defaultdict, namedtuple
import re

sys.path.append( os.path.join( os.path.dirname( __file__), 
                               "../file_types/fast_gtf_parser/" ) )
from gtf import load_gtfs, load_gtf

sys.path.append( os.path.join( os.path.dirname( __file__), "..", "file_types" ) )
from gtf_file import create_gtf_line, GenomicInterval, parse_gff_line
from wiggle import Wiggle

BlackLists = namedtuple( 'blacklists', ['regions', 'intron_chains'] )
RawData = namedtuple( 'raw_data', ['jns', 'cage', 'tss_exons'] )

TSS_EXT_MAX = 500

class HexavijezimalCntr( object ):
    def __init__( self ):
        self.value = 0
    
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

def get_introns_by_gene_name( genes ):
    raw_int_bndry_to_gene_name = defaultdict(set)
    trans_to_ref_name = defaultdict(list)
    trans_id_to_ext_coords = {}
    for ( gene_name, chrm, strand, start, stop, transcripts ) in genes:
        for trans_id, exons in transcripts:
            if len( exons ) <= 2: continue
            
            trans_id_to_ext_coords[ trans_id ] = (exons[0], exons[-1], gene_name)
            trans_to_ref_name[ ( chrm, strand, tuple(exons[1:-1]) ) ].append( trans_id )
            for bndry in exons[1:-1]:
                raw_int_bndry_to_gene_name[ ( chrm, strand, bndry ) ].add( gene_name )
    
    int_bndry_to_gene_name = {}
    for key, gene_ids in raw_int_bndry_to_gene_name.iteritems():
        int_bndry_to_gene_name[key] = tuple( sorted( gene_ids ) )
    
    return int_bndry_to_gene_name, trans_to_ref_name, trans_id_to_ext_coords

def find_associated_genes( exons, int_bndry_to_gene_name, chrm, strand ):
    genes = set()
    for bndry in exons[1:-1]:
        if (chrm, strand, bndry) in int_bndry_to_gene_name:
            genes.add( tuple(int_bndry_to_gene_name[ (chrm, strand, bndry) ] ) )
    
    return tuple(sorted(genes))

def is_blacklisted( region, bl_regions ):
    start, stop = region
    for bl_start, bl_stop in bl_regions:
        if not( stop < bl_start or start > bl_stop ):
            return True
    
    return False

def get_matching_ends( exons, trans_ids, trans_id_to_ext_coords ):
    closest_match_diff = None
    for trans_id in trans_ids:
        start, stop, gene_id = trans_id_to_ext_coords[trans_id]
        if closest_match_diff == None:
            matching_gene_trans_names = [(gene_id, trans_id),]
            closest_match_diff = abs(start - exons[0]) + abs(stop - exons[-1])
        elif abs(start - exons[0]) + abs(stop - exons[-1]) < closest_match_diff:
            matching_gene_trans_names = [(gene_id, trans_id),]
            closest_match_diff = abs(start - exons[0]) + abs(stop - exons[-1])
        elif abs(start - exons[0]) + abs(stop - exons[-1]) == closest_match_diff:
            matching_gene_trans_names.append( (gene_id, trans_id) )
    
    return matching_gene_trans_names

def name_transcripts( ann_genes, ref_genes, good_gene_merges, blacklists ):
    int_bndry_to_gene_name, trans_to_ref_name, trans_id_to_ext_coords \
        = get_introns_by_gene_name(ref_genes)
    
    ann_ofp = open( "filtered.renamed.transcripts.gtf", "w" )
    
    novel_gene_names_set = set()
    novel_gene_cntr = 0
    novel_trans_cntr = defaultdict( HexavijezimalCntr )
    # even if the transcript match is identical, we may have alternate 5' and 3'
    # ends which means multiple 'identical' transcripts. So keep a counter. If 
    # the counter is at 0 *for that ref trans id*, keep the ref name. Otherwise,
    # increment it by one.
    ident_cntr = defaultdict( HexavijezimalCntr )
    # for novel splices, againb, index by gene name
    alternative_splice_cntr = defaultdict( HexavijezimalCntr )
    
    # keep track of the genes we have seen
    observed_reference_trans = set()
    
    trans_id_to_names = defaultdict( list )
    
    def set_novel_trans_name( gene_name, trans_id, novel_gene_cntr ):
        if gene_name not in novel_gene_names_set:
            novel_gene_cntr += 1
            novel_gene_names_set.add( gene_name )
        
        new_gene_name = "mgn%s" % str(novel_gene_cntr).zfill( 5 )
        if new_gene_name not in novel_trans_cntr:
            novel_trans_cntr[ new_gene_name ].increment()
        novel_trans_cntr[ new_gene_name ].increment()
        
        new_trans_name = new_gene_name + '.' + str( novel_trans_cntr[ new_gene_name ] )
        trans_id_to_names[trans_id].append( ( new_gene_name, new_trans_name ) )
        
        return novel_gene_cntr
    
    def set_ident_trans_name_s( trans_id, gene_and_trans_names ):
        for ref_gene_name, ref_trans_name in gene_and_trans_names:
            observed_reference_trans.add( ref_trans_name )
            
            if ident_cntr[ ref_trans_name ].value > 0:
                new_trans_name = ref_trans_name + "." + str(ident_cntr[ ref_trans_name ])
            else:
                ident_cntr[ ref_trans_name ].increment()
                new_trans_name = ref_trans_name
            
            ident_cntr[ ref_trans_name ].increment()
            
            trans_id_to_names[ trans_id ].append( ( 
                    ref_gene_name, new_trans_name ) )
            
        return
    
    def find_and_set_ident_trans_name_s( chrm, strand, exons, assoc_genes, trans_id ):
        # check to see if this is an exact match
        if (chrm, strand, tuple(exons[1:-1])) in trans_to_ref_name:
            # get all of the reference transcripts with this internal structure
            ref_trans_names = trans_to_ref_name[(chrm, strand, tuple(exons[1:-1]))]
            
            # if there is only one transcript with this internal structure
            if len(ref_trans_names) == 1:
                set_ident_trans_name_s( trans_id, ((assoc_genes[0][0], ref_trans_names[0]),) )
            # if this is a disisctronic or transcript with alternative ends
            else:
                # find the transcript(s) with the closest distal end to the annotation exons
                ref_gene_trans_names = get_matching_ends( 
                    exons, ref_trans_names, trans_id_to_ext_coords )
                set_ident_trans_name_s( trans_id, ref_gene_trans_names )
            
        # if it's not, it must be a partial match
        else:
            # just pick the first accosiated ref gene in the case of dicistronics
            assoc_ref_gene = assoc_genes[0][0]
            
            # update the counters and set new transcript names
            if assoc_ref_gene not in alternative_splice_cntr:
                alternative_splice_cntr[ assoc_ref_gene ].increment()
            alternative_splice_cntr[ assoc_ref_gene ].increment()
            new_trans_name = "{0}.{1}".format( 
                assoc_genes[0][0], alternative_splice_cntr[ assoc_ref_gene ])
            trans_id_to_names[ trans_id ].append( ( 
                    assoc_ref_gene, new_trans_name ) )
        
        return
    
    def find_and_set_merge_trans_name_s( assoc_genes, chrm, strand, exons, trans_id ):
        # if there are any shared bndries by ref genes then
        # this is likely a merge between a 5.45 gene dicistronic and another gene
        # or this is a complex dicistronic locus: e.g.:ABCB7
        if any( len(i_assoc_genes) > 1 for i_assoc_genes in assoc_genes ):
            # Deal with complex dicistronic loci by:
            #
            # 1) find all genes that share all of the boundaries in this transcript
            # 2) if this is only one gene then see if it is an exact match to a transcript
            actual_matchs = set( assoc_genes[0] )
            for i_assoc_genes in assoc_genes:
                actual_matchs = actual_matchs.intersection( i_assoc_genes )
            if len( actual_matchs ) != 1: return
            if (chrm, strand, tuple(exons[1:-1])) not in trans_to_ref_name: return
            
            ref_trans_names = trans_to_ref_name[(chrm, strand, tuple(exons[1:-1]))]
            # if there is only one transcript with this internal structure
            if len(ref_trans_names) == 1:
                set_ident_trans_name_s( trans_id, ((assoc_genes[0][0], ref_trans_names[0]),) )
            # if this is a disisctronic or transcript with alternative ends
            else:
                # find the transcript(s) with the closest distal end to the annotation exons
                ref_gene_trans_names = get_matching_ends( 
                    exons, ref_trans_names, trans_id_to_ext_coords )
                set_ident_trans_name_s( trans_id, ref_gene_trans_names )
            
            return
        
        assoc_genes = tuple(sorted( i_assoc_genes[0] for i_assoc_genes in assoc_genes))
        
        # make sure that it's white listed
        if assoc_genes not in good_gene_merges: return
        
        new_gene_name = "/".join( assoc_genes )
        alternative_splice_cntr[ new_gene_name ].increment()
        new_trans_name = "%s.%s" % ( 
            new_gene_name, alternative_splice_cntr[ new_gene_name ] )
        trans_id_to_names[ trans_id ].append( ( new_gene_name, new_trans_name ) )
        
        return
    
    
    for ( gene_name, chrm, strand, start, stop, transcripts ) in ann_genes:
        # check to see if this overlaps a black list region
        if is_blacklisted( (start, stop), blacklists.regions[(chrm, strand)] ): 
            continue
        
        for trans_id, exons in transcripts:
            if len( exons ) == 2: continue
            # if this is a blacklisted intron chain
            #if ( chrm, strand, tuple(exons[1:-1]) ) in blacklists.intron_chains:
            #    continue
            
            assoc_genes = find_associated_genes(
                exons, int_bndry_to_gene_name, chrm, strand )
            
            # if this is a completely novel transcript
            if len( assoc_genes ) == 0:
                novel_gene_cntr = set_novel_trans_name( 
                    gene_name, trans_id, novel_gene_cntr )
            # if there is only 1 corresponding ref gene (or a dicistronic)
            elif len( assoc_genes ) == 1:
                find_and_set_ident_trans_name_s( chrm, strand, exons, assoc_genes, trans_id )
            # if this is a gene merge
            elif len( assoc_genes ) > 1:
                find_and_set_merge_trans_name_s( assoc_genes, chrm, strand, exons, trans_id )
            
            # write all new gene and trans names for this transcript if any have been assigned
            if trans_id not in trans_id_to_names: continue
            for new_gene_name, new_trans_name in trans_id_to_names[ trans_id ]:
                for start, stop in zip( exons[0::2], exons[1::2] ):
                    region = GenomicInterval( chrm, strand, start, stop )
                    ann_ofp.write( create_gtf_line( 
                            region, new_gene_name, new_trans_name, {}, 
                            feature='exon', source='grit' ) + "\n" )
    
    ann_ofp.close()
    
    return observed_reference_trans

def process_missed_transcripts( observed_reference_trans,
                                raw_data, ref_genes, ann_genes ):
    ann_tss_exons = defaultdict( set )
    ann_int_bndries_to_tss_exon = defaultdict( set )
    for ( gene_name, chrm, strand, start, stop, transcripts ) in ann_genes:
        for trans_id, exons in transcripts:
            int_tss_coord = exons[1] if strand == '+' else exons[-2]
            ext_tss_coord = exons[0] if strand == '+' else exons[-1]
            ann_tss_exons[(chrm, strand, int_tss_coord)].add(ext_tss_coord)
            
            int_exons = exons[2:-1] if strand == '+' else exons[1:-2]
            tss_exon = exons[:2] if strand == '+' else exons[-2:]
            ann_int_bndries_to_tss_exon[(chrm, strand, tuple(int_exons))].add(
                tuple(tss_exon))
    
    def shift_fp_bndry( exons, chrm, strand, int_tss_coord ):
        ext_tss_coords = raw_data.tss_exons[(chrm, strand, int_tss_coord)]
        ref_ext_tss_coord = exons[0] if strand == '+' else exons[-1]
        
        if ((min( ext_tss_coords ) > ref_ext_tss_coord and strand == "+") or \
                    (max( ext_tss_coords ) < ref_ext_tss_coord and strand == "-")) :
            return None
        
        min_tss_extention = None
        for ext_tss_coord in ext_tss_coords:
            # skip tss exons that are shorter
            if (ext_tss_coord >  ref_ext_tss_coord and strand == "+") or \
                    (ext_tss_coord < ref_ext_tss_coord and strand == "-"):
                continue
            
            # find the closest tss external coord that is longer than the ref
            if min_tss_extention == None or \
                    abs( ext_tss_coord - ref_ext_tss_coord ) < min_tss_extention:
                min_tss_extention = abs( ext_tss_coord -  ref_ext_tss_coord )
        
        if min_tss_extention < TSS_EXT_MAX:
            # convert exons to use updated 5' boundary
            exons = [(ref_ext_tss_coord - min_tss_extention),] + exons[1:] if strand == '+' else \
                exons[:-1] + [(ref_ext_tss_coord + min_tss_extention),]
            
            return exons
        
        return None
    
    
    unobserved_ofp = open( 'unobserved_refernce_transcripts.gtf', 'w' )
    shifted_unobserved_ofp = open( 'unobserved_shifted_refernce_transcripts.gtf', 'w' )
    
    def write_trans( exons, chrm, strand, gene_name, trans_id, fp ):
        for start, stop in zip( exons[0::2], exons[1::2] ):
            region = GenomicInterval( chrm, strand, start, stop )
            fp.write( create_gtf_line( 
                    region, gene_name, trans_id, {}, 
                    feature='exon', source='grit' ) + "\n" )
        return
    
    
    for ( gene_name, chrm, strand, start, stop, transcripts ) in ref_genes:
        for trans_id, exons in transcripts:
            # if we found this transcript then we don't need to handle it here
            if trans_id in observed_reference_trans:
                continue
            
            # all single exon genes should be printed out as not found
            if len( exons ) == 2:
                write_trans( exons, chrm, strand, gene_name, trans_id, unobserved_ofp )
                continue
            
            int_tss_coord = exons[1] if strand == '+' else exons[-2]
            if (chrm, strand, int_tss_coord) in ann_tss_exons:
                ref_ext_tss_coord = exons[0] if strand == '+' else exons[-1]
                closest_tss_match = None
                for ext_tss_coord in ann_tss_exons[(chrm, strand, int_tss_coord)]:
                    if closest_tss_match == None:
                        closest_tss_match = ext_tss_coord
                    elif abs(ext_tss_coord - ref_ext_tss_coord) < \
                            abs(closest_tss_match - ref_ext_tss_coord):
                        closest_tss_match = ext_tss_coord
                
                # convert exons to use updated 5' boundary
                exons = [ext_tss_coord,] + exons[1:] if strand == '+' else \
                    exons[:-1] + [ext_tss_coord,]
                write_trans( exons, chrm, strand, 
                             gene_name, trans_id, shifted_unobserved_ofp )
                continue
            
            elif (chrm, strand, int_tss_coord) in raw_data.tss_exons:
                tmp_exons = shift_fp_bndry( exons, chrm, strand, int_tss_coord )
                if tmp_exons != None:
                    write_trans( tmp_exons, chrm, strand, 
                                 gene_name, trans_id, shifted_unobserved_ofp )
                    continue
            
            # try to use cage to redefine remove flybase transcripts
            int_exons = exons[2:-1] if strand == '+' else exons[1:-2]
            if (chrm, strand, tuple(int_exons)) in ann_int_bndries_to_tss_exon:
                ann_tss_exons = ann_int_bndries_to_tss_exon[
                    (chrm, strand, tuple(int_exons))]
                ref_tss_exon = exons[:2] if strand == '+' else exons[-2:]
                
                #tss_cvrg = '{0:.1f}'.format( \
                #    rd.tss_wig[chrm_strand][cage_tss_exon[0]:cage_tss_exon[1]].sum()) \
                #    if chrm_strand in rd.tss_wig else 0
            
            # write out flybase models
            write_trans( exons, chrm, strand, gene_name, trans_id, unobserved_ofp )
            
    
    unobserved_ofp.close()
    shifted_unobserved_ofp.close()
    
    return

def get_good_gene_merges( fname ):
    good_gene_merges = set()
    with open( fname ) as fp:
        for line in fp:
            gene_names = tuple(sorted(line.split()))
            good_gene_merges.add( gene_names )
            # add all of the comboinations
            for i in xrange( 2, len(gene_names) ):
                for new_gene_names in combinations(gene_names, i):
                    good_gene_merges.add( tuple(sorted(new_gene_names)) )
    
    return good_gene_merges

def get_blacklist_regions( fname, ref_genes ):
    gene_regions = {}
    for gene_name, chrm, strand, start, stop, transcripts in ref_genes:
        gene_regions[gene_name] = ( chrm, strand, start, stop )
    
    blacklist_regions = defaultdict( set )
    with open( fname ) as fp:
        for line in fp:
            gene_name = line.strip()
            try:
                chrm, strand, start, stop = gene_regions[ gene_name ]
            except KeyError:
                print >> sys.stderr, "WARNING: Can't find black list gene '%s'" % gene_name
                continue
            
            blacklist_regions[(chrm, strand)].add( (start, stop ) )
    
    for key in blacklist_regions:
        blacklist_regions[key] = sorted( blacklist_regions[key] )
    
    return blacklist_regions

def get_blacklist_intron_chains( bl_trans_fn ):
    bl_trans = load_gtf( bl_trans_fn )
    
    bl_intron_chains = set()
    for ( gene_name, chrm, strand, start, stop, transcripts ) in bl_trans:
        for trans_id, exons in transcripts:
            bl_intron_chains.add( (chrm, strand, tuple(exons[1:-1]) ) )
    
    return bl_intron_chains

def parse_jns_file( jns_fp ):
    all_jns = set()
    for line in jns_fp:
        gffl = parse_gff_line( line )
        if gffl == None: continue
        
        all_jns.add( 
            (gffl.region.chr, gffl.region.strand, 
             gffl.region.start, gffl.region.stop) )
    
    jns_fp.close()
    return all_jns

def parse_tss_exons_file( tss_exons_fp ):
    tss_exons = defaultdict( set )
    for line in tss_exons_fp:
        gffl = parse_gff_line( line )
        if gffl == None: continue
        (int_coord, ext_coord) = (gffl.region.stop, gffl.region.start) \
            if gffl.region.strand == '+' else (gffl.region.start, gffl.region.stop)
        tss_exons[(gffl.region.chr, gffl.region.strand, int_coord)].add( ext_coord )
    
    return tss_exons

def build_objects( ann_fn, ref_fn, good_gene_merges_fn,
                   bl_regions_fn, cage_fps, chrm_sizes_fp, tss_exons_fp ):#, bl_trans_fn, jns_fp ):
    ann_genes, ref_genes =  load_gtfs( ( ann_fn, ref_fn ), 2 )
    
    good_gene_merges = get_good_gene_merges( good_gene_merges_fn )
    
    # get and store blacklist regions
    bl_regions = get_blacklist_regions( bl_regions_fn, ref_genes )
    #bl_intron_chains = get_blacklist_intron_chains( bl_trans_fn )
    blacklists = BlackLists( bl_regions, None )#bl_intron_chains )
    
    # parse raw data files
    #all_jns = parse_jns_file( jns_fp )
    #cage_wig = Wiggle( chrm_sizes_fp, cage_fps )
    cage_wig = None
    tss_exons = parse_tss_exons_file( tss_exons_fp )
    tss_exons_fp.close()
    
    # pack all raw data
    raw_data = RawData( None, cage_wig, tss_exons )
    
    return ann_genes, ref_genes, good_gene_merges, blacklists, raw_data

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(
        description = '' )
    parser.add_argument(
        'annotation',
        help='.' )
    parser.add_argument(
        'reference',
        help=".")
    parser.add_argument(
        'good_gene_merges',
        help=".")
    parser.add_argument(
        'blacklist_gene_names',
        help=".")
    parser.add_argument(
        'tss_exons_gff', type=file,
        help=".")
    parser.add_argument(
        '--chrm_sizes', type=file,
        help=".")
    parser.add_argument(
        '--cage_wigs', type=file, nargs='*',
        help=".")
    #parser.add_argument(
    #    'blacklist_trans_gtf',
    #    help=".")
    #parser.add_argument(
    #    'jns_gtf', type=file, 
    #    help=".")
    
    parser.add_argument( 
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.annotation, args.reference, args.good_gene_merges, \
        args.blacklist_gene_names, args.tss_exons_gff, args.cage_wigs, args.chrm_sizes#, args.blacklist_trans_gtf, args.jns_gtf

def main():
    ann_fn, ref_fn, good_gene_merges_fn, bl_regions_fn, tss_exons_fp, cage_fps, \
        chrm_sizes_fp = parse_arguments()#, bl_trans_fn, jns_fp \
    
    ann_genes, ref_genes, good_gene_merges, blacklists, raw_data \
        = build_objects( ann_fn, ref_fn, good_gene_merges_fn, \
                             bl_regions_fn, cage_fps, chrm_sizes_fp, \
                             tss_exons_fp)#, bl_trans_fn, jns_fp )
    
    observed_reference_trans = name_transcripts( ann_genes, ref_genes, good_gene_merges, blacklists )
    
    process_missed_transcripts( observed_reference_trans, raw_data, ref_genes, ann_genes )
    
    return

if __name__ == "__main__":
    main()
