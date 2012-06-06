import os, sys
import numpy
from collections import defaultdict

sys.path.append(os.path.join(os.path.dirname(__file__), "../../", "file_types"))
from junctions_file import parse_jn_gff

sys.path.append(os.path.join(os.path.dirname(__file__), 
                             "../../", "file_types", "fast_gtf_parser"))
from gtf import load_gtf

MIN_CNT_RATIO= 0.1
USE_MIN_ACCEPTOR = False

# store string to give the status of the gene merge 
FAILS_UNIQ = ' gene_merge_status "FAILS_UNIQ";'
FAILS_SCORE = ' gene_merge_status "FAILS_CNT_RATIO %.2e";'
PASSES = ' gene_merge_status "GOOD";'

def is_gene_merge_jn( gene_bndrys, gene_names, jn ):
    if (jn.chrm, jn.strand) not in gene_bndrys:
        return False, []

    # find the bndrys that it intersects
    starts, stops = gene_bndrys[ (jn.chrm, jn.strand) ]
    i = starts.searchsorted( jn.start )-1

    # if this isn't even in a gene, return False
    if jn.start < starts[i] or jn.start > stops[i]:
        return False, []
    start_gene = gene_names[(jn.chrm, jn.strand, starts[i], stops[i])]

    i = starts.searchsorted( jn.stop ) - 1

    # if this isn't even in a gene, return False
    if jn.stop < starts[i] or jn.stop > stops[i]:
        return False, []
    stop_gene = gene_names[(jn.chrm, jn.strand, starts[i], stops[i])]
    
    return start_gene != stop_gene, [start_gene, stop_gene]

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Find and analyze gene merge junctions.')
    
    parser.add_argument('jns-gff', type=file, help="GFF file to load jns from.")

    parser.add_argument('ref-gtf', type=file, help="GTF file to load genes from.")
    
    parser.add_argument( '--write-good-gene-merge-jns', 
                         default=False, action='store_true',
                         help='Write out gene merge jns that pass the filter.')
    parser.add_argument( '--write-bad-gene-merge-jns', 
                         default=False, action='store_true',
                         help='Write out gene merge jns that dont pass the filter.')

    parser.add_argument( '--gene-merges-summary-fname', type=str,
                         help='Summarize merges at teh gene level, and write them to this file..')
    
    parser.add_argument( 
        '--write-non-gene-merge-jns', default=False, action='store_true',
        help='Write out gene merge and non gene merge jns that pass the filter.')

    args = parser.parse_args()
        
    return getattr(args, 'jns-gff'), getattr(args, 'ref-gtf'), \
        args.write_good_gene_merge_jns, args.write_bad_gene_merge_jns, \
        args.write_non_gene_merge_jns, args.gene_merges_summary_fname


def main():
    jns_fp, ref_fp, write_good_gm_jns, write_bad_gm_jns, write_non_gm_jns, \
        gene_merges_summary_fname = parse_arguments()
    
    jns = parse_jn_gff( jns_fp )
    ref_genes = load_gtf( ref_fp.name )

    gene_bndrys = defaultdict( list )
    bndry_to_gene_mapping = {}

    for name, chrm, strand, start, stop, trans in ref_genes:
        gene_bndrys[ (chrm, strand) ].append( (start, stop) )
        bndry_to_gene_mapping[ (chrm, strand, start, stop) ] = name
    
    for key, bndrys in gene_bndrys.iteritems():
        starts, stops = zip( *sorted(gene_bndrys[key]) )
        gene_bndrys[key] = ( numpy.array(starts), numpy.array(stops) )

    # build a mapping from donor site to acceptor
    donor_acceptor_mapping = {}
    for jn in jns:
        donor = jn.start if jn.strand == '+' else jn.stop
        acceptor = jn.stop if jn.strand == '+' else jn.start
        key = ( jn.chrm, jn.strand, donor )
        if key in donor_acceptor_mapping:
            if USE_MIN_ACCEPTOR:
                other_acceptor = donor_acceptor_mapping[key][0]
                assert other_acceptor != acceptor
                if abs(other_acceptor-donor) > abs(acceptor-donor):
                    donor_acceptor_mapping[key] = ( acceptor, jn.cnt )
            else:
                donor_acceptor_mapping[key] = \
                    max( jn.cnt, donor_acceptor_mapping[key] )
        else:
            if USE_MIN_ACCEPTOR:
                donor_acceptor_mapping[key] = ( acceptor, jn.cnt )
            else:
                donor_acceptor_mapping[key] = jn.cnt
    
    gene_merges = defaultdict( lambda: [0,0,0] )
    
    for jn in jns:
        is_gene_merge, genes = is_gene_merge_jn(gene_bndrys, bndry_to_gene_mapping, jn)
        if is_gene_merge:
            gene_merge_line = ' merged_genes "%s";' % ','.join(genes)
            if USE_MIN_ACCEPTOR:
                min_acceptor_cnt = donor_acceptor_mapping[
                    (jn.chrm, jn.strand, jn.start if jn.strand == '+' else jn.stop) ][1]
            else:
                min_acceptor_cnt = donor_acceptor_mapping[
                    (jn.chrm, jn.strand, jn.start if jn.strand == '+' else jn.stop) ]
                
            if float(jn.uniq_cnt)/jn.cnt < 0.5:
                gene_merges[ tuple(genes) ][0] += 1
                if write_bad_gm_jns:
                    print jn.build_gff_line() + gene_merge_line + FAILS_UNIQ
            # if it passes the uniqueness, check to see if it passes the ratio threshold
            elif float(jn.cnt)/min_acceptor_cnt < MIN_CNT_RATIO:
                gene_merges[ tuple(genes) ][1] += 1
                if write_bad_gm_jns:
                    print jn.build_gff_line() + gene_merge_line \
                        + FAILS_SCORE % (float(jn.cnt)/min_acceptor_cnt)
            # if it passes both tests
            else:
                gene_merges[ tuple(genes) ][2] += 1
                if write_good_gm_jns:
                    print jn.build_gff_line() + gene_merge_line + PASSES
        
        elif write_non_gm_jns:
            print jn.build_gff_line()
    
    if gene_merges_summary_fname != None:
        gene_merge_of = open( gene_merges_summary_fname, "w" )
        print >> gene_merge_of, "\t".join(
            ("Gene1".ljust(30), "Gene2".ljust(30), 
             "F_UNIQ", "F_SIG", "PASS", "STATUS"))
                                            
        for key in sorted( gene_merges ):
            args = list( key )
            args[0] = args[0].ljust(30)
            args[1] = args[1].ljust(30)
            args.extend( gene_merges[key] )
            if args[-2] == 0 and args[-3] == 0:
                args.append( "GOOD" )
            else:
                args.append( "BAD" )

            print >> gene_merge_of, "\t".join( map( str, args ) )
        
        gene_merge_of.close()
    
if __name__ == "__main__":
    main()
