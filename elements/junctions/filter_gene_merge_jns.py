import os, sys
import numpy
from collections import defaultdict

sys.path.append(os.path.join(os.path.dirname(__file__), "../../", "file_types"))
from junctions_file import parse_jn_gff

sys.path.append(os.path.join(os.path.dirname(__file__), 
                             "../../", "file_types", "fast_gtf_parser"))
from gtf import load_gtf

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
    
    parser.add_argument( 
        '--write-non-gene-merge-jns', default=False, action='store_true',
        help='Write out gene merge and non gene merge jns that pass the filter.')

    args = parser.parse_args()
        
    return getattr(args, 'jns-gff'), getattr(args, 'ref-gtf'), \
        args.write_good_gene_merge_jns, args.write_bad_gene_merge_jns, \
        args.write_non_gene_merge_jns


def main():
    jns_fp, ref_fp, write_good_gm_jns, write_bad_gm_jns, write_non_gm_jns \
        = parse_arguments()
    
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

    for jn in jns:
        is_gene_merge, genes = is_gene_merge_jn(gene_bndrys, bndry_to_gene_mapping, jn)
        if is_gene_merge:
            gene_merge_line = ' merged_genes "%s";' % ','.join(genes)
            if float(jn.uniq_cnt)/jn.cnt > 0.5:
                if write_good_gm_jns:
                    print jn.build_gff_line() + gene_merge_line
            else:
                if write_bad_gm_jns:
                    print jn.build_gff_line() + gene_merge_line
        elif write_non_gm_jns:
            print jn.build_gff_line()
    
if __name__ == "__main__":
    main()
