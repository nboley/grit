import sys, os
from itertools import izip

VERBOSE = False

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../file_types/fast_gtf_parser/" ) )
from gtf import load_gtf
sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../file_types/" ) )
from genomic_intervals import GenomicInterval
from gtf_file import iter_gff_lines

def get_elements_from_gene( gene, get_tss=True, get_jns=True, \
                                get_tes=True, get_exons=False ):
    tss_exons = set()
    tes_exons = set()
    introns = set()
    exons = set()
    
    chrm, strand = gene[1], gene[2]
    transcripts = gene[-1]
    
    for trans in transcripts:
        bndries = trans[1]

        fp_region = GenomicInterval(chrm, strand, bndries[0], bndries[1])
        tp_region = GenomicInterval(chrm, strand, bndries[-2], bndries[-1])
        if strand == '+':
            if get_tss:
                tss_exons.add( fp_region )
            if get_tes:
                tes_exons.add( tp_region )
        else:
            if strand != '-':
                print >> sys.stderr, "BADBADBAD", strand
                continue
            assert strand == '-'
            if get_tss:
                tss_exons.add( tp_region )
            if get_tes:
                tes_exons.add( fp_region )
        
        if get_jns:
            for start, stop in izip( bndries[1:-2:2], bndries[2:-1:2] ):
                # add and subtract 1 to ge tthe inclusive intron boundaries,
                # rather than the exon boundaries
                if start >= stop:
                    continue
                introns.add( GenomicInterval(chrm, strand, start+1, stop-1) )

        if get_exons:
            for start, stop in izip( bndries[::2], bndries[1::2] ):
                exons.add( GenomicInterval(chrm, strand, start, stop) )
    
    return tss_exons, introns, tes_exons, exons

def get_element_sets( genes, get_tss=True, get_jns=True, \
                          get_tes=True, get_exons=True ):
    tss_exons = set()
    introns = set()
    tes_exons = set()
    exons = set()
    for gene in genes:
        i_tss_exons, i_introns, i_tes_exons, i_exons = \
            get_elements_from_gene( gene, get_tss, get_jns, get_tes, get_exons )
        
        tss_exons.update( i_tss_exons )
        introns.update( i_introns  )
        tes_exons.update( i_tes_exons )
        exons.update( i_exons )

    return sorted(tss_exons), sorted(introns), sorted( exons ), sorted(tes_exons)

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(\
        description='Extract tss, tes, exons, junctions from a gtf file.')
    parser.add_argument( 'gtf', help='A gtf file containing transcripts.')
    
    parser.add_argument( '--tss-fname', \
        help='Output filename for tss exons. Default: ignore tss exons' )
    parser.add_argument( '--tes-fname', \
        help='Output filename for tes exons. Default: ignore tes exons' )
    parser.add_argument( '--jn-fname', \
        help='Output filename for introns ( junctions ). Default: ignore junctions' )
    parser.add_argument( '--exon-fname', \
        help='Output filename for exons. Default: ignore exons' )

    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
        help='Whether or not to print status information.')
    
    args = parser.parse_args()
    
    global VERBOSE
    VERBOSE = args.verbose
    
    if all(i == None for i in ( \
            args.tss_fname, args.tes_fname, args.jn_fname, args.exon_fname)):
        e_temp = "%s: error: You must choose at least one output element type."
        parser.print_help()
        print >> sys.stderr, e_temp % __file__
        sys.exit( 1  )
    
    return args.gtf, args.tss_fname, args.tes_fname, \
        args.jn_fname, args.exon_fname
                      
def main():
    gtf_fname, tss_fname, tes_fname, jn_fname, exon_fname = parse_arguments()

    # load the genes and build sorted, unique lists
    genes = load_gtf( gtf_fname )
    tss_exons, introns, exons, tes_exons = get_element_sets( \
        genes, tss_fname != None, jn_fname != None, tes_fname != None )

    def write_to_disk( elements_iter, op_fname, feature_name ):
        """write the elements to disk """
        # if we didn't provide a filename, dont do anything
        if None == op_fname:
            return
        
        with open( op_fname, "w" ) as op_fp:
            lines_iter = iter_gff_lines( elements_iter, feature=feature_name )
            op_fp.write( "\n".join( lines_iter ) )
        
        return
    
    write_to_disk( tss_exons, tss_fname, "tss" )
    write_to_disk( introns, jn_fname, "intron" )
    write_to_disk( exons, exon_fname, "exons" )
    write_to_disk( tes_exons, tes_fname, "tes" )
    
    return

if __name__ == "__main__":
    main()

