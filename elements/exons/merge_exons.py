import sys, os

VERBOSE = False

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "..", "..", "file_types" ) )
from gtf_file import parse_gff_line, iter_gff_lines

def load_exons( exon_fp ):
    exon_regions = set()
    for line in exon_fp:
        data = parse_gff_line( line )
        if data == None: continue
        exon_regions.add( data.region )
    return exon_regions

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description='Merge exon files.')
    parser.add_argument( 'exons', nargs="+", type=file, \
                             help='A gff file containing exons.')
    
    parser.add_argument( '--out-fname', '-o', \
        help='Output file will be written to default. default: stdout')
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
        help='Whether or not to print status information.')
    parser.add_argument( '--feature-type', default="exon", \
        help='The feature type to put in the output gff line.')

    args = parser.parse_args()
    
    out_fp = open( args.out_fname, "w" ) if args.out_fname else sys.stdout
    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.exons, out_fp, args.feature_type


def main():
    exons_fps, out_fp, feature_type = parse_arguments()
    all_exons = set()
    for exons_fp in exons_fps:
        if VERBOSE: print >> sys.stderr, "Loaded", exons_fp.name
        all_exons.update( load_exons( exons_fp ) )
    
    regions_iter = sorted( all_exons )
    line_iter = iter_gff_lines( regions_iter, feature=feature_type )
    out_fp.write( "\n".join( line_iter  ) + "\n" )
    
    return

if __name__ == "__main__":
    main()

