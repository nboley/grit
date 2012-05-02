import sys, os

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../../file_types/" ) )
from gtf_file import parse_gff_line, iter_gff_lines

def parse_exon_file( exons_fname ):
    exons = set()
    for line in open(exons_fname):
        gff = parse_gff_line( line )
        if gff == None: continue
        exons.add( gff.region )
    return exons

def main( fn_1, fn_2 ):
    r1 = parse_exon_file( fn_1 )
    r2 = parse_exon_file( fn_2 )
    print "\n".join( (iter_gff_lines( sorted(r2 - r1) ) ) )
    return len(r1.intersection( r2 )), len( r1 ), len( r2 )

if __name__ == "__main__":
    print >> sys.stderr, main( sys.argv[1], sys.argv[2] )
