# Copyright (c) 2011-2012 Nathan Boley

import sys, os
sys.path.append( os.path.join( os.path.dirname(__file__), "../.." ) )
from gene_models import parse_gff_line, GenomicInterval

def get_tsss( fname ):
    tsss = set()
    with open( fname ) as fp:
        for line in fp:
            interval = parse_gff_line( line )[-1]
            if interval.strand == '-':
                internal_coord = interval.start
            else:
                internal_coord = interval.stop
            tsss.add( ( interval.chr, interval.strand, internal_coord ) )
    
    return tsss

ann_tsss = get_tsss( sys.argv[1] )
other_tsss = get_tsss( sys.argv[2] )
print len( ann_tsss ), len( other_tsss ), len( ann_tsss.intersection( other_tsss ) )
