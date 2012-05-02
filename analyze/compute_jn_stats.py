import os
import sys

sys.path.append( os.path.join( os.path.dirname(__file__), "../.." ) )
from gene_models import parse_gff_line, GenomicInterval

def get_ann_jns( ann_fname ):
    ann_jns = set()
    with open( ann_fname ) as ann_fp:
        for line in ann_fp:
            ann_jns.add( parse_gff_line( line )[-1] )
    return ann_jns

def cnt_jns( fname, ann_jns ):
    total_jns = 0
    jns_in_ann = 0
    with open( fname ) as fp:
        for line in fp:
            total_jns += 1
            interval = parse_gff_line( line )[-1]
            if interval in ann_jns:
                jns_in_ann += 1
    
    return jns_in_ann, total_jns

ann_fname = sys.argv[1]
ann_jns = get_ann_jns( ann_fname )
overlap, total = cnt_jns( ann_fname, ann_jns )
print overlap, total, ann_fname
for fname in sys.argv[2:]:
    overlap, total = cnt_jns( fname, ann_jns )
    print overlap, total, fname

