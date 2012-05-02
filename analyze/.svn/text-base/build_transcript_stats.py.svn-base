import sys
import re

from collections import defaultdict

fp = open( sys.argv[1] )

transcripts = defaultdict( set ) 

for line in fp:
    data = line.strip().split("\t")
    meta_data = data[-1].split(';')
    if len(meta_data) >= 2:
        gene_id = meta_data[0].split(" ")[1][1:-1]
        transcript_id = meta_data[1].split()[-1][1:-1]
        transcripts[ gene_id ].add( transcript_id )

with open( "transcript_cnts.txt", "w" ) as cnt_fp:
    for entry in transcripts.values():
        print >> cnt_fp, len( entry )


with open( "small_transcript_locis.txt", "w" ) as small_loci_fp:
    for key, entry in transcripts.iteritems():
        if len( entry ) > 100: continue
        print >> small_loci_fp, key


