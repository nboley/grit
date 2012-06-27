# Copyright (c) 2011-2012 Nathan Boley

import sys
import re

cluster_pat = re.compile( 'gene_id "(cluster_\d+)";' )

clusters_fp = open( sys.argv[1] )
transcripts_fp = open( sys.argv[2] )

# build the set of clusters
clusters = set( line.strip() for line in clusters_fp )

for line in transcripts_fp:
    # get the cluster id
    cluster_id = cluster_pat.findall( line )[0]
    if cluster_id in clusters:
        print line.strip()


