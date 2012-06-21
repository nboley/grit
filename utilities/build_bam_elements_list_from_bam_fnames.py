# Copyright (c) 2011-2012 Nathan Boley

import sys, os
from collections import defaultdict

samples = defaultdict(list)

for fname in sys.argv[1:]:
    # skip non bam files
    if not fname[-4:] == ".bam":
        continue
    # extract the sample id
    base_name = os.path.basename( fname )
    sample_type = ".".join(base_name.split(".")[:-3])
    samples[ sample_type ].append( os.path.abspath( fname ) )

pad_len = len( "-                   " )

for sample_type, fnames in samples.iteritems():
    for i, fname in enumerate( fnames ):
        line = "\t".join( ( "rnaseq_unstranded_bam", sample_type.ljust( pad_len ), "sample_%i" % (i+1), ".", fname  ) )
        print line
