# Copyright (c) 2011-2012 Nathan Boley

import sys


myfid = open(sys.argv[1])
cr = sys.argv[2:]

def parse_gtf( fid, ch, st, en, stran ):
    all_names = dict()
    for line in fid:       
        data = line.strip().split('\t')
        chrm = data[0]
        #print data
        strand = data[6]
        if not chrm == ch or not strand == stran:
            continue
        start = int(data[3])
        if start > en:
            continue        
        stop = int(data[4])
        if stop < start:
            continue
        gene_name = data[8].split(' ')[1][1:-2]
        all_names[gene_name] = ''

    fid.seek(0)
    for line in fid:
        data = line.strip().split('\t')
        gene_name = data[8].split(' ')[1][1:-2]
        if all_names.has_key(gene_name):
            print line,


parse_gtf(myfid, cr[0],int(cr[1]),int(cr[2]),cr[3])

    
