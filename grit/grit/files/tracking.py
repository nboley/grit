import os, sys

from collections import namedtuple
from itertools import izip

from reads import fix_chrm_name_for_ucsc

import grit.config

FPKMTrackingLine = namedtuple( 
    'FPKMTrackingLine', 
    ["tracking_id", "class_code", "nearest_ref_id", 
     "gene_id", "gene_short_name", "tss_id", "locus", 
     "length"] )

ExpressionTrackingLine = namedtuple( 
    'ExpressionTrackingLine', 
    ["tracking_id", "gene_id", 
     "coverage", "FPKM", "FPKM_lo", "FPKM_hi", "status"])



def build_gene_tracking_lines(gene):
    rv = []
    contig_name = gene.chrm
    if config.FIX_CHRM_NAMES_FOR_UCSC:
        contig_name = fix_chrm_name_for_ucsc(contig_name)
    
    for (mle, ub, lb, t) in izip( 
            mles, ubs, lbs, gene.transcripts):
        line = FPKMTrackingLine(
            t.id, # tracking ID
            '-' if t.class_code == None else t.class_code,  # class code
            '-' if t.ref_trans == None else t.ref_trans, # nearest ref id
            t.gene_id, # gene unique id
            '-' if t.ref_gene == None else t.ref_gene, # nearest reference gene
            '-', # TSS ID
            "%s:%s:%i-%i" % (contig_name, t.strand, t.start, t.stop), 
            str(t.calc_length()) # transcript length
            )
        rv.append(line)
    
    return rv

def load_expression_tracking_data(fp):
    rv = {}
    for line_num, line in enumerate(fp):
        # skip the header
        if line_num == 0: continue
        data = line.split()
        for i, val in enumerate(data[2:6]):
            if val == '-': val = None
            else: val = float(val)
            data[i+2]  = val
        data = ExpressionTrackingLine(*data)
        rv[data.tracking_id] = data
    
    return rv
