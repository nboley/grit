import argparse

VEBOSE = False
log_statement = None
from files.gtf import parse_gtf_line, load_gtf
from files.reads import MergedReads, RNAseqReads, CAGEReads, RAMPAGEReads, PolyAReads

def load_polya_reads( polya_bams, polya_strands_need_to_be_reversed, 
                      ref_genes=None ):
    if VERBOSE: log_statement( 'Loading polyA reads bams' )        
    if polya_bams == None: return None

    all_reads = []
    for i, polya_bam_fp in enumerate(polya_bams):
        reads = PolyAReads(polya_bam_fp.name)
        rev_reads = None if None == polya_strands_need_to_be_reversed else\
            polya_strands_need_to_be_reversed[i]
        reads.init(pairs_are_opp_strand=True,
                   reverse_read_strand=rev_reads, 
                   ref_genes=ref_genes)
        all_reads.append(reads)
    
    return MergedReads(all_reads)

def load_promoter_reads(cage_bams, cage_strands_need_to_be_reversed,
                        rampage_bams, rampage_strands_need_to_be_reversed,
                        ref_genes = None):
    assert cage_bams == None or rampage_bams == None, \
        "Can not use both RAMPAGE and CAGE reads"
    if cage_bams != None: 
        bams = cage_bams
        reads_class = CAGEReads
        strands_need_to_be_reversed = cage_strands_need_to_be_reversed
    elif rampage_bams != None: 
        bams =rampage_bams
        reads_class = RAMPAGEReads
        strands_need_to_be_reversed = rampage_strands_need_to_be_reversed
    else: return None

    all_reads = []
    for i, bam_fp in enumerate(bams):
        reads = reads_class(bam_fp.name)
        rev_reads = None if strands_need_to_be_reversed == None else \
            strands_need_to_be_reversed[i]
        reads.init(reverse_read_strand=rev_reads, ref_genes=ref_genes)
        all_reads.append(reads)
    
    return MergedReads(all_reads)

def load_rnaseq_reads(rnaseq_bams, rnaseq_strands_need_to_be_reversed, 
                      ref_genes=None):
    if VERBOSE: log_statement( 'Loading RNAseq read bams' )
    all_reads = []
    for i, rnaseq_bam_fp in enumerate(rnaseq_bams):
        reads = RNAseqReads(rnaseq_bam_fp.name)
        rev_reads = None if None == rnaseq_strands_need_to_be_reversed else\
            rnaseq_strands_need_to_be_reversed[i]
        reads.init(reverse_read_strand=rev_reads, ref_genes=ref_genes)
        all_reads.append(reads)
    
    rnaseq_reads = MergedReads(all_reads)
    
    return rnaseq_reads

def initialize_reads_from_args(args):
    """Initialize the various reads objects.
    
    Given an argparse arg object, parse the argument, do error checking, 
    and then load and initialize the read objects. We put this in lib because
    find_elements and build_transcripts needs the same code. 
    
    returns: promoter_reads, rnaseq_reads, polya_reads, genes
    """
    if ((    args.cage_read_type == None 
            or any(x=='auto' for x in args.cage_read_type) )
        and None == args.reference ):
        raise ValueError, "--reference must be set if --cage-read-type is not set or set to 'auto' (GRIT needs an annotation to determine the read type)"
    
    if (( args.rampage_read_type == None 
            or any(x=='auto' for x in args.rampage_read_type) )
        and None == args.reference ):
        raise ValueError, "--reference must be set if --rampage-read-type is not set or set to 'auto' (GRIT needs an annotation to determine the read type)"

    if (( args.rnaseq_read_type == None 
            or any(x=='auto' for x in args.rnaseq_read_type) )
        and None == args.reference ):
        raise ValueError, "--reference must be set if --rnaseq-read-type is not set or set to 'auto' (GRIT needs an annotation to determine the read type)"

    if (( args.polya_read_type == None 
            or any(x=='auto' for x in args.polya_read_type) )
        and None == args.reference ):
        raise ValueError, "--reference must be set if --polya-read-type is not set or set to 'auto' (GRIT needs an annotation to determine the read type)"

    if args.cage_read_type == None:
        cage_strands_need_to_be_reversed = ['auto']*(0 if args.cage_reads == None else len(args.cage_reads))
    else:
        cage_strands_need_to_be_reversed = [ 
            bool(read_type.lower() == 'backward')
            for read_type in args.cage_read_type ]
    
    if args.rampage_read_type == None:
        rampage_strands_need_to_be_reversed = ['auto']*(0 if args.rampage_reads == None else len(args.rampage_reads))
    else:
        rampage_strands_need_to_be_reversed = [ 
            bool(read_type.lower() == 'backward')
            for read_type in args.rampage_read_type ]

    if args.rnaseq_read_type == None:
        rnaseq_strands_need_to_be_reversed = ['auto']*len(args.rnaseq_reads)
    else:
        rnaseq_strands_need_to_be_reversed = [ 
            bool(read_type.lower() == 'backward')
            for read_type in args.rnaseq_read_type ]

    if args.polya_read_type == None:
        polya_strands_need_to_be_reversed = ['auto']*len(args.rnaseq_reads)
    else:
        polya_strands_need_to_be_reversed = [ 
            bool(read_type.lower() == 'backward')
            for read_type in args.polya_read_type ]

    if args.reference != None:
        if VERBOSE: log_statement("Loading annotation file.")
        ref_genes = load_gtf( args.reference )
    else:
        ref_genes = []

    # load and initialize the rnaseq read bam, including the global
    # TOTAL_MAPPED_READS variable
    rnaseq_reads = load_rnaseq_reads(
        args.rnaseq_reads, rnaseq_strands_need_to_be_reversed, ref_genes)

    if VERBOSE: log_statement( 'Loading promoter reads bams' )        
    promoter_reads = load_promoter_reads(
        args.cage_reads, cage_strands_need_to_be_reversed,
        args.rampage_reads, rampage_strands_need_to_be_reversed,
        ref_genes=ref_genes)

    polya_reads = load_polya_reads(
        args.polya_reads, polya_strands_need_to_be_reversed, ref_genes=ref_genes)

    return promoter_reads, rnaseq_reads, polya_reads, ref_genes
