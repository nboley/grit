import os, sys
import subprocess

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(\
        description='Build transcripts and quantify expression levels from RNAseq, CAGE, and poly(A) assays.')

    parser.add_argument( 
        '--rnaseq-reads', type=argparse.FileType('rb'), required=True, 
        help='BAM file containing mapped RNAseq reads.')
    parser.add_argument( '--rnaseq-read-type', required=True,
        choices=["forward", "backward"],
        help='Whether or not the first RNAseq read in a pair needs to be reversed to be on the correct strand.')
    
    parser.add_argument( '--cage-reads', type=argparse.FileType('rb'),
        help='BAM file containing mapped cage reads.')
    parser.add_argument( '--rampage-reads', type=argparse.FileType('rb'),
        help='BAM file containing mapped rampage reads.')
    
    parser.add_argument( '--polya-reads', type=argparse.FileType('rb'), 
        help='BAM file containing mapped polya reads.')

    parser.add_argument( '--build-bedgraphs', default=False,action='store_true',
                         help='Build read coverage bedgraphs.')
    
    parser.add_argument( '--reference', help='Reference GTF')
    parser.add_argument( '--fasta', type=file,
        help='Fasta file containing the genome sequence - if provided the ORF finder is automatically run.')
    parser.add_argument( '--use-reference-genes', 
                         help='Use genes boundaries from the reference annotation.', 
                         default=False, action='store_true')
    parser.add_argument( '--use-reference-junctions', 
                         help='Include junctions from the reference annotation.',
                         default=False, action='store_true')
    parser.add_argument( '--use-reference-tss', 
                         help='Use TSS\'s taken from the reference annotation.',
                         default=False, action='store_true')
    parser.add_argument( '--use-reference-tes', 
                         help='Use TES\'s taken from the reference annotation.',
                         default=False, action='store_true')
    parser.add_argument( '--use-reference-promoters', 
                         help='Use promoters\'s inferred from the start of reference transcripts.',
                         default=False, action='store_true')
    parser.add_argument( '--use-reference-polyas', 
                         help='Use polya sites inferred from the end of reference transcripts.',
                         default=False, action='store_true')
    
    parser.add_argument( '--ofprefix', '-o', default="discovered",
        help='Output files prefix. (default: discovered)')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    parser.add_argument( '--debug-verbose', default=False, action='store_true',
        help='Whether or not to print debugging information.')
    
    parser.add_argument( '--batch-mode', '-b', 
        default=False, action='store_true',
        help='Disable the ncurses frontend, and just print status messages to stderr.')
    
    parser.add_argument( '--threads', '-t', default=1, type=int,
        help='The number of threads to use.')
        
    args = parser.parse_args()
    
    if None == args.reference and args.use_reference_genes:
        raise ValueError, "--reference must be set if --use-reference-genes is set"
    if None == args.reference and args.use_reference_junctions:
        raise ValueError, "--reference must be set if --use-reference-junctions is set"
    if None == args.reference and args.use_reference_tss:
        raise ValueError, "--reference must be set if --use-reference-tss is set"
    if None == args.reference and args.use_reference_tes:
        raise ValueError, "--reference must be set if --use-reference-tes is set"
    if None == args.reference and args.use_reference_promoters:
        raise ValueError, "--reference must be set if --use-reference-promoters is set"
    if None == args.reference and args.use_reference_polyas:
        raise ValueError, "--reference must be set if --use-reference-polyas is set"
     
    if ( args.cage_reads == None 
         and args.rampage_reads == None 
         and not args.use_reference_tss
         and not args.use_reference_promoters):
        raise ValueError, "--cage-reads or --rampage-reads or --use-reference-tss or --use-reference-promoters must be set"
    if args.cage_reads != None and args.rampage_reads != None:
        raise ValueError, "--cage-reads and --rampage-reads may not both be set"
    
    if ( args.polya_reads == None 
         and not args.use_reference_tes
         and not args.use_reference_polyas ):
        raise ValueError, "Either --polya-reads or --use-reference-tes or --use-reference-polyas must be set"
    
    return args

def run_find_elements( args ):
    elements_ofname = args.ofprefix + ".elements.bed"
    
    command = [ "python", 
                os.path.join( os.path.dirname(__file__), 
                              "..", "grit/find_elements.py" ) ]
    command.extend( ("--rnaseq-reads", args.rnaseq_reads.name) )
    command.extend( ("--rnaseq-read-type", args.rnaseq_read_type) )
    if args.cage_reads != None: 
        command.extend( ("--cage-reads", args.cage_reads.name) )
    if args.rampage_reads != None:
        command.extend( ("--rampage-reads", args.rampage_reads.name) )
    if args.polya_reads != None:
        command.extend( ("--polya-reads", args.polya_reads.name) )
    if args.reference != None: command.extend( ("--reference", args.reference) )
    if args.use_reference_genes: command.append( "--use-reference-genes" )
    if args.use_reference_junctions: command.append("--use-reference-junctions")
    if args.use_reference_tss: command.append("--use-reference-tss")
    if args.use_reference_tes: command.append("--use-reference-tes")
    if args.use_reference_promoters: command.append("--use-reference-promoters")
    if args.use_reference_polyas: command.append("--use-reference-polyas")

    command.extend( ("--ofname", elements_ofname) )
    if args.batch_mode: command.append( "--batch-mode" )
    command.extend( ("--threads", str(args.threads)) )
    
    subprocess.check_call(command)
    
    return  elements_ofname

def run_build_transcripts(args, elements_fname):
    transcripts_ofname = args.ofprefix + ".transcripts.gtf"
    expression_ofname = args.ofprefix + ".isoforms.fpkm_tracking"
    
    command = [ "python", 
                os.path.join( os.path.dirname(__file__), 
                              "..", "grit/build_transcripts.py" ) ]
    command.extend( ("--ofname", transcripts_ofname) )
    command.extend( ("--expression-ofname", expression_ofname) )
    
    command.extend( ("--elements", elements_fname) )    
    
    command.extend( ("--rnaseq-reads", args.rnaseq_reads.name) )
    command.extend( ("--rnaseq-read-type", args.rnaseq_read_type) )
    
    command.append( "--estimate-confidence-bounds" )
    
    if args.cage_reads != None: 
        command.extend( ("--cage-reads", args.cage_reads.name) )
    if args.rampage_reads != None:
        command.extend( ("--rampage-reads", args.rampage_reads.name) )
    
    if args.polya_reads != None:
        command.extend( ("--polya-reads", args.polya_reads.name) )
    
    if args.fasta != None: command.extend( ("--fasta", args.fasta.name) )
    
    if args.batch_mode: command.append( "--batch-mode" )
    command.extend( ("--threads", str(args.threads)) )
    
    subprocess.check_call(command)
    
    return

def run_bam2wig(fname, op_prefix, assay, 
                nthreads, reverse_read_strand, verbose):
    print >> sys.stderr, "Building bedgraph for %s" % fname
    assert assay in ["r", "c", "p"], "Unrecognized assay '%s'" % assay
    command = ["python", os.path.join(os.path.dirname(__file__), "bam2wig.py" )]
    command.extend( ("--mapped-reads-fname", op_prefix ))
    command.extend( ("--out-fname-prefix", fname ))
    command.extend( ("--assay",  assay))
    command.extend( ("--threads",  str(nthreads)))
    if reverse_read_strand:
        command.append( "--reverse-read-strand" )
    if verbose: command.append( "--verbose" )
    subprocess.check_call(command)

def run_all_bam2wigs( args ):
    run_bam2wig(args.rnaseq_reads.name, 
                os.path.basename(args.rnaseq_reads.name),
                'r', args.threads, bool(args.rnaseq_read_type=='backward'),
                args.verbose)
    if args.cage_reads != None:
        run_bam2wig(args.cage_reads.name, 
                    os.path.basename(args.cage_reads.name),
                    'c', args.threads, True, args.verbose)
    if args.rampage_reads != None:
        raise NotImplemented, "Cannot build rampage read coverage bedgraph."
    if args.polya_reads != None:
        run_bam2wig(args.polya_reads.name, 
                    os.path.basename(args.polya_reads.name),
                    'p', args.threads, False, args.verbose)
    
    return

def main():
    args = parse_arguments()
    if args.build_bedgraphs:
        run_all_bam2wigs(args)
    
    elements_fname = run_find_elements(args)
    run_build_transcripts( args, elements_fname )
    
if __name__ == '__main__':
    main()
