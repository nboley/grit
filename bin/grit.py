import os, sys
import subprocess
from collections import defaultdict, namedtuple
from itertools import chain

ControlFileEntry = namedtuple('ControlFileEntry', [
        'sample_type', 'rep_id', 
        'assay', 'paired', 'stranded', 'read_type', 
        'filename'])

InputData = namedtuple('InputData', [
        'sample_type', 'rep_ids',
        'rnaseq_reads', 'rnaseq_read_types', 'num_mapped_rnaseq_reads',
        'cage_reads', 'rampage_reads', 'polya_reads' ])

def run_find_elements( all_rnaseq_reads, all_rnaseq_read_types, 
                       all_num_mapped_rnaseq_reads,
                       all_cage_reads, all_rampage_reads, all_polya_reads, 
                       ofprefix, args ):
    elements_ofname = ofprefix + ".elements.bed"
    
    assert len(all_cage_reads) <= 1
    assert len(all_rampage_reads) <= 1
    assert len(all_polya_reads) <= 1
    
    command = [ "python", 
                os.path.join( os.path.dirname(__file__), 
                              "..", "grit/find_elements.py" ) ]
    command.extend( chain(("--rnaseq-reads",), (
               rnaseq_reads.name for rnaseq_reads in all_rnaseq_reads) ))
    command.extend( chain(("--rnaseq-read-type",), (
               rnaseq_read_type for rnaseq_read_type in all_rnaseq_read_types)))
    if all_num_mapped_rnaseq_reads != None:
        all_num_mapped_rnaseq_reads = sum(all_num_mapped_rnaseq_reads)
        command.extend( ("--num-mapped-rnaseq-reads",
                         str(num_mapped_rnaseq_reads) ))
    if len(all_cage_reads) > 0:
        command.extend( ("--cage-reads", all_cage_reads[0].name) )
    if len(all_rampage_reads) > 0:
        command.extend( ("--rampage-reads", all_rampage_reads[0].name) )
    if len(all_polya_reads) > 0:
        command.extend( ("--polya-reads", all_polya_reads[0].name) )
    
    if args.reference != None: command.extend( 
        ("--reference", args.reference.name) )
    if args.use_reference_genes: command.append( "--use-reference-genes" )
    if args.use_reference_junctions: command.append("--use-reference-junctions")
    if args.use_reference_tss: command.append("--use-reference-tss")
    if args.use_reference_tes: command.append("--use-reference-tes")
    if args.use_reference_promoters: command.append("--use-reference-promoters")
    if args.use_reference_polyas: command.append("--use-reference-polyas")

    if args.ucsc: command.append("--ucsc")

    command.extend( ("--ofname", elements_ofname) )
    if args.batch_mode: command.append( "--batch-mode" )
    if args.region != None: command.append( "--region %s" % args.region )
    command.extend( ("--threads", str(args.threads)) )
    if args.verbose: command.append( "--verbose" )
    
    subprocess.check_call(command)
    
    return elements_ofname

def run_build_transcripts(elements_fname, ofprefix,
                          rnaseq_reads,rnaseq_read_type,num_mapped_rnaseq_reads,
                          cage_reads, rampage_reads, polya_reads, args):
    transcripts_ofname = ofprefix + ".transcripts.gtf"
    expression_ofname = ofprefix + ".isoforms.fpkm_tracking"
    
    command = [ "python", 
                os.path.join( os.path.dirname(__file__), 
                              "..", "grit/build_transcripts.py" ) ]
    command.extend( ("--ofname", transcripts_ofname) )
    command.extend( ("--expression-ofname", expression_ofname) )
    
    command.extend( ("--elements", elements_fname) )    
    
    command.extend( ("--rnaseq-reads", rnaseq_reads.name) )
    command.extend( ("--rnaseq-read-type", rnaseq_read_type) )
    if num_mapped_rnaseq_reads != None: 
        command.extend( ("--num-mapped-rnaseq-reads",
                         str(num_mapped_rnaseq_reads) ))
    
    command.append( "--estimate-confidence-bounds" )
    
    if cage_reads != None: 
        command.extend( ("--cage-reads", cage_reads.name) )
    if rampage_reads != None:
        command.extend( ("--rampage-reads", rampage_reads.name) )
    
    if polya_reads != None:
        command.extend( ("--polya-reads", polya_reads.name) )
    
    if args.fasta != None: command.extend( ("--fasta", args.fasta.name) )
    
    if args.batch_mode: command.append( "--batch-mode" )
    command.extend( ("--threads", str(args.threads)) )
    if args.verbose: command.append( "--verbose" )

    if args.ucsc: command.append("--ucsc")
    
    subprocess.check_call(command)
    
    return

def run_bam2wig(fname, op_prefix, assay, 
                nthreads, reverse_read_strand, verbose):
    print >> sys.stderr, "Building bedgraph for %s" % fname
    assert assay in ["rnaseq", "cage", "rampage", "polya"], \
        "Unrecognized assay '%s'" % assay
    command = ["python", os.path.join(os.path.dirname(__file__), "bam2wig.py" )]
    command.extend( ("--mapped-reads-fname", op_prefix ))
    command.extend( ("--out-fname-prefix", fname ))
    command.extend( ("--assay",  assay))
    command.extend( ("--threads",  str(nthreads)))
    if reverse_read_strand:
        command.append( "--reverse-read-strand" )
    if verbose: command.append( "--verbose" )
    subprocess.check_call(command)

def run_all_bam2wigs( all_rnaseq_reads, all_rnaseq_read_types, 
                      all_cage_reads, all_rampage_reads, all_polya_reads, 
                      args ):
    for rnaseq_reads, rnaseq_reads_type in zip(
            all_rnaseq_reads, all_rnaseq_read_types):
        run_bam2wig(rnaseq_reads.name, os.path.basename(rnaseq_reads.name),
                    'rnaseq', args.threads, bool(rnaseq_reads_type=='backward'),
                    args.verbose)
    for rampage_reads in all_rampage_reads:
        run_bam2wig(rampage_reads.name, os.path.basename(rampage_reads.name),
                    'rampage', args.threads, True, args.verbose)    
    for cage_reads in all_cage_reads:
        run_bam2wig(cage_reads.name, os.path.basename(cage_reads.name),
                    'cage', args.threads, True, args.verbose)    
    if args.polya_reads != None:
        run_bam2wig(args.polya_reads.name, 
                    os.path.basename(args.polya_reads.name),
                    'polya', args.threads, False, args.verbose)
    
    return

def truth_value(string):
    if string.lower() in ('true', '1', 't'): return True
    elif string.lower() in ('false', '0', 'f'): return False
    else: raise ValueError, "Unrecognized boolean '%s'" % string

def parse_control_file(control_fp):
    samples = defaultdict(list)
    for line in control_fp:
        if line.strip().startswith("#"): continue
        if line.strip() == '': continue
        line = ControlFileEntry(*(line.split()))
        samples[line.sample_type].append( line )
        
    rv = []
    for sample_type, lines in samples.iteritems():
        rnaseq_reads = []
        rnaseq_read_types = []
        num_mapped_rnaseq_reads = None
        cage_reads = []
        rampage_reads = []
        polya_reads = []
        rep_ids = []
        for line in lines:
            if line.assay == 'rnaseq':
                rnaseq_reads.append( open(line.filename) )
                assert truth_value(line.paired) == True
                assert truth_value(line.stranded) == True
                rnaseq_read_types.append(line.read_type)
                rep_ids.append(line.rep_id)
            elif line.assay == 'cage':
                assert truth_value(line.paired) == False
                assert truth_value(line.stranded) == True
                assert line.read_type == 'backward'
                cage_reads.append(open(line.filename))
            elif line.assay == 'rampage':
                raise NotImplemented, "RAMPAGE not implemented in control file"
            elif line.assay == 'polya':
                assert truth_value(line.paired) == False
                assert truth_value(line.stranded) == True
                assert line.read_type == 'backward'
                polya_reads.append(open(line.filename))
        rv.append( InputData( sample_type, rep_ids,
                              rnaseq_reads, rnaseq_read_types, 
                              num_mapped_rnaseq_reads,
                              cage_reads, rampage_reads, polya_reads ) )
    
    return rv

def verify_args_are_sufficient(args, rnaseq_reads, 
                               cage_reads, rampage_reads, polya_reads):
    if ( len(cage_reads) == 0 
         and len(rampage_reads) == 0
         and not args.use_reference_tss
         and not args.use_reference_promoters):
        raise ValueError, "--cage-reads or --rampage-reads or --use-reference-tss or --use-reference-promoters must be set"
    
    if len(cage_reads) > 0 and len(rampage_reads) > 0:
        raise ValueError, "--cage-reads and --rampage-reads may not both be set"
    
    if ( len(polya_reads) == 0
         and not args.use_reference_tes
         and not args.use_reference_polyas ):
        raise ValueError, "Either --polya-reads or --use-reference-tes or --use-reference-polyas must be set"

    return

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='Build transcripts and quantify expression levels from RNAseq, CAGE, and poly(A) assays.')

    parser.add_argument( '--control', type=file, 
        help='GRIT control file. Allows better control over the types of input files.')

    parser.add_argument( '--rnaseq-reads', type=argparse.FileType('rb'), 
        help='BAM file containing mapped RNAseq reads.')
    parser.add_argument( '--rnaseq-read-type', choices=["forward", "backward"],
        help='Whether or not the first RNAseq read in a pair needs to be reversed to be on the correct strand.')
    parser.add_argument( '--num-mapped-rnaseq-reads', type=int,
        help="The total number of mapped rnaseq reads ( needed to calculate the FPKM ). This only needs to be set if it isn't found by a call to samtools idxstats." )
    
    parser.add_argument( '--cage-reads', type=argparse.FileType('rb'),
        help='BAM file containing mapped cage reads.')
    parser.add_argument( '--rampage-reads', type=argparse.FileType('rb'),
        help='BAM file containing mapped rampage reads.')
    
    parser.add_argument( '--polya-reads', type=argparse.FileType('rb'), 
        help='BAM file containing mapped polya reads.')

    parser.add_argument( '--build-bedgraphs', default=False,action='store_true',
                         help='Build read coverage bedgraphs.')
    
    parser.add_argument( '--reference', help='Reference GTF', type=file)
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

    parser.add_argument( '--ucsc', default=False, action='store_true',
        help='Try to format contig names in the ucsc format (typically by prepending a chr).')    
    parser.add_argument( '--batch-mode', '-b', 
        default=False, action='store_true',
        help='Disable the ncurses frontend, and just print status messages to stderr.')
    parser.add_argument( '--region', 
        help='Only use the specified region ( currently only accepts a contig name ).')
    
    parser.add_argument( '--threads', '-t', default=1, type=int,
        help='The number of threads to use.')
        
    args = parser.parse_args()

    if None == args.control and None == args.rnaseq_reads:
        raise ValueError, "--control or --rnaseq-reads must be set"
    
    #if args.rampage_reads != None:
    #    raise NotImplemented, "RAMPAGE is not currently implemented in grit.py."
    
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

    if args.control != None:
        assert args.num_mapped_rnaseq_reads == None
    
    return args

def main():
    args = parse_arguments()
    
    # build the element sets. If there is no control file, then 
    # we simply take the passed in files. If there is a control
    # file, then build the sets
    if args.control == None:
        all_input_data = [ InputData(
                None, None,
                (args.rnaseq_reads,),
                (args.rnaseq_read_type,),
                (args.num_mapped_rnaseq_reads,) \
                    if args.num_mapped_rnaseq_reads != None else None,
                (args.cage_reads,) if args.cage_reads != None else (),
                (args.rampage_reads,) if args.rampage_reads != None else (),
                (args.polya_reads,) if args.polya_reads != None else (),
        ),]
    else:
        all_input_data = parse_control_file(args.control)
    
    if args.build_bedgraphs:
        for input_data in all_input_data:
            run_all_bam2wigs(input_data.rnaseq_reads, 
                             input_data.rnaseq_read_types, 
                             input_data.cage_reads, 
                             input_data.rampage_reads, 
                             input_data.polya_reads, 
                             args)

    for ipd in all_input_data:
            verify_args_are_sufficient( args, ipd.rnaseq_reads, ipd.cage_reads, 
                                        ipd.rampage_reads, ipd.polya_reads )
            ofprefix = args.ofprefix
            if ipd.sample_type != None: ofprefix += ("." + ipd.sample_type)

            elements_fname = run_find_elements( 
                ipd.rnaseq_reads, ipd.rnaseq_read_types, 
                ipd.num_mapped_rnaseq_reads,
                ipd.cage_reads, ipd.rampage_reads, ipd.polya_reads, 
                ofprefix, args)

            for i, (rnaseq_reads, rnaseq_read_type) in enumerate(zip(
                    ipd.rnaseq_reads, ipd.rnaseq_read_types)):
                if ipd.num_mapped_rnaseq_reads != None:
                    num_mapped_rnaseq_reads = ipd.num_mapped_rnaseq_reads[i]
                else: num_mapped_rnaseq_reads = None
                
                assert len(ipd.cage_reads) <= 1
                cage_reads = ipd.cage_reads[0] \
                    if len(ipd.cage_reads) == 1 else None
                assert len(ipd.rampage_reads) <= 1
                rampage_reads = ipd.rampage_reads[0] \
                    if len(ipd.rampage_reads) == 1 else None
                assert len(ipd.polya_reads) <= 1
                polya_reads = ipd.polya_reads[0] \
                    if len(ipd.polya_reads) == 1 else None
                
                ofprefix = args.ofprefix
                if ipd.sample_type != None: ofprefix += ("." + ipd.sample_type)
                if ipd.rep_ids != None: ofprefix += ("." + ipd.rep_ids[i])
                
                run_build_transcripts( 
                    elements_fname, ofprefix, 
                    rnaseq_reads, rnaseq_read_type, 
                    num_mapped_rnaseq_reads,
                    cage_reads, rampage_reads, polya_reads, args)

if __name__ == '__main__':
    main()
