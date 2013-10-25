import os, sys
import subprocess
from collections import defaultdict, namedtuple
from itertools import chain
import sqlite3

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
        
    command = [ "python", 
                os.path.join( os.path.dirname(__file__), 
                              "..", "grit/find_elements.py" ) ]
    command.extend( chain(("--rnaseq-reads",), (
               rnaseq_reads for rnaseq_reads in all_rnaseq_reads) ))
    command.extend( chain(("--rnaseq-read-type",), (
               rnaseq_read_type for rnaseq_read_type in all_rnaseq_read_types)))
    if all_num_mapped_rnaseq_reads != None:
        all_num_mapped_rnaseq_reads = sum(all_num_mapped_rnaseq_reads)
        command.extend( ("--num-mapped-rnaseq-reads",
                         str(num_mapped_rnaseq_reads) ))
    if len(all_cage_reads) > 0:
        command.extend( chain(("--cage-reads",), (
               cage_reads for cage_reads in all_cage_reads) ) )
    if len(all_rampage_reads) > 0:
        command.extend( chain(("--rampage-reads",), (
               rampage_reads for rampage_reads in all_rampage_reads) ) )
    if len(all_polya_reads) > 0:
        command.extend( chain(("--polya-reads",), (
               polya_reads for polya_reads in all_polya_reads) ) )
    
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
    if args.region != None: command.extend( ("--region", "%s" % args.region) )
    command.extend( ("--threads", str(args.threads)) )
    if args.verbose: command.append( "--verbose" )

    #print " ".join(command)
    subprocess.check_call(command)
    
    return elements_ofname

def run_build_transcripts(elements_fname, ofprefix,
                          rnaseq_reads, rnaseq_read_types, 
                          num_mapped_rnaseq_reads,
                          cage_reads, rampage_reads, polya_reads, args):
    transcripts_ofname = ofprefix + ".transcripts.gtf"
    expression_ofname = ofprefix + ".isoforms.fpkm_tracking"
    
    assert len(rnaseq_reads) == 1
    assert len(rnaseq_read_types) == 1
    assert len(cage_reads) <= 1
    assert len(rampage_reads) <= 1
    assert len(polya_reads) <= 1
    
    command = [ "python", 
                os.path.join( os.path.dirname(__file__), 
                              "..", "grit/build_transcripts.py" ) ]
    command.extend( ("--ofname", transcripts_ofname) )
    command.extend( ("--expression-ofname", expression_ofname) )
    
    command.extend( ("--elements", elements_fname) )    
    

    if args.only_build_candidate_transcripts:
        command.append( "--only-build-candidate-transcripts" )
    else:
        command.extend( ("--rnaseq-reads", rnaseq_reads[0]) )
        command.extend( ("--rnaseq-read-type", rnaseq_read_types[0]) )
        if num_mapped_rnaseq_reads != None: 
            command.extend( ("--num-mapped-rnaseq-reads",
                             str(num_mapped_rnaseq_reads) ))    
        if len(cage_reads) > 0: 
            command.extend( ("--cage-reads", cage_reads[0]) )
        if len(rampage_reads) > 0:
            command.extend( ("--rampage-reads", rampage_reads[0]) )

        if len(polya_reads) > 0:
            command.extend( ("--polya-reads", polya_reads[0]) )
        
        command.append( "--estimate-confidence-bounds" )

    
    if args.fasta != None: command.extend( ("--fasta", args.fasta.name) )
    
    if args.batch_mode: command.append( "--batch-mode" )
    command.extend( ("--threads", str(args.threads)) )
    if args.verbose: command.append( "--verbose" )

    if args.ucsc: command.append("--ucsc")
    
    #print " ".join(command)
    subprocess.check_call(command)
    
    return

def run_bam2wig(fname, op_prefix, assay, region,
                nthreads, reverse_read_strand, verbose):
    print >> sys.stderr, "Building bedgraph for %s" % fname
    assert assay in ["rnaseq", "cage", "rampage", "polya"], \
        "Unrecognized assay '%s'" % assay
    command = ["python", os.path.join(os.path.dirname(__file__), "bam2wig.py" )]
    command.extend( ("--mapped-reads-fname", fname ))
    command.extend( ("--out-fname-prefix", op_prefix ))
    command.extend( ("--assay",  assay))
    command.extend( ("--threads",  str(nthreads)))
    if reverse_read_strand:
        command.append( "--reverse-read-strand" )
    if verbose: command.append( "--verbose" )
    if region != None: command.extend( ("--region", "%s" % region) )
    subprocess.check_call(command)

def run_all_bam2wigs(conn, args):
    for rnaseq_reads, rnaseq_reads_type in get_elements( 
            conn, ('filename', 'read_type'), 'rnaseq'):
        run_bam2wig(rnaseq_reads, os.path.basename(rnaseq_reads),
                    'rnaseq', args.region,
                    args.threads, bool(rnaseq_reads_type=='backward'),
                    args.verbose)
    for data_type in ('cage', 'rampage', 'polya'):
        for reads, in get_elements( conn, ('filename', ), data_type):
            run_bam2wig(reads, os.path.basename(reads),
                        data_type, args.region,
                        args.threads, False, args.verbose)
    return

def truth_value(string):
    if string.lower() in ('true', '1', 't'): return True
    elif string.lower() in ('false', '0', 'f'): return False
    else: raise ValueError, "Unrecognized boolean '%s'" % string


def get_elements( conn, column_names, assay, sample_type=None, rep_id=None ):
    """Get values of the specified column name, optionally filtered by
       sample_type and rep_id
    """
    query = "SELECT group_concat({}, ' ')".format(column_names[0])
    for column_name in column_names[1:]:
        query += ", group_concat({}, ' ')".format(column_name)
    query += " FROM data WHERE assay='{}' ".format(assay)
    
    if rep_id != None: assert sample_type != None, "rep_id can't be filtered without a sample type filter"
    if sample_type != None:
        query += " AND (sample_type = '{}' OR sample_type = '*') ".format(
            sample_type)
    if rep_id != None:
        query += " AND (rep_id = '{}' OR rep_id = '*') ".format(rep_id)
    with conn:
        return conn.execute( query  ).fetchall()
    
def initialize_sample_db():
    conn = sqlite3.connect(':memory:')
    with conn:
        conn.execute("""
        CREATE TABLE data (
           sample_type text,
           rep_id text,
           assay text,
           paired text,
           stranded text,
           read_type text,
           filename text
        );""")
    return conn

def parse_control_file(conn, control_fp):
    lines = []
    for line in control_fp:
        if line.strip().startswith("#"): continue
        if line.strip() == '': continue
        lines.append( ControlFileEntry(*(line.split())) )
    return lines

def verify_args_are_sufficient(args, rnaseq_reads, 
                               cage_reads, rampage_reads, polya_reads):
    if ( cage_reads == None
         and rampage_reads == None
         and not args.use_reference_tss
         and not args.use_reference_promoters):
        raise ValueError, "--cage-reads or --rampage-reads or --use-reference-tss or --use-reference-promoters must be set"
    
    if cage_reads != None and rampage_reads != None:
        raise ValueError, "--cage-reads and --rampage-reads may not both be set"
    
    if ( polya_reads == None
         and not args.use_reference_tes
         and not args.use_reference_polyas ):
        raise ValueError, "Either --polya-reads or --use-reference-tes or --use-reference-polyas must be set"

    return

def get_run_data(conn, args, sample_type, rep_id=None):
    rnaseq_reads, rnaseq_read_types = get_elements( 
        conn, ('filename', 'read_type'), 'rnaseq', sample_type, rep_id )[0]
    cage_reads, = get_elements( 
        conn, ('filename',), 'cage', sample_type, rep_id )[0]
    rampage_reads, = get_elements( 
        conn, ('filename',), 'rampage', sample_type, rep_id )[0]
    polya_reads, = get_elements( 
        conn, ('filename',), 'polya', sample_type, rep_id )[0]

    verify_args_are_sufficient( args, rnaseq_reads, cage_reads, 
                                rampage_reads, polya_reads )
    return ( rnaseq_reads.split(" ") if rnaseq_reads != None else [], 
             rnaseq_read_types.split(" ") if rnaseq_read_types != None else [], 
             cage_reads.split(" ") if cage_reads != None else [], 
             rampage_reads.split(" ") if rampage_reads != None else [], 
             polya_reads.split(" ") if polya_reads != None else [],  )

def parse_single_sample_args(args):
    """Parse read data passed in as arguments.

    """
    lines = []
    lines.append( ControlFileEntry( 
            None, None, "rnaseq", True, True, 
            args.rnaseq_read_type, args.rnaseq_reads.name ) )
    if args.cage_reads != None:
        lines.append( ControlFileEntry( 
                None, None, "cage", True, True, 'backward',
                args.cage_reads.name ) )
    if args.rampage_reads != None:
        lines.append( ControlFileEntry( 
                None, None, "rampage", True, True, 'backward',
                args.rampage_reads.name ) )
    if args.polya_reads != None:
        lines.append( ControlFileEntry( 
                None, None, "polya", True, True, 'backward',
                args.polya_reads.name ) )
    return lines

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

    parser.add_argument( '--only-build-candidate-transcripts', 
                         help='Do not estiamte transcript frequencies - just build trnascripts.',
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

    if None != args.num_mapped_rnaseq_reads:
        assert args.control == None, \
            "Can't currently provide read counts with a control file. Try adding them into the BAM header."
    
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
    conn = initialize_sample_db()
    if args.control == None:
        lines = parse_single_sample_args( args )
    else:
        lines = parse_control_file(conn, args.control)
    # insert the various data sources into the database
    with conn:
        conn.executemany( """
        INSERT INTO data VALUES(?, ?, ?, ?, ?, ?, ?)
        """, lines )
    
    if args.build_bedgraphs:
        run_all_bam2wigs(conn, args)
    
    # get the set of distinct sample types
    query = "SELECT DISTINCT sample_type FROM data"
    sample_types = [x[0] for x in conn.execute(query).fetchall()]
    for sample_type in sample_types:
        ofprefix = args.ofprefix
        if sample_type != None: ofprefix += ("." + sample_type)
        
        ( rnaseq_reads, rnaseq_read_types, cage_reads, rampage_reads,
          polya_reads ) = get_run_data(conn, args, sample_type)
        elements_fname = run_find_elements( 
            rnaseq_reads, rnaseq_read_types, args.num_mapped_rnaseq_reads,
            cage_reads, rampage_reads, polya_reads,
            ofprefix, args)
        
        # if we used a control file, and thus have sample types, then find 
        # the unqiue repids for this sample
        if sample_type != None:
            query = "SELECT DISTINCT rep_id FROM data \
                     WHERE sample_type = '{}' AND rep_id != '*'"
            rep_ids = [ x[0] for x in 
                        conn.execute(query.format(sample_type)).fetchall()]
        # otherwise, use everything by setting hte rep id to None
        else: rep_ids = [None,]
        
        for rep_id in rep_ids:
            ( rnaseq_reads, rnaseq_read_types, cage_reads, rampage_reads,
              polya_reads ) = get_run_data(conn, args, sample_type, rep_id)

            ofprefix = args.ofprefix
            if sample_type != None: ofprefix += ("." + sample_type)
            if rep_id != None: ofprefix += ("." + rep_id)
            
            run_build_transcripts( 
                elements_fname, ofprefix, 
                rnaseq_reads, rnaseq_read_types, 
                args.num_mapped_rnaseq_reads,
                cage_reads, rampage_reads, polya_reads, args)

if __name__ == '__main__':
    main()
