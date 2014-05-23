#!/usr/bin/python
import os, sys
import subprocess
from collections import defaultdict, namedtuple
from itertools import chain
import sqlite3
import cPickle as pickle

sys.path.insert(0, "/home/nboley/grit/grit/")

from grit.files.gtf import load_gtf
from grit.files.reads import (
    MergedReads, clean_chr_name,
    RNAseqReads, CAGEReads, RAMPAGEReads, PolyAReads )

from grit.lib.logging import Logger

import grit.find_elements
import grit.build_transcripts
import grit.estimate_transcript_expression
import grit.frag_len

import grit.config as config

ControlFileEntry = namedtuple('ControlFileEntry', [
        'sample_type', 'rep_id', 
        'assay', 'paired', 'stranded', 'read_type', 
        'filename'])

import numpy

def convert_elements_to_arrays(all_elements):
    # convert into array
    all_array_elements = defaultdict( 
        lambda: defaultdict(lambda: numpy.zeros(0)) )
    for key, elements in all_elements.iteritems():
        for element_type, contig_elements in elements.iteritems():
            all_array_elements[key][element_type] \
                = numpy.array( sorted( contig_elements ) )

    return all_array_elements

def load_elements( fp ):
    all_elements = defaultdict( lambda: defaultdict(set) )
    for line in fp:
        if line.startswith( 'track' ): continue
        chrm, start, stop, element_type, score, strand = line.split()[:6]
        # subtract 1 from stop becausee beds are closed open, and we 
        # wnat everything in 0-based closed-closed
        all_elements[(chrm, strand)][element_type].add( 
            (int(start), int(stop)-1) )
    
    return convert_elements_to_arrays(all_elements)

def extract_elements_from_genes( genes ):
    all_elements = defaultdict( lambda: defaultdict(set) )
    for gene in genes:
        for key, val in gene.extract_elements().iteritems():
            all_elements[(gene.chrm, gene.strand)][key].update(val)

    
    return convert_elements_to_arrays( all_elements )

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
            # unfortunately, sqllite returns none for an empty query
            if reads == None: continue
            run_bam2wig(reads, os.path.basename(reads),
                        data_type, args.region,
                        args.threads, False, args.verbose)
    return

class Samples(object):
    """Store and retrieve sample information.

    """
    def parse_control_file(self, control_fp):
        lines = []
        for line in control_fp:
            if line.strip().startswith("#"): continue
            if line.strip() == '': continue
            lines.append( ControlFileEntry(*(line.split())) )
        return lines
    
    def parse_single_sample_args(self, args):
        """Parse read data passed in as arguments.

        """
        lines = []
        lines.append( ControlFileEntry( 
                None, None, "rnaseq", 'true', 'true', 
                args.rnaseq_read_type, args.rnaseq_reads.name ) )
        if args.cage_reads != None:
            lines.append( ControlFileEntry( 
                    None, None, "cage", 'false', 'true', 
                    args.cage_read_type, args.cage_reads.name ) )
        if args.rampage_reads != None:
            lines.append( ControlFileEntry( 
                    None, None, "rampage", 'true', 'true', 
                    args.rampage_read_type, args.rampage_reads.name ) )
        if args.polya_reads != None:
            lines.append( ControlFileEntry( 
                    None, None, "polya", 'false', 'true',
                    args.polya_read_type, args.polya_reads.name ) )
        return lines
    
    def verify_args_are_sufficient(self, rnaseq_reads, promoter_reads, polya_reads):
        if ( len(promoter_reads) == 0
             and not self.args.use_reference_tss_exons
             and not self.args.use_reference_promoters):
            raise ValueError, "Either (cage reads or rampage reads) must be provided for each sample or (--use-reference-tss-exons or --use-reference-promoters) must be set"
        
        if ( len(polya_reads) == 0
             and not self.args.use_reference_tes_exons
             and not self.args.use_reference_polyas ):
            raise ValueError, "Either polya-reads must be provided or (--use-reference-tes-exons or --use-reference-polyas) must be set"
        
        return
    
    def initialize_sample_db(self):
        self.conn = sqlite3.connect(':memory:')
        with self.conn:
            self.conn.execute("""
            CREATE TABLE data (
               sample_type text,
               rep_id text,
               assay text,
               paired text,
               stranded text,
               read_type text,
               filename text
            );""")
    
    def __init__(self, args):
        # the args object returned by parse arguments
        self.args = args
        # cache mapped reads objects that have already been initialized
        self.mapped_reads_cache = {}
        # store parsed reference genes, if necessary
        self.ref_genes = None
        # initialize a sqlite db to store samples
        self.initialize_sample_db()
        # parse the control file, if it exists
        if args.control != None:
            self.control_entries = self.parse_control_file(args.control)
        # if there any reads (bnut no control file from above) then 
        # load them
        elif (args.rnaseq_reads != None or args.polya_reads != None \
                  or args.rampage_reads != None or args.cage_reads != None):
            self.control_entries = self.parse_single_sample_args( args )
        # if there are no provided reads (ie, if we're jsut building transcripts
        # from elements) then there is nothing to load
        else: 
            self.control_entries = []
        
        # if any of the read_type arguments are 'auto', then load the
        # reference genome
        if any(x.read_type == 'auto' for x in self.control_entries):
            if args.reference == None:
                raise ValueError, "One of the read_type entries is set to 'auto' but a reference was not provided"
            if config.VERBOSE: config.log_statement("Loading annotation file.")
            self.ref_genes = load_gtf( args.reference )
        
        # insert the various data sources into the database
        with self.conn:
            self.conn.executemany( """
            INSERT INTO data VALUES(?, ?, ?, ?, ?, ?, ?)
            """, self.control_entries )
        
        return
    
    def __str__(self):
        header = "#%s\n#\n" % ("\t".join([
                    'sample_type', 'rep_id', 'assay', 
                    'paired', 'stranded', 'read_type', 'filename']))
        return header + "\n".join("\t".join(x) for x in self.control_entries)

    def get_elements( self, assay, sample_type=None, rep_id=None, 
                      include_merged=False ):
        """Get values of the specified column name, optionally filtered by
           sample_type and rep_id
        """
        # get the data
        query = "SELECT * FROM data WHERE assay='{a}'".format(a=assay)
        
        if rep_id != None: 
            assert sample_type != None, \
                "rep_id can't be filtered without a sample type filter"
        if sample_type != None:
            query += " AND (sample_type = '{s}'".format(s=sample_type)
            if include_merged: query += "OR sample_type = '*') "
            else: query += ") "
        if rep_id != None:
            query += " AND (rep_id = '{r}' ".format(r=rep_id)
            if include_merged: query += "OR rep_id = '*') "
            else: query += ") "
        
        with self.conn:
            return [ ControlFileEntry(*x) for x in 
                     self.conn.execute( query  ).fetchall() ]

    def _load_rnaseq_reads(self, sample_type, rep_id, include_merged):
        if config.VERBOSE: config.log_statement(
            "Initializing RNAseq reads for sample '%s'." % sample_type)
        all_reads = []
        for data in self.get_elements( 
                'rnaseq', sample_type, rep_id, include_merged=include_merged):
            
            if data.filename in self.mapped_reads_cache:
                reads = self.mapped_reads_cache[data.filename]
                reads.reload()
            else:
                assert data.paired == 'true', "RNASeq reads must be paired"
                assert data.stranded == 'true', "RNASeq reads must be stranded"
                rev_reads = {'forward':False, 'backward':True, 'auto': None}[
                    data.read_type]
                reads = RNAseqReads(data.filename)
                reads.init(reverse_read_strand=rev_reads, ref_genes=self.ref_genes)
                self.mapped_reads_cache[data.filename] = reads
            all_reads.append(reads)
        
        return all_reads

    def _load_promoter_reads(self, sample_type, rep_id, include_merged):
        if config.VERBOSE: config.log_statement(
            "Initializing promoter reads for sample '%s'." % sample_type)
        cage_elements = self.get_elements( 
            'cage', sample_type, rep_id, include_merged=include_merged )
        rampage_elements = self.get_elements( 
            'rampage', sample_type, rep_id, include_merged=include_merged )

        assert len(cage_elements) == 0 or len(rampage_elements) == 0, \
            "Can not use both RAMPAGE and CAGE reads in a single sample"
        if len(cage_elements) > 0: 
            elements = cage_elements
            reads_class = CAGEReads
        elif len(rampage_elements) > 0:
            elements = rampage_elements
            reads_class = RAMPAGEReads
        else: 
            return []

        promoter_reads = []
        for data in elements:
            if data.filename in self.mapped_reads_cache:
                reads = self.mapped_reads_cache[data.filename]
                reads.reload()
            else:
                rev_reads = {'forward':False, 'backward':True, 'auto': None}[
                    data.read_type]
                reads = reads_class(data.filename)
                reads.init(reverse_read_strand=rev_reads, ref_genes=self.ref_genes)
                self.mapped_reads_cache[data.filename] = reads
            promoter_reads.append(reads)
        
        return promoter_reads

    def _load_polya_reads(self, sample_type, rep_id, include_merged):
        if config.VERBOSE: config.log_statement(
            "Initializing poly(A) reads for sample '%s'." % sample_type)
        all_reads = []
        for data in self.get_elements( 
                'polya', sample_type, rep_id, include_merged=include_merged ):
            if data.filename in self.mapped_reads_cache:
                reads = self.mapped_reads_cache[data.filename]
                reads.reload()
            else:
                assert data.stranded == 'true', "polya-site-seq reads must be stranded"
                rev_reads = {'forward':False, 'backward':True, 'auto': None}[
                    data.read_type]
                reads = PolyAReads(data.filename)
                reads.init(reverse_read_strand=rev_reads, ref_genes=self.ref_genes)
                self.mapped_reads_cache[data.filename] = reads
            all_reads.append(reads)
        
        return all_reads
    
    def get_reads(self, sample_type=None, rep_id=None, 
                  verify_args=True, include_merged=False):
        rnaseq_reads = self._load_rnaseq_reads(
            sample_type, rep_id, include_merged)
        promoter_reads = self._load_promoter_reads(
            sample_type, rep_id, include_merged)
        polya_reads = self._load_polya_reads(
            sample_type, rep_id, include_merged)
        if verify_args:
            self.verify_args_are_sufficient( 
                rnaseq_reads, promoter_reads, polya_reads )
        return (None if len(promoter_reads)==0 else MergedReads(promoter_reads),
                None if len(rnaseq_reads) == 0 else MergedReads(rnaseq_reads), 
                None if len(polya_reads) == 0 else MergedReads(polya_reads) )
    
    def get_sample_types(self):
        query = "SELECT DISTINCT sample_type FROM data"
        with self.conn:
            return [x[0] for x in self.conn.execute(query).fetchall()]

    def get_rep_ids(self, sample_type):
        query = "SELECT DISTINCT rep_id FROM data \
                 WHERE sample_type = '{sample}' AND rep_id != '*'"
        with self.conn:
            return [ x[0] for x in 
                     self.conn.execute(query.format(sample=sample_type)).fetchall()]

def load_ref_elements_to_include(args):
    if None == args.reference and args.use_reference_genes:
        raise ValueError, "--reference must be set if --use-reference-genes is set"
    if None == args.reference and (
            args.use_reference_junctions or args.only_use_reference_junctions):
        raise ValueError, "--reference must be set if --only_use-reference-junctions is set"
    if None == args.reference and args.use_reference_tss_exons:
        raise ValueError, "--reference must be set if --use-reference-tss is set"
    if None == args.reference and args.use_reference_tes_exons:
        raise ValueError, "--reference must be set if --use-reference-tes is set"
    if None == args.reference and args.use_reference_promoters:
        raise ValueError, "--reference must be set if --use-reference-promoters is set"
    if None == args.reference and args.use_reference_polyas:
        raise ValueError, "--reference must be set if --use-reference-polyas is set"
    RefElementsToInclude = namedtuple(
        'RefElementsToInclude', 
        ['genes', 'junctions', 'TSS', 'TES', 'promoters', 'polya_sites'])
    return RefElementsToInclude( args.use_reference_genes, 
                                 ( args.use_reference_junctions 
                                   or args.only_use_reference_junctions),
                                 args.use_reference_tss_exons, 
                                 args.use_reference_tes_exons,
                                 args.use_reference_promoters,
                                 args.use_reference_polyas )

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='Build transcripts and quantify expression levels from RNAseq, CAGE, and poly(A) assays.')

    parser.add_argument( '--control', type=file, 
        help='GRIT control file. Allows better control over the types of input files.')

    parser.add_argument( '--GTF', "-G", type=file, default=None, 
        help='Estimate transcript abundance using the supplied gtf file (do not discover elements or build transcripts).')
    parser.add_argument( '--elements', type=file, default=None, 
        help='Do not find elements - use the provided elements to build transcript models.')
    
    parser.add_argument( '--rnaseq-reads', type=argparse.FileType('rb'), 
        help='BAM file containing mapped RNAseq reads.')
    parser.add_argument( '--rnaseq-read-type', 
                         choices=["forward", "backward", "auto"],
                         default='auto',
        help="If 'forward' then the first RNAseq read in a pair that maps to the genome without being reverse complemented is assumed to be on the correct strand. default: auto")
    parser.add_argument( '--num-mapped-rnaseq-reads', type=int,
        help="The total number of mapped rnaseq reads ( needed to calculate the FPKM ). This only needs to be set if it isn't found by a call to samtools idxstats." )
    
    parser.add_argument( '--cage-reads', type=argparse.FileType('rb'),
        help='BAM file containing mapped cage reads.')
    parser.add_argument( '--cage-read-type', 
                         choices=["forward", "backward", "auto"],
                         default='auto',
        help="If 'forward' then the reads that maps to the genome without being reverse complemented are assumed to be on the '+'. default: auto")

    parser.add_argument( '--rampage-reads', type=argparse.FileType('rb'),
        help='BAM file containing mapped rampage reads.')
    parser.add_argument( '--rampage-read-type', 
                         choices=["forward", "backward", "auto"],
                         default='auto',
        help="If 'forward' then the first read in a pair that maps to the genome without being reverse complemented are assumed to be on the '+' strand. default: auto")
    
    parser.add_argument( '--polya-reads', type=argparse.FileType('rb'), 
        help='BAM file containing mapped polya reads.')
    parser.add_argument( '--polya-read-type', choices=["forward", "backward", "auto"],
                         default='auto',
        help="If 'forward' then the reads that maps to the genome without being reverse complemented are assumed to be on the '+'. default: auto")

    parser.add_argument( '--fasta', type=file,
        help='Fasta file containing the genome sequence - if provided the ORF finder is automatically run.')
    
    parser.add_argument( '--reference', 
                         help='Reference GTF - this option alone only uses the reference for read type inference and gene renaming. No reference elements are used unless explicitly specified (e.g. --use-reference-junctions)',
                         type=file)
    parser.add_argument( '--use-reference-genes', 
                         default=False, action='store_true',
        help='Use genes boundaries from the reference annotation.')
    parser.add_argument( '--use-reference-junctions', 
                         default=False, action='store_true',
        help='Include junctions from the reference annotation.')
    parser.add_argument( '--only-use-reference-junctions', 
                         default=False, action='store_true',
        help='Dont identify novel junctions.')
    parser.add_argument( '--use-reference-tss-exons', 
                         default=False, action='store_true',
        help='Use TSS\'s taken from the reference annotation.')
    parser.add_argument( '--use-reference-tes-exons', 
                         default=False, action='store_true',
        help='Use TES\'s taken from the reference annotation.')
    parser.add_argument( '--use-reference-promoters', 
                         default=False, action='store_true',
        help='Use promoters\'s inferred from the start of reference transcripts.')
    parser.add_argument( '--use-reference-polyas', 
                         default=False, action='store_true',
        help='Use polya sites inferred from the end of reference transcripts.')

    parser.add_argument( '--build-models-with-retained-introns', 
                         default=False, action='store_true',
        help='Include identified retained introns in transcript models.')
    parser.add_argument( '--max-num-candidate-transcripts', 
                         default=1000000, type=int,
        help='Max number of transcript models to build in any gene locus (DEFAULT 1e6).')
    parser.add_argument( '--max-num-transcripts-to-quantify', 
                         default=1000, type=int,
        help='Maximum number of transcript in which to produce quantifications (DEFAULT 1000) (in loci with more than this, we greddily remove low expression elements until we are below the limit).')
    
    parser.add_argument( '--build-bedgraphs', 
                         default=False, action='store_true',
        help='Build read coverage bedgraphs.')
    parser.add_argument( '--only-build-elements', 
                         default=False, action='store_true',
        help='Only build transcript elements - do not build transcript models.')
    parser.add_argument( '--only-build-candidate-transcripts', 
                         default=False, action='store_true',
        help='Do not estimate transcript frequencies - just build transcript models.')
    parser.add_argument( '--dont-estimate-confidence-bounds',
                         default=False,
                         action="store_true",
        help='If set, do not estimate confidence bounds.')
    parser.add_argument( '--dont-estimate-upper-confidence-bounds',
                         default=False,
                         action="store_true",
        help='If set, do not estimate upper confidence bounds.')

    parser.add_argument( '--output-dir', '-o', default="discovered",
        help='Write all output files to this directory. (default: discovered)')
    parser.add_argument( '--continue-run', default=False, action='store_true',
        help='Continue a previously started run in --output-dir')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    parser.add_argument( '--debug-verbose', default=False, action='store_true',
        help='Whether or not to print debugging information.')
    parser.add_argument('--write-debug-data',default=False,action='store_true',
        help='Whether or not to print out gff files containing intermediate exon assembly data.')

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
    if None == args.control and None == args.rnaseq_reads and \
            False == args.only_build_candidate_transcripts:
        raise ValueError, "--control or --rnaseq-reads must be set"

    config.BUILD_MODELS_WITH_RETAINED_INTRONS = \
        args.build_models_with_retained_introns

    config.MAX_NUM_CANDIDATE_TRANSCRIPTS = \
        args.max_num_candidate_transcripts    
 
    config.MAX_NUM_TRANSCRIPTS_TO_QUANTIFY = \
        args.max_num_transcripts_to_quantify   
   
    config.VERBOSE = args.verbose
    
    config.DEBUG_VERBOSE = args.debug_verbose
    
    config.WRITE_DEBUG_DATA = args.write_debug_data
    
    config.NTHREADS = args.threads
    
    config.TOTAL_MAPPED_READS = 1e6
    
    config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS = \
        args.only_build_candidate_transcripts
    assert not (config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS 
                and args.only_build_elements)

    config.ONLY_USE_REFERENCE_JUNCTIONS = \
        args.only_use_reference_junctions

    if ( args.dont_estimate_confidence_bounds
         and args.dont_estimate_upper_confidence_bounds ):
        raise ValueError, "It doesn't make sense to set --dont-estimate-confidence-bounds and --dont-estimate-upper-confidence-bounds"
    
    if (args.only_build_candidate_transcripts and 
            args.dont_estimate_upper_confidence_bounds):
        raise ValueError, "--only-build-candidate-transcripts implies --dont-estimate-confidence-upper-bounds (don't set both)"        

    if (args.only_build_candidate_transcripts and 
            args.dont_estimate_confidence_bounds):
        raise ValueError, "--only-build-candidate-transcripts implies --dont-estimate-confidence-upper-bounds (don't set both)"        

    config.ESTIMATE_UPPER_CONFIDENCE_BOUNDS = True
    config.ESTIMATE_LOWER_CONFIDENCE_BOUNDS = True
    
    if args.dont_estimate_confidence_bounds:
        config.ESTIMATE_UPPER_CONFIDENCE_BOUNDS = False
        config.ESTIMATE_LOWER_CONFIDENCE_BOUNDS = False

    if args.dont_estimate_upper_confidence_bounds:
        config.ESTIMATE_UPPER_CONFIDENCE_BOUNDS = False
    
    config.FIX_CHRM_NAMES_FOR_UCSC = args.ucsc

    args.output_dir = os.path.abspath(args.output_dir)
    config.tmp_dir = os.path.join(args.output_dir, "./.tmp_files/")
    try: 
        os.mkdir(args.output_dir)
        os.mkdir(config.tmp_dir)
    except OSError:
        if not args.continue_run:
            raise OSError, "Directory '%s' already exists: either remove it, change the --output-dir option, or use the --continue-run option" % args.output_dir
    os.chdir(args.output_dir)
    
    log_ofstream = open( "log.txt", "w" )
    log_statement = Logger(
        nthreads=config.NTHREADS,
        use_ncurses=(not args.batch_mode), 
        log_ofstream=log_ofstream)
    config.log_statement = log_statement
    
    if args.region != None:
        args.region = clean_chr_name(args.region)

    args.ref_elements_to_include = load_ref_elements_to_include(args)
    return args

class discover_elements(dict):
    """A lazy elements dictionary. 
    
    It initializes with nones and then builds elements as they're requested.
    This allows elements and transcript building to proceed in sync.
    """
    def _build_elements(self, sample_type):
        if config.VERBOSE: config.log_statement("Initializing read objects.")
        promoter_reads, rnaseq_reads, polya_reads = \
            self.sample_data.get_reads(sample_type, include_merged=True)
        elements_fname = "%s.elements.bed" % sample_type
        try: 
            ofp = open(elements_fname)
        except IOError:
            grit.find_elements.find_elements(
                promoter_reads, rnaseq_reads, polya_reads,
                elements_fname, self.sample_data.ref_genes, 
                self.args.ref_elements_to_include, 
                region_to_use=self.args.region)
            ofp = open(elements_fname)
        else:
            msg = "WARNING: '%s' already exists - using existing file."
            config.log_statement(msg % elements_fname, log=True)
        return  (ofp, None)

    def __init__(self, sample_data, args):
        """Discover elements for all samples

        """
        self.sample_data = sample_data
        self.args = args
        for sample_type in sample_data.get_sample_types():
            self[sample_type] = None
        return
    
    def iteritems(self):
        for sample_type in self.keys():
            if self[sample_type] == None:
                self[sample_type] = self._build_elements(sample_type)
            yield sample_type, self[sample_type]
        return

def main():
    args = parse_arguments()

    # load the samples into database, and the reference genes if necessary
    sample_data = Samples(args)

    # if the reference genes weren't loaded while parsing the data, and we
    # need the reference elements, then load the reference genes now
    if args.reference != None and sample_data.ref_genes == None:
        if config.VERBOSE: 
            config.log_statement("Loading annotation file.", log=True)
        sample_data.ref_genes = load_gtf(args.reference)
    
    # find elements if necessary, load the gtf if we are running in 
    # quantification mode
    if args.GTF != None:
        assert not config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS
        assert not args.only_build_elements
        assert args.elements == None, "It doesn't make sense to set --GTF and --elements - did you mean to set --reference instead of --GTF?"
        sample_types = sample_data.get_sample_types()
        if len(sample_types) != 1: 
            raise ValueError, "Can not provide --elements and data from multiple sample types."
        elements = {sample_types[0]: (None, args.GTF)}
    elif args.elements !=  None:
        assert not args.only_build_elements
        sample_types = sample_data.get_sample_types()
        if len(sample_types) > 1: 
            raise ValueError, "Can not provide --elements and data from multiple sample types."
        elif len(sample_types) == 1:
            sample_type = sample_types[0]
        elif len(sample_types) == 0:
            sample_type = os.path.basename(args.elements.name)
        elements = {sample_type: (args.elements, None)}
    else:
        elements = discover_elements(sample_data, args)

    # if we are only building elements, then we are done
    if args.only_build_elements:
        return 

    # build transcripts for each sample
    for sample_type, (elements_fp, gtf_fp) in elements.iteritems():
        gtf_fname = sample_type + ".gtf"
        tracking_fname = sample_type + ".transcript_tracking"
        assert elements_fp == None or gtf_fp == None
        # try to load an already built gtf
        if args.continue_run and gtf_fp == None:
            try: gtf_fp = open(gtf_fname)
            except IOError: pass
            else:
                msg = "WARNING: '%s' already exists - using existing file."
                config.log_statement(msg % gtf_fname, log=True)
        if gtf_fp != None:
            config.log_statement( "Loading %s" % gtf_fp.name )
            genes_fnames = []
            genes = load_gtf(gtf_fp)
            elements = extract_elements_from_genes(genes)
            for gene in genes:
                gene_fname = config.get_gene_tmp_fname(
                    gene.id, sample_type)
                genes_fnames.append(
                    (gene.id, 
                     len(gene.transcripts), 
                     gene.write_to_file(gene_fname)))
            config.log_statement("Finished Loading %s" % gtf_fp.name)
        else:
            elements = load_elements(elements_fp)
            genes_fnames = grit.build_transcripts.build_transcripts(
                elements_fp, gtf_fname, tracking_fname, 
                args.fasta, sample_data.ref_genes)
        if config.ONLY_BUILD_CANDIDATE_TRANSCRIPTS: continue
        rep_ids = sample_data.get_rep_ids(sample_type)
        
        if len(rep_ids) == 0: rep_ids = [None,]
        for rep_id in rep_ids:
            config.log_statement("Loading reads for %s-%s" % (
                    sample_type, rep_id))
            (promoter_reads, rnaseq_reads, polya_reads) = sample_data.get_reads(
                sample_type, rep_id, 
                verify_args=False, include_merged=False)
            
            # find what to prefix the output files with
            if rep_id == None:
                exp_ofname = "%s.expression_tracking" % sample_type
            else:
                exp_ofname = "%s.%s.expression_tracking" % (sample_type, rep_id)
            
            fldist_fname = "%s.%s.fldists" % (sample_type, rep_id)
            try: 
                with open(fldist_fname) as fp:
                    config.log_statement("Loading fldists")
                    fl_dists = pickle.load(fp)
            except IOError:
                config.log_statement(
                    "Estimating fragment length distribution", log=True)
                fl_dists = grit.frag_len.build_fl_dists(elements, rnaseq_reads)
                with open(fldist_fname, "w") as ofp:
                    pickle.dump(fl_dists, ofp)
            
            grit.estimate_transcript_expression.quantify_transcript_expression(
                promoter_reads, rnaseq_reads, polya_reads,
                genes_fnames, fl_dists,
                exp_ofname )
    
if __name__ == '__main__':
    try: main()
    finally: 
        try:config.log_statement.close()
        except: pass
