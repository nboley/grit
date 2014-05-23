import os, sys

VERSION = "1.1.3-dev"

# config options shared across modules

# the maximum number of transcripts to produce quantification estimates for
MAX_NUM_TRANSCRIPTS_TO_QUANTIFY = 1000
# the maximum number of candidate transcripts to build in a particular gene locus
MAX_NUM_CANDIDATE_TRANSCRIPTS = 1000000

CB_SIG_LEVEL = 0.025

# log statement is set in the main init, and is a global
# function which facilitates smart, ncurses based logging
log_statement = None
VERBOSE = None
DEBUG_VERBOSE = None
WRITE_DEBUG_DATA = None

FIX_CHRM_NAMES_FOR_UCSC = None
ONLY_BUILD_CANDIDATE_TRANSCRIPTS = None

NTHREADS = None
TOTAL_MAPPED_READS = None

ESTIMATE_UPPER_CONFIDENCE_BOUNDS = True
ESTIMATE_LOWER_CONFIDENCE_BOUNDS = True

# the minimum number of flanking bases we need to identify a junction
# from a gapped read alignmento
MIN_INTRON_FLANKING_SIZE = 12 
MIN_ENTROPY = 0.0
MIN_INTRON_SIZE = 40
MIN_EXON_SIZE = 40 #MIN_INTRON_SIZE
MAX_EMPTY_REGION_SIZE = 200 #MIN_INTRON_SIZE
MAX_INTRON_SIZE = int(1e6)
MIN_GENE_LENGTH = 400
# the maximum number of bases to expand gene boundaries from annotated genes
MAX_GENE_EXPANSION = 1000
NOISE_JN_FILTER_FRAC = 0.05

ONLY_USE_REFERENCE_JUNCTIONS = False
MIN_EXON_BPKM = 0.1
EXON_EXT_CVG_RATIO_THRESH = 5
BUILD_MODELS_WITH_RETAINED_INTRONS = False
#POLYA_MERGE_SIZE = 100

CAGE_PEAK_WIN_SIZE = 15
CAGE_FILTER_ALPHA = 0.9
MIN_NUM_CAGE_TAGS = 5
MAX_CAGE_FRAC = 0.10
NUM_TSS_BASES_TO_SKIP = 200

MIN_NUM_POLYA_TAGS = 2
NUM_TES_BASES_TO_SKIP = 300

TES_EXON_MERGE_DISTANCE = 500
TSS_EXON_MERGE_DISTANCE = 500

tmp_dir = None

def get_gene_tmp_fname(gene_id, sample_type=None, rep_id=None):
    return os.path.join(tmp_dir, "%s.gene" % gene_id)

def get_fmat_tmp_fname(gene_id, sample_type=None, rep_id=None):
    return os.path.join(tmp_dir, "%s.fmat" % gene_id)

def log_statement(*args):
    print args[0]
