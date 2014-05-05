VERSION = "1.1.3-dev"

# config options shared across modules

# the maximum number of transcripts to build a design matrix for
MAX_NUM_TRANSCRIPTS = 5000
# the maximum number of candidate transcripts to build a particualr gene locus
MAX_NUM_CANDIDATE_TRANSCRIPTS = 25000

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

MIN_INTRON_SIZE = 50
MIN_EXON_SIZE = 50 #MIN_INTRON_SIZE
MAX_EMPTY_REGION_SIZE = 200 #MIN_INTRON_SIZE
MAX_INTRON_SIZE = int(1e6)
MIN_GENE_LENGTH = 400
# the maximum number of bases to expand gene boundaries from annotated genes
MAX_GENE_EXPANSION = 1000

MIN_EXON_BPKM = 0.1
EXON_EXT_CVG_RATIO_THRESH = 3
#POLYA_MERGE_SIZE = 100

CAGE_PEAK_WIN_SIZE = 15
CAGE_FILTER_ALPHA = 0.9
MIN_NUM_CAGE_TAGS = 5
MAX_CAGE_FRAC = 0.10
NUM_TSS_BASES_TO_SKIP = 200

MIN_NUM_POLYA_TAGS = 2
NUM_TES_BASES_TO_SKIP = 300

def log_statement(*args):
    print args[0]
