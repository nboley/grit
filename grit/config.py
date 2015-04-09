"""
Copyright (c) 2011-2015 Nathan Boley

This file is part of GRIT.

GRIT is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GRIT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GRIT.  If not, see <http://www.gnu.org/licenses/>.
"""

import os, sys

# config options shared across modules

# the maximum number of transcripts to produce quantification estimates for
MAX_NUM_TRANSCRIPTS_TO_QUANTIFY = 1000
# thxe maximum number of candidate transcripts to build in a particular gene locus
MAX_NUM_CANDIDATE_TRANSCRIPTS = 50000

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

ELEMENT_FILTER_FRAC = 0.05

# the minimum number of flanking bases we need to identify a junction
# from a gapped read alignmento
MIN_INTRON_FLANKING_SIZE = 12 
MIN_ENTROPY = 0
MIN_INTRON_SIZE = 40
MAX_EMPTY_REGION_SIZE = 200 #MIN_INTRON_SIZE
MAX_INTRON_SIZE = int(1e6)
MIN_GENE_LENGTH = 400
# the maximum number of bases to expand gene boundaries from annotated genes
MAX_GENE_EXPANSION = 1000

NOISE_JN_FILTER_FRAC = 0.01
MAX_JN_OFFSET_FILTER = 15

MAX_FRAGMENT_LENGTH = 1000
MIN_FRAGMENT_LENGTH = 75

ONLY_USE_REFERENCE_JUNCTIONS = False
MIN_EXON_FPKM = 1
MAX_EXPRESSION_RATIO = 10
BUILD_MODELS_WITH_RETAINED_INTRONS = False

MIN_NUM_CAGE_TAGS = 5
NUM_TSS_BASES_TO_SKIP = 200

MIN_NUM_POLYA_TAGS = 2
NUM_TES_BASES_TO_SKIP = 300

TES_EXON_MERGE_DISTANCE = 0
TSS_EXON_MERGE_DISTANCE = 0
MAX_DISTAL_SIZE_FOR_MATCH_OFFSET = 100

tmp_dir = None

def get_gene_tmp_fname(gene_id, sample_type=None, rep_id=None):
    rv = os.path.join(tmp_dir, "%s" % gene_id )
    if sample_type != None: rv += ".%s" % sample_type
    if rep_id != None: rv += ".%s" % rep_id
    return rv + ".gene"


def get_fmat_tmp_fname(gene_id, sample_type=None, rep_id=None):
    rv = os.path.join(tmp_dir, "%s" % gene_id )
    if sample_type != None: rv += ".%s" % sample_type
    if rep_id != None: rv += ".%s" % rep_id
    return rv + ".fmat"

def log_statement(*args, **kwargs):
    print args[0]
