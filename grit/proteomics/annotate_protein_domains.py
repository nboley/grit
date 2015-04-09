"""
Copyright (c) 2011-2015 Nathan Boley, Marcus Stoiber

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
import subprocess

from tempfile import NamedTemporaryFile

from pysam import Fastafile

from proteomics import find_coding_sequence, format_into_80_char_lines
from ..files.gtf import load_gtf, Transcript, partition_coding_and_utr_segments, \
    create_gff_line, GenomicInterval

CDDB="/home/nboley/conserved_domain_analysis/Cdd_ALL/Smart"
CDDBmapping_fname="/home/nboley/conserved_domain_analysis/Cdd_ALL/cddid_all.tbl"

def find_domains(protein_seq, min_eval=1e-5):
    # write the protein to a file
    ofp = NamedTemporaryFile(suffix='.fa')
    ofp.write( ">1\n" )
    for line in format_into_80_char_lines( protein_seq ):
        ofp.write( line + "\n" )
    ofp.flush()
    
    cmd = "rpsblast -query %s -db %s -evalue %e -outfmt 7" % ( 
        ofp.name, CDDB, min_eval )
    res = subprocess.check_output( cmd, shell=True )
    ofp.close()
    rv = []
    for line in res.splitlines():
        if line.startswith("#"): continue
        data = line.split()
        start, stop = int(data[6]), int(data[7])
        code = int(data[1].split("|")[-1])
        rv.append( (code, start, stop) )
    
    return rv

def load_id_mapping():
    mapping = {}
    with open( CDDBmapping_fname ) as fp:
        for line in fp:
            data = line.strip().split("\t")
            mapping[int(data[0])] = data[2]
    
    return mapping

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(
        description = 'Annotate functional domains in proteins.' )
    parser.add_argument(
        'gtf', type=file,
        help='GTF file to search for ORFs.' )
    parser.add_argument(
        'fasta',
        help='Fasta file with reference sequence.' )
        
    parser.add_argument(
        '--output-filename', '-o',
        help='Output file. (default: stdout)')
    
    parser.add_argument( 
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    args = parser.parse_args()
        
    # create default if no prefix provided or if same as gtf filename
    if args.output_filename == None:
        ofp = sys.stdout
    else:
        ofp = open( args.output_filename, 'w' )
    
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.gtf, args.fasta, ofp

def convert_to_genome_coords( t, protein_seq, domains ):
    for d_id, start, stop in domains:
        if t.strand == '-': 
            utr_len = sum( b-a+1 for a,b in t.tp_utr_exons )
            start, stop = len(protein_seq)-stop, len(protein_seq)-start
        else:
            utr_len = sum( b-a+1 for a,b in t.fp_utr_exons )
        
        start, stop = t.genome_pos(utr_len+3*start),t.genome_pos(utr_len+3*stop)
        us, mid, ds = partition_coding_and_utr_segments( t.exons, start, stop )
        for start, stop in mid:
            yield d_id, ( t.chrm, t.strand, start, stop )
    
    return

def main():
    gtf_fp, fasta_fn, ofp = parse_arguments()
    mapping = load_id_mapping()
    genes = load_gtf( gtf_fp.name )
    fa = Fastafile( fasta_fn )
    for gene in genes:
        for i, t in enumerate(gene.transcripts):
            if not t.is_protein_coding: continue
            protein_seq = find_coding_sequence( t, fa )
            domains = find_domains(protein_seq, min_eval=1e-5)
            for d_id, region in convert_to_genome_coords( 
                    t, protein_seq, domains ):
                name = "%s.%i.%s" % (gene.meta_data['gene_name'], i, mapping[d_id])
                print create_gff_line( GenomicInterval(*region), name )
                                       
                                       
    
    return
    

if __name__ == '__main__':
    main()
