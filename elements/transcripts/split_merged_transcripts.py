import sys, os

sys.path.append( os.path.join( os.path.dirname(__file__), 
                               '..', '..', 'file_types', 'fast_gtf_parser' ) )
from gtf import load_gtf, Transcript

def parse_sources( sources_fp ):
    sources = {}
    for line in sources_fp:
        gene_id, trans_id, trans_sources = line.split()
        trans_sources = trans_sources.split(',')
        
        sources[(gene_id, trans_id)] = trans_sources
    
    return sources

def find_sources( sources, gene_id, trans_id ):
    """Find the op file pointers corresponding to (gene_id, trans_id).

    This would be trivial, except that if we've passed the transcripts
    through ORF finder, then, if a transcript had multiple ORFs, we may
    have multiple similar ids. In such cases, since orf finder uses
    trans_id_ORF%i, we can simply parse the name.
    """
    if ( gene_id, trans_id ) in sources:
        return sources[ (gene_id, trans_id) ]
    
    # if this is not, assume that it is a CDS modification.
    base_trans_id = trans_id.split("_")[:-1].join( "_" )
    return sources[ (gene_id, base_trans_id) ]

def get_ofp_from_source_fname(sources_fp_mapping, source_fname, opdir, suffix):
    # if we've already observed this, then just return it
    if source_fname in sources_fp_mapping:
        return sources_fp_mapping[ source_fname ]
    
    ofname = os.path.join( opdir, source_fname + suffix )
    ofp = open( ofname, "w" )
    sources_fp_mapping[ source_fname ] = ofp
    
    return ofp

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(
        description = 'Recover sample by sample gtfs from a merged gtf file ' \
            + 'and a sources file. Deal with transcripts that code for ' \
            + 'multiple proteins intelligently.' )
    parser.add_argument(
        'merged_transcripts', type=file,
        help='GTF file containing meged transcripts.' )
    parser.add_argument(
        'sources_file', type=file,
        help='sources file from merge_transcripts.py.')

    parser.add_argument(
        '--output-dir', type=str, default='.',
        help='The directory in which to place the output files (default . )')
    parser.add_argument(
        '--suffix', type=str, default="",
        help='Suffix to append to the output files.')
    parser.add_argument( 
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    
    args = parser.parse_args()
    
    # set flag args
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.merged_transcripts, args.sources_file, \
        args.output_dir, args.suffix

def main():
    trans_fp, sources_fp, opdir, suffix = parse_arguments()

    if VERBOSE: print >> sys.stderr, 'parsing sources'
    sources = parse_sources( sources_fp )
    sources_fp.close()

    if VERBOSE: print >> sys.stderr, 'parsing transcripts gtf'
    genes = load_gtf( trans_fp.name )
    trans_fp.close()
    
    sources_fp_mapping = {}
    
    for gene in genes:
        for trans in gene.transcripts:
            source_fnames = find_sources( sources, gene.id, trans.id )
            for source_fn in source_fnames:
                ofp = get_ofp_from_source_fname( 
                    sources_fp_mapping, source_fn, opdir, suffix )
                ofp.write( trans.build_gtf_lines( gene.id, {} ) + "\n" )

    for fp in sources_fp_mapping.values():
        fp.close()
    
    return

if __name__ == '__main__':
    main()
