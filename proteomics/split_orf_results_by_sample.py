import sys, os

from collections import defaultdict

DO_PROFILE = False

def write_orfs_for_all_sources( disc_orfs, sources, out_prefix, out_suffix ):
    sources_fps = {}
    for key, lines in disc_orfs.iteritems():
        for source in sources[key]:
            if source in sources_fps:
                sources_fps[source].write( lines )
            else:
                sources_fps[source] = open( out_prefix + source + out_suffix, 'w' )
                sources_fps[source].write( lines )
    
    for source_fp in sources_fps.itervalues():
        source_fp.close()
    
    return

def parse_sources( sources_fp ):
    sources = {}
    for line in sources_fp:
        gene_id, trans_id, trans_sources = line.split()
        trans_sources = trans_sources.split(',')
        
        sources[(gene_id, trans_id)] = trans_sources
    
    return sources

def parse_discovered_orfs( disc_orfs_fp ):
    disc_orfs = defaultdict( str )
    for line in disc_orfs_fp:
        fields = line.split()
        if len( fields ) < 12: continue
        if fields[2] not in ( 'CDS', '5UTR', '3UTR' ): continue
        
        gene_id = fields[9][1:-2]
        trans_id = '_'.join( fields[11][1:-2].split('_')[:-1]  )
        disc_orfs[(gene_id, trans_id)] += line
    
    return disc_orfs

def build_objects( disc_orfs_fp, sources_fp ):    
    if VERBOSE: print >> sys.stderr, 'parsing discovered orfs'
    disc_orfs = parse_discovered_orfs( disc_orfs_fp )
    disc_orfs_fp.close()
    
    if VERBOSE: print >> sys.stderr, 'parsing sources'
    sources = parse_sources( sources_fp )
    sources_fp.close()
    
    return disc_orfs, sources

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(
        description = 'Find discovered alternatively splice ORFs that do not '
        + 'shift the frame of the known ORF.' )
    parser.add_argument(
        'orf_finder_output', type=file,
        help='GTF file with CDS and UTR regions.' )
    parser.add_argument(
        'source_file', type=file,
        help='sources_file from merge_transcripts.py.')
    
    parser.add_argument(
        '--out-prefix',
        help='Prefix or path of output files.')
    parser.add_argument(
        '--out-suffix', default='.annotated.gtf',
        help='Suffix of output file. (default: %(default)s)')
    parser.add_argument( 
        '--verbose', '-v', default=False, action='store_true',
        help='Whether or not to print status information.')
    args = parser.parse_args()
    
    # set flag args
    global VERBOSE
    VERBOSE = args.verbose
    
    return args.orf_finder_output, args.source_file, args.out_prefix, args.out_suffix

def main():
    disc_orfs_fp, sources_fp, out_prefix, out_suffix = parse_arguments()
    disc_orfs, sources = build_objects( disc_orfs_fp, sources_fp )
    
    write_orfs_for_all_sources( disc_orfs, sources, out_prefix, out_suffix )
    
    return

if __name__ == '__main__':
    if DO_PROFILE:
        import cProfile
        cProfile.run( 'main()' )
    else:
        main()
