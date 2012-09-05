import sys, os

sys.path.append(os.path.join(os.path.dirname(__file__), "../", 'file_types'))
from wiggle import Wiggle, load_wiggle, guess_strand_from_fname



def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(\
        description='Find exons from wiggle and junctions files.')

    parser.add_argument( 'junctions', type=file, \
        help='GTF format file of junctions(introns).')
    parser.add_argument( 'chrm_sizes_fname', type=file, \
        help='File with chromosome names and sizes.')

    parser.add_argument( 'wigs', type=file, nargs="+", \
        help='wig files over which to search for exons.')
    
    parser.add_argument( '--cage-wigs', type=file, nargs='+', \
        help='wig files with cage reads, to identify tss exons.')
    parser.add_argument( '--polya-candidate-sites', type=file, nargs='*', \
        help='files with allowed polya sites.')
    
    parser.add_argument( '--out-file-prefix', '-o', default="discovered_exons",\
        help='Output file name. (default: discovered_exons)')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
        help='Whether or not to print status information.')
    parser.add_argument( '--threads', '-t', default=1, type=int,
        help='The number of threads to use.')
        
    args = parser.parse_args()

    global num_threads
    num_threads = args.threads

    ofps_prefixes = [ "single_exon_genes", 
                      "tss_exons", "internal_exons", "tes_exons", 
                      "all_exons"]
    
    fps = []
    for field_name in ofps_prefixes:
        fps.append(open("%s.%s.gff" % (args.out_file_prefix, field_name), "w"))
    ofps = OrderedDict( zip( ofps_prefixes, fps ) )
    
    # set flag args
    global VERBOSE
    VERBOSE = args.verbose
    
    rd1_plus_wigs = [ fp for fp in args.wigs 
                      if fp.name.endswith("rd1.plus.bedGraph") ]
    rd1_minus_wigs = [ fp for fp in args.wigs 
                      if fp.name.endswith("rd1.minus.bedGraph") ]
    rd2_plus_wigs = [ fp for fp in args.wigs 
                      if fp.name.endswith("rd2.plus.bedGraph") ]
    rd2_minus_wigs = [ fp for fp in args.wigs 
                      if fp.name.endswith("rd2.minus.bedGraph") ]
    
    grpd_wigs = [ rd1_plus_wigs, rd1_minus_wigs, rd2_plus_wigs, rd2_minus_wigs ]
    
    return grpd_wigs, args.junctions, args.chrm_sizes_fname, \
        args.cage_wigs, args.polya_candidate_sites, ofps


def main():
    wigs, jns_fp, chrm_sizes_fp, cage_wigs, polya_candidate_sites_fps, out_fps \
        = parse_arguments()
    
    # set up all of the file processing calls
    p = multiprocessing.Pool( processes=num_threads )
    
    if VERBOSE: print >> sys.stderr,  'Loading read pair 1 wiggles'    
    args = [ chrm_sizes_fp.name, [wigs[0][0].name, wigs[1][0].name], ['+','-'] ]
    rd1_cov_proc = p.apply_async( load_wiggle, args )
    
    if VERBOSE: print >> sys.stderr,  'Loading merged read pair wiggles'    
    fnames =[wigs[0][0].name, wigs[1][0].name, wigs[2][0].name, wigs[3][0].name]
    rd_cov_proc = p.apply_async( load_wiggle,
        [ chrm_sizes_fp.name, fnames, ['+', '-', '+', '-'] ] )
    
    if VERBOSE: print >> sys.stderr,  'Loading read pair 2 wiggles'    
    rd2_cov_proc = p.apply_async( load_wiggle,
        [ chrm_sizes_fp.name, [wigs[2][0].name, wigs[3][0].name], ['+','-'] ] )

    if VERBOSE: print >> sys.stderr,  'Loading CAGE.'
    cage_cov_proc = p.apply_async( 
        load_wiggle, [ chrm_sizes_fp.name, [x.name for x in cage_wigs] ] )
    
    if VERBOSE: print >> sys.stderr,  'Loading junctions.'
    jns_proc = p.apply_async( 
        parse_junctions_file_dont_freeze, [ jns_fp.name, ] )
    
    if VERBOSE: print >> sys.stderr,  'Loading candidate polyA sites'
    polya_sites_proc = p.apply_async( 
        find_polya_sites,  [[x.name for x in polya_candidate_sites_fps],] )
    
    # now, all of the async calls have been made so we collect the results
    
    polya_sites = polya_sites_proc.get()
    for fp in polya_candidate_sites_fps: fp.close()
    if VERBOSE: print >> sys.stderr, 'Finished loading candidate polyA sites'

    rd1_cov = rd1_cov_proc.get()
    if VERBOSE: print >> sys.stderr, 'Finished loading read pair 1 wiggles'    
    
    rd2_cov = rd2_cov_proc.get()
    if VERBOSE: print >> sys.stderr, 'Finished loading read pair 2 wiggles'    

    read_cov = rd_cov_proc.get()
    if VERBOSE: print >> sys.stderr, 'Finished loading merged read pair wiggles'
        
    # open the cage data
    cage_cov = cage_cov_proc.get()
    cage_sum = sum( cage_cov.apply( lambda a: a.sum() ).values() )
    for cage_fp in cage_wigs: cage_fp.close()
    if VERBOSE: print >> sys.stderr, 'Finished loading CAGE data'
    
    jns = jns_proc.get() 
    if VERBOSE: print >> sys.stderr,  'Finished loading junctions.'
    
    if read_cov_sum == 0:
        [ fp.close() for fp in out_fps ]
        return


    all_regions_iters = [ [], [], [], [], [] ]

    keys = sorted( set( jns ) )
    for chrm, strand in keys:        
        if VERBOSE: print >> sys.stderr, \
                'Processing chromosome %s strand %s.' % ( chrm, strand )
        
        if VERBOSE: print >> sys.stderr, \
                'Loading read coverage data (%s,%s)' % ( chrm, strand )
        
        read_cov_obj = ReadCoverageData( \
            read_cov.zero_intervals[(chrm, strand)], 
            read_cov[(chrm, strand)] )
        
        introns = jns[(chrm, strand)]
        scores = jns._scores[(chrm, strand)]
        jns_and_scores = zip( introns, scores )
        
        disc_grpd_exons = find_exons_in_contig( \
           strand, 
           read_cov_obj, 
           jns_and_scores,
           cage_cov[ (chrm, strand) ], 
           polya_sites[ (chrm, strand) ] )


        if VERBOSE: print >> sys.stderr, \
                'Building output (%s,%s)' % ( chrm, strand )
        
        for container, exons in zip( all_regions_iters, disc_grpd_exons ):
            regions_iter = ( GenomicInterval( chrm, strand, start, stop) \
                                 for start, stop in exons                \
                                 if stop - start < MAX_EXON_SIZE )
            container.extend( regions_iter )

    # process each chrm, strand combination separately
    for track_name, out_fp in out_fps.iteritems():
        out_fp.write( "track name=%s\n" % track_name )
    
    if VERBOSE: print >> sys.stderr,  'Writing output.'
    for regions_iter, ofp in zip( all_regions_iters, out_fps.itervalues() ):
        ofp.write( "\n".join(iter_gff_lines( sorted(regions_iter) )) + "\n")
        
    return
        
if __name__ == "__main__":
    main()
