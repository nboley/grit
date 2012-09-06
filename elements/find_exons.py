import sys, os
import numpy
import multiprocessing

from collections import OrderedDict, defaultdict
from itertools import izip

sys.path.append(os.path.join(os.path.dirname(__file__), "../", 'file_types'))
from wiggle import Wiggle, load_wiggle, guess_strand_from_fname
from junctions_file import parse_jn_gff, Junctions

def build_empty_array():
    return numpy.array(())

def find_polya_sites( polya_sites_fnames ):
    locs = defaultdict( list )
    for fname in polya_sites_fnames:
        strand = guess_strand_from_fname( fname )
        with open( fname ) as fp:
            for line in fp:
                if line.startswith( "track" ): continue
                data = line.split()
                chrm, start, stop, value = \
                    data[0], int(data[1]), int(data[2]), float(data[3])
                assert start == stop
                assert value == 1
                locs[(chrm, strand)].append( start )
    
    # convert to a dict of sorted numpy arrays
    numpy_locs = defaultdict( build_empty_array )

    for (chrm, strand), polya_sites in locs.iteritems():
        # make sure they're unique
        assert len( polya_sites ) == len( set( polya_sites ) )

        polya_sites.sort()
        if chrm.startswith( 'chr' ):
            chrm = chrm[3:]
        
        numpy_locs[(chrm, strand)] = numpy.array( polya_sites )
    
    return numpy_locs

class Bins( list ):
    def __init__( self, chrm, strand, iter=[] ):
        self.chrm = chrm
        self.strand = strand
        self.extend( iter )
        self._bed_template = "\t".join( ["chr"+chrm, '{start}', '{stop}', '{name}', 
                                         '1000', strand, '{start}', '{stop}', 
                                         '{color}']  ) + "\n"
        
    def writeBedgraph( self, ofp ):
        """
            chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
        """
        ofp.write('track name="bins" visibility=2 itemRgb="On"\n')
        for bin in self:
            length = max( (bin.stop - bin.start)/4, 1)
            colors = bin._find_colors()
            if isinstance( colors, str ):
                op = self._bed_template.format(
                    start=bin.start,stop=bin.stop-1,color=colors, 
                    name="%s_%s"%(bin.left_label, bin.right_label) )
                ofp.write( op )
            else:
                op = self._bed_template.format(
                    start=bin.start,stop=(bin.start + length),color=colors[0], 
                    name=bin.left_label)
                ofp.write( op )
                
                op = self._bed_template.format(
                    start=(bin.stop-1-length),stop=bin.stop-1,color=colors[1],
                    name=bin.right_label)
                ofp.write( op )
        
        return

class Bin( object ):
    def __init__( self, start, stop, left_label, right_label ):
        self.start = start
        self.stop = stop
        self.left_label = left_label
        self.right_label = right_label
        
    def __repr__( self ):
        return "%i-%i\t%s\t%s" % ( self.start, self.stop, self.left_label,
                                   self.right_label )

    _bndry_color_mapping = {
        'START': '0,0,0',
        'STOP': '0,0,0',
        
        'POLYA': '255,255,0', #
        
        'D_JN': '173,255,47', # 
        'R_JN': '0,255,0',
        
        'E_ST': '0,255,255',
        'E_SP': '0,0,255'    
    }
    
    def find_bndry_color( self, bndry ):
        return self._bndry_color_mapping[ bndry ]
    
    def _find_colors( self ):
        if self.left_label == 'E_ST' and self.right_label  == 'E_SP':
            return '0,0,0'
        if self.left_label == 'D_JN' and self.right_label  == 'R_JN':
            return '108,108,108'
        if self.left_label == 'D_JN' and self.right_label  == 'D_JN':
            return '135,206,250'
        if self.left_label == 'R_JN' and self.right_label  == 'R_JN':
            return '135,206,250'
        if self.left_label == 'R_JN' and self.right_label  == 'D_JN':
            return '0,0,255'
        if self.left_label == 'R_JN' and self.right_label  == 'POLYA':
            return '255,0,0'
        if self.left_label == 'POLYA' and self.right_label  == 'POLYA':
            return ' 240,128,128'
        if self.left_label == 'POLYA' and self.right_label  == 'D_JN':
            return '147,112,219'
        if self.left_label == 'POLYA' and self.right_label  == 'R_JN':
            return '159,153,87'
        
        return ( self.find_bndry_color(self.left_label), 
                 self.find_bndry_color(self.right_label) )
        


def find_exons_in_contig( ( chrm, strand ), rnaseq_cov, jns, cage_cov, polya_sites ):
    assert strand == '+'
    
    locs = {}
    for polya in polya_sites:
        locs[ polya ] = "POLYA"
    
    for start, stop in jns:
        locs[start-1] = "D_JN"
        locs[stop+1] = "R_JN"

    """
    for start, stop in find_empty_regions( rnaseq_cov ):
        locs[start-1] = "E_ST"
        locs[stop+1] = "E_SP"
    """
    
    poss = sorted( locs.iteritems() )
    bins = Bins( chrm, strand, [ Bin(0, poss[0][0], "START", poss[0][1]), ])
    for index, ((start, left_label), (stop, right_label)) in \
            enumerate(izip(poss[:-1], poss[1:])):
        mean_cov = rnaseq_cov[start:stop].mean()
        bins.append( Bin(start, stop, left_label, right_label) )
        
    bins.append( Bin( poss[-1][0], len(rnaseq_cov), poss[-1][1], "STOP" ) )
    
    bins.writeBedgraph( sys.stdout )
    
    return

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
    
    if VERBOSE: print >> sys.stderr,  'Loading merged read pair wiggles'    
    fnames =[wigs[0][0].name, wigs[1][0].name, wigs[2][0].name, wigs[3][0].name]
    rd_cov_proc = p.apply_async( load_wiggle,
        [ chrm_sizes_fp.name, fnames, ['+', '-', '+', '-'] ] )
    
    if VERBOSE: print >> sys.stderr,  'Loading CAGE.'
    cage_cov_proc = p.apply_async( 
        load_wiggle, [ chrm_sizes_fp.name, [x.name for x in cage_wigs] ] )
    
    if VERBOSE: print >> sys.stderr,  'Loading junctions.'
    jns_proc = p.apply_async( 
        parse_jn_gff, [ jns_fp.name, ] )
    
    if VERBOSE: print >> sys.stderr,  'Loading candidate polyA sites'
    polya_sites_proc = p.apply_async( 
        find_polya_sites,  [[x.name for x in polya_candidate_sites_fps],] )
    
    # now, all of the async calls have been made so we collect the results
    
    polya_sites = polya_sites_proc.get()
    for fp in polya_candidate_sites_fps: fp.close()
    if VERBOSE: print >> sys.stderr, 'Finished loading candidate polyA sites'

    read_cov = rd_cov_proc.get()
    if VERBOSE: print >> sys.stderr, 'Finished loading merged read pair wiggles'
        
    # open the cage data
    cage_cov = cage_cov_proc.get()
    cage_sum = sum( cage_cov.apply( lambda a: a.sum() ).values() )
    for cage_fp in cage_wigs: cage_fp.close()
    if VERBOSE: print >> sys.stderr, 'Finished loading CAGE data'
    
    jns = Junctions( jns_proc.get() )
    if VERBOSE: print >> sys.stderr,  'Finished loading junctions.'
    
    keys = sorted( set( jns ) )
    for chrm, strand in keys:        
        if VERBOSE: print >> sys.stderr, \
                'Processing chromosome %s strand %s.' % ( chrm, strand )
        
        if strand != '+': continue
        
        disc_grpd_exons = find_exons_in_contig( \
           ( chrm, strand ),
           read_cov[ (chrm, strand) ],
           jns[ (chrm, strand) ],
           cage_cov[ (chrm, strand) ], 
           polya_sites[ (chrm, strand) ] )
    
    return
        
if __name__ == "__main__":
    main()
