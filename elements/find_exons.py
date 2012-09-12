import sys, os
import numpy
import multiprocessing

from collections import OrderedDict, defaultdict
from itertools import izip

sys.path.append(os.path.join(os.path.dirname(__file__), "../", 'file_types'))
from wiggle import Wiggle, load_wiggle_asynchronously, guess_strand_from_fname
from junctions_file import parse_jn_gff, Junctions

class ThreadSafeFile( file ):
    def __init__( self, fname, mode, trackname=None ):
        file.__init__( self, fname, mode )
        self._writelock = multiprocessing.Lock()
        if trackname != None:
            self.write('track name="%s" visibility=2 itemRgb="On"\n'%trackname)
    
    def write( self, data ):
        self._writelock.acquire()
        file.write( self, data )
        file.flush( self )
        self._writelock.release()


class GroupIdCntr( object ):
    def __init__( self, start_val=0 ):
        self.value = multiprocessing.RawValue( 'i', start_val )
        self._lock = multiprocessing.RLock()

    def __iadd__(self, val_to_add):
        self._lock.acquire()
        self.value.value += val_to_add
        self._lock.release()
        return self
    
    def lock(self):
        self._lock.acquire()

    def release(self):
        self._lock.release()
    
    def __str__( self ):
        self._lock.acquire()
        rv = str(self.value.value)
        self._lock.release()
        return rv


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
    
    def reverse_strand( self, contig_len ):
        rev_bins = Bins( self.chrm, self.strand )
        for bin in reversed(self):
            rev_bins.append( bin.reverse_strand( contig_len ) )
        return rev_bins
    
    def writeBedgraph( self, ofp, contig_len ):
        """
            chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
        """
        if self.strand == '-':
            writetable_bins = self.reverse_strand( contig_len )
        else:
            writetable_bins = self
        
        for bin in writetable_bins:
            length = max( (bin.stop - bin.start)/4, 1)
            colors = bin._find_colors( self.strand )
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
    def __init__( self, start, stop, left_label, right_label, bin_type=None ):
        self.start = start
        self.stop = stop
        self.left_label = left_label
        self.right_label = right_label
        self.type = bin_type
    
    def reverse_strand(self, contig_len):
        return Bin(contig_len-self.stop, contig_len-self.start, 
                   self.right_label, self.left_label, self.type)
    
    def __repr__( self ):
        if self.type == None:
            return "%i-%i\t%s\t%s" % ( self.start, self.stop, self.left_label,
                                       self.right_label )

        return "%i-%i\t%s" % ( self.start, self.stop, self.type )
    
    _bndry_color_mapping = {
        'CONTIG_BNDRY': '0,0,0',
        
        'POLYA': '255,255,0',
        
        'D_JN': '173,255,47',
        'R_JN': '0,255,0',
    }
    
    def find_bndry_color( self, bndry ):
        return self._bndry_color_mapping[ bndry ]
    
    def _find_colors( self, strand ):
        if self.type != None:
            if self.type =='GENE':
                return '0,0,0'
            if self.type =='CAGE_PEAK':
                return '0,0,0'

        if strand == '+':
            left_label, right_label = self.left_label, self.right_label
        else:
            assert strand == '-'
            left_label, right_label = self.right_label, self.left_label
        
        
        if left_label == 'D_JN' and right_label  == 'R_JN':
            return '108,108,108'
        if left_label == 'D_JN' and right_label  == 'D_JN':
            return '135,206,250'
        if left_label == 'R_JN' and right_label  == 'R_JN':
            return '135,206,250'
        if left_label == 'R_JN' and right_label  == 'D_JN':
            return '0,0,255'
        if left_label == 'R_JN' and right_label  == 'POLYA':
            return '255,0,0'
        if left_label == 'POLYA' and right_label  == 'POLYA':
            return ' 240,128,128'
        if left_label == 'D_JN' and right_label  == 'POLYA':
            return '240,128,128'
        if left_label == 'POLYA' and right_label  == 'D_JN':
            return '147,112,219'
        if left_label == 'POLYA' and right_label  == 'R_JN':
            return '159,153,87'
        
        return ( self.find_bndry_color(left_label), 
                 self.find_bndry_color(right_label) )
        

def reverse_contig_data( rnaseq_cov, jns, cage_cov, polya_sites ):
    assert len( rnaseq_cov ) == len( cage_cov )
    genome_len = len( rnaseq_cov )
    rev_jns = dict( ((genome_len-stop-1, genome_len-start-1), value) 
                for (start, stop), value in jns.iteritems() )
    return ( rnaseq_cov[::-1], rev_jns, cage_cov[::-1], genome_len-polya_sites )

def check_for_se_gene( start, stop, cage_cov, rnaseq_cov ):
    # check to see if the rnaseq signal after the cage peak
    # is much greater than before
    split_pos = cage_cov[start:stop+1].argmax() + start
    if cage_cov[split_pos-1:split_pos+2].sum() < 20:
        return False
    window_size = min( split_pos - start, stop - split_pos, 200 )
    after_cov = rnaseq_cov[split_pos:split_pos+window_size].sum()
    before_cov = rnaseq_cov[split_pos-window_size:split_pos].sum()
    if (after_cov/(before_cov+1)) > 5:
        return True
 
    return False

def find_gene_boundaries( (chrm, strand), bins, cage_cov, rnaseq_cov, jns ):
    # find regions that end with a polya, and look like genic regions
    # ( ie, they are a double polya that looks like a single exon gene, or
    # they start after a polya leading into a donor exon
    gene_starts_indices = []
    for i, bin in enumerate( bins ):
        if bin.left_label == 'POLYA' and bin.right_label  == 'D_JN':
            gene_starts_indices.append( i )
        if bin.left_label == 'POLYA' and bin.right_label  == 'POLYA':
            if check_for_se_gene( bin.start, bin.stop, cage_cov, rnaseq_cov ):
                gene_starts_indices.append( i )
    
    # build a bins object from the initial labels
    gene_bndry_bins = Bins( chrm, strand )
    for start_i, stop_i in \
            zip( gene_starts_indices[:-1], gene_starts_indices[1:] ):
        start = bins[start_i].start
        left_label = bins[start_i].left_label
        stop = bins[stop_i-1].stop
        right_label = bins[stop_i-1].right_label
        gene_bndry_bins.append(
            Bin(start, stop, left_label, right_label, 'GENE'))

    if 0 == len( gene_bndry_bins  ):
        return gene_bndry_bins
    
    gene_bndry_bins[0].start = 1
    gene_bndry_bins[-1].stop = len( rnaseq_cov )-1
    
    # find the junctions that overlap multiple gene bins
    sorted_jns = sorted( jns.keys() )
    start_i = 0
    merge_jns = []
    for bin_i, bin in enumerate( gene_bndry_bins ):
        # increment the start pointer
        for jn in sorted_jns[start_i:]:
            if jn[1] >= bin.start:
                break
            start_i += 1
        
        # find matching junctions
        for jn in sorted_jns[start_i:]:
            if jn[0] > bin.stop:
                break
            if jn[1] >= bin.stop:
                merge_jns.append( ( jn, bin_i ) )
    
    joined_bins_mapping = defaultdict( int )
    for jn, bin_i in merge_jns:
        # find the bin index that the jn merges into
        for end_bin_i, bin in enumerate( gene_bndry_bins[bin_i:] ):
            if jn[1] < bin.start:
                break
        # take the largest bin merge. This means that we will merge over genes 
        # if they exist int he middle
        joined_bins_mapping[bin_i] = max( 
            bin_i+end_bin_i-1, joined_bins_mapping[bin_i] )
    
    for start_bin, stop_bin in reversed(sorted(joined_bins_mapping.iteritems())):
        gene_bndry_bins[start_bin].stop = gene_bndry_bins[stop_bin].stop
        for x in xrange( start_bin+1, stop_bin+1 ):
            del gene_bndry_bins[x]
    
    return gene_bndry_bins

def find_peaks( cov, window_len, min_score, max_score_frac, max_num_peaks ):
    cumsum_cvg_array = \
        numpy.append(0, numpy.cumsum( cov ))
    scores = cumsum_cvg_array[window_len:] - cumsum_cvg_array[:-window_len]
    indices = numpy.argsort( scores )
    
    def overlaps_prev_peak( new_loc ):
        for start, stop in peaks:
            if not( new_loc > stop or new_loc + window_len < start ):
                return True
        return False
    
    # merge the peaks
    def grow_peak( start, stop, grow_size=max(3, window_len/4), min_grow_ratio=0.5 ):
        # grow a peak at most max_num_peaks times
        for i in xrange(max_num_peaks):
            curr_signal = cov[start:stop+1].sum()
            downstream_sig = cov[max(0, start-grow_size):start].sum()
            upstream_sig = cov[stop+1:stop+1+grow_size].sum()
            exp_factor = float( stop - start + 1 )/grow_size
            
            # if neither passes the threshold, then return the current peak
            if float(max( upstream_sig, downstream_sig ))*exp_factor \
                    < curr_signal*min_grow_ratio: return (start, stop)
            
            # otherwise, we know one does
            if upstream_sig > downstream_sig:
                stop += grow_size
            else:
                start = max(0, start - grow_size )
        
        if VERBOSE:
            print "Warning: reached max peak iteration at %i-%i ( signal %.2f )"\
                % (start, stop, cov[start:stop+1].sum() )
            print 
        return (start, stop )
    
    peaks = []
    peak_scores = []
    
    for index in reversed(indices):
        if not overlaps_prev_peak( index ):
            score = scores[ index ]
            new_peak = grow_peak( index, index + window_len )
            # if we are below the minimum score, then we are done
            if score < min_score:
                break

            # if we have observed peaks, and the ratio between the highest
            # and the lowest is sufficeintly high, we are done
            if len( peak_scores ) > 0:
                if float(score)/peak_scores[0] < max_score_frac:
                    break
                        
            peaks.append( new_peak ) 
            peak_scores.append( score )
    
    if len( peaks ) == 0:
        return []
    
    # merge cage peaks together
    def merge_peaks( peaks_and_scores ):
        merged_peaks = set()
        new_peaks = []
        new_scores = []
        for pk_i, (peak, score) in enumerate(peaks_and_scores):
            if pk_i in merged_peaks: continue
            curr_pk = list( peak )
            curr_score = score
            for i_pk_i, (i_peak, i_score) in enumerate(peaks_and_scores):
                if i_pk_i in merged_peaks: continue
                if i_peak[0] < curr_pk[0]: continue
                if i_peak[0] - curr_pk[1] < max( window_len, 
                                                 curr_pk[1]-curr_pk[0] ):
                    curr_pk[1] = i_peak[1]
                    curr_score += i_score
                    merged_peaks.add( i_pk_i )
                else:
                    break

            new_peaks.append( curr_pk )
            new_scores.append( curr_score )
        return zip( new_peaks, new_scores )
    
    peaks_and_scores = sorted( zip(peaks, peak_scores) )
    old_len = len( peaks_and_scores )
    for i in xrange( 99 ):
        if i == 100: assert False
        peaks_and_scores = merge_peaks( peaks_and_scores )
        if len( peaks_and_scores ) == old_len: break
    
    max_score = max( s for p, s in peaks_and_scores )
    return [ pk for pk, score in peaks_and_scores \
                 if score/max_score > max_score_frac ]


def find_cage_peaks_in_gene( ( chrm, strand ), gene, cage_cov, rnaseq_cov ):
    raw_peaks = find_peaks( cage_cov[gene.start:gene.stop+1], 
                            window_len=20, min_score=20, 
                            max_score_frac=0.10, max_num_peaks=20 )
    if len( raw_peaks ) == 0:
        print >> sys.stderr, "WARNING: Can't find peak in region %s:%i-%i" % \
            ( chrm, gene.start, gene.stop )
        return []
    
    cage_peaks = Bins( chrm, strand )
    for peak_st, peak_sp in raw_peaks:
        cage_peaks.append( Bin( peak_st+gene.start, peak_sp+gene.start+1,
                                "CAGE_PEAK_START", "CAGE_PEAK_STOP", "CAGE_PEAK") )
    return cage_peaks


def find_bins_in_contig( ( chrm, strand ), rnaseq_cov, jns, cage_cov, polya_sites ):
    if strand == '-':
        rnaseq_cov, jns, cage_cov, polya_sites = reverse_contig_data( 
            rnaseq_cov, jns, cage_cov, polya_sites )
        
    locs = {}
    for polya in polya_sites:
        locs[ polya ] = "POLYA"
    
    for start, stop in jns:
        locs[start-1] = "D_JN"
        locs[stop+1] = "R_JN"

    # build all of the bins
    poss = sorted( locs.iteritems() )
    bins = Bins( chrm, strand, [ Bin(1, poss[0][0], "CONTIG_BNDRY", poss[0][1]), ])
    for index, ((start, left_label), (stop, right_label)) in \
            enumerate(izip(poss[:-1], poss[1:])):
        bins.append( Bin(start, stop, left_label, right_label) )
    
    bins.append( Bin( poss[-1][0], len(rnaseq_cov)-1, poss[-1][1], "CONTIG_BNDRY" ) )
    bins.writeBedgraph( binsFps[strand], len(rnaseq_cov) )
    
    # find gene regions
    # first, find the gene bin indices
    gene_bndry_bins = find_gene_boundaries( 
        (chrm, strand), bins, cage_cov, rnaseq_cov, jns )
    
    # find cage peaks
    refined_gene_bndry_bins = Bins( chrm, strand, [] )
    cage_peaks = Bins( chrm, strand )
    for gene_bin in gene_bndry_bins:
        gene_cage_peaks = find_cage_peaks_in_gene( 
                (chrm, strand), gene_bin, cage_cov, rnaseq_cov )
        cage_peaks.extend( gene_cage_peaks )
        if len( gene_cage_peaks ) > 0:
            start, stop = gene_bin.start, gene_bin.stop
            start = gene_cage_peaks[0].start
            refined_gene_bndry_bins.append( 
                Bin( start, stop, 
                     gene_bin.left_label, gene_bin.right_label, "GENE") )
    
    cage_peaks.writeBedgraph( cagePeaksFps[strand], len(rnaseq_cov) )
    refined_gene_bndry_bins.writeBedgraph( 
        geneBoundariesFps[strand], len(rnaseq_cov) )
    
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
    
    # prepare the intermediate output objects
    global binsFps, geneBoundariesFps, cagePeaksFps
    binsFps = { "+": ThreadSafeFile( 
            args.out_file_prefix + "bins.plus.gff", "w", "bins_plus" ),
                "-": ThreadSafeFile( 
            args.out_file_prefix + "bins.minus.gff", "w", "bins_minus" ) }
    geneBoundariesFps = { 
        "+": ThreadSafeFile( 
            args.out_file_prefix + "gene_boundaries.plus.gff", "w", "gene_bndrys_plus" ),
        "-": ThreadSafeFile( 
            args.out_file_prefix + "gene_boundaries.minus.gff", "w", "gene_bndrys_minus" ) 
    }
    cagePeaksFps = { 
        "+": ThreadSafeFile( 
            args.out_file_prefix + "cage_peaks.plus.gff", "w", "cage_peaks_plus" ),
        "-": ThreadSafeFile(
            args.out_file_prefix + "cage_peaks.minus.gff", "w", "cage_peaks_minus" ) 
    }


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
    
    rnaseq_grpd_wigs = [ rd1_plus_wigs, rd1_minus_wigs, rd2_plus_wigs, rd2_minus_wigs ]
    
    cage_plus_wigs = [ fp for fp in args.cage_wigs 
                      if fp.name.endswith("+.bedGraph") ]
    cage_minus_wigs = [ fp for fp in args.cage_wigs 
                      if fp.name.endswith("-.bedGraph") ]
    cage_grpd_wigs = [ cage_plus_wigs, cage_minus_wigs ]
    
    return rnaseq_grpd_wigs, args.junctions, args.chrm_sizes_fname, \
        cage_grpd_wigs, args.polya_candidate_sites, ofps


def main():
    wigs, jns_fp, chrm_sizes_fp, cage_wigs, polya_candidate_sites_fps, out_fps \
        = parse_arguments()
    
    group_id_starts = [ GroupIdCntr(1) for fp in out_fps ]
    
    # set up all of the file processing calls
    worker_pool = multiprocessing.BoundedSemaphore( num_threads )
    ps = []
    
    if VERBOSE: print >> sys.stderr,  'Loading merged read pair wiggles'    
    fnames =[wigs[0][0].name, wigs[1][0].name, wigs[2][0].name, wigs[3][0].name]
    strands = ["+", "-", "+", "-"]
    read_cov, new_ps = load_wiggle_asynchronously( 
        chrm_sizes_fp.name, fnames, strands, worker_pool )
    ps.extend( new_ps )    
    
    if VERBOSE: print >> sys.stderr,  'Loading CAGE.'
    assert all( len(fps) == 1 for fps in cage_wigs )
    cage_cov, new_ps = load_wiggle_asynchronously( 
        chrm_sizes_fp.name, [fps[0].name for fps in cage_wigs], ["+", "-"], worker_pool)
    ps.extend( new_ps )
            
    if VERBOSE: print >> sys.stderr,  'Loading candidate polyA sites'
    polya_sites = find_polya_sites([x.name for x in polya_candidate_sites_fps])
    for fp in polya_candidate_sites_fps: fp.close()
    if VERBOSE: print >> sys.stderr, 'Finished loading candidate polyA sites'

    if VERBOSE: print >> sys.stderr,  'Loading junctions.'
    jns = Junctions( parse_jn_gff( jns_fp.name ) )
    
    for p in ps:
        if VERBOSE:
            print "Joining process: ", p
        p.join()
        if VERBOSE:
            print "Joined process: ", p
    
    assert all( x.min() >= 0 and x.max() >= 0 for x in cage_cov.values() )
    assert all( x.min() >= 0 and x.max() >= 0 for x in read_cov.values() )

    ##### join all of the wiggle processes
    
    # get the joined read data, and calculate stats
    if VERBOSE: print >> sys.stderr, 'Finished loading merged read pair wiggles'
    
    # join the cage data, and calculate statistics
    for cage_fps in cage_wigs: [ cage_fp.close() for cage_fp in cage_fps ]
    if VERBOSE: print >> sys.stderr, 'Finished loading CAGE data'

    if VERBOSE: print >> sys.stderr, 'Finished loading read pair 1 wiggles'

    if VERBOSE: print >> sys.stderr, 'Finished loading read pair 2 wiggles'
    
    keys = sorted( set( jns ) )
    for chrm, strand in keys:        
        if VERBOSE: print >> sys.stderr, \
                'Processing chromosome %s strand %s.' % ( chrm, strand )
        
        bins = find_bins_in_contig( \
           ( chrm, strand ),
           read_cov[ (chrm, strand) ],
           jns[ (chrm, strand) ],
           cage_cov[ (chrm, strand) ], 
           polya_sites[ (chrm, strand) ] )        


    for fps_dict in (binsFps, geneBoundariesFps, cagePeaksFps):
        for fp in fps_dict.values():
            fp.close()
    
    return
        
if __name__ == "__main__":
    main()
