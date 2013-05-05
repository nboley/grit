import sys
import numpy
from collections import defaultdict
from itertools import chain, izip

from igraph import Graph

from files.reads import RNAseqReads, CAGEReads, clean_chr_name, guess_strand_from_fname, iter_coverage_intervals_for_read
from files.junctions import extract_junctions_in_contig

def flatten( regions ):
    new_regions = []
    curr_start = regions[0][0]
    curr_end = regions[0][1]
    for i,(start,end) in enumerate(regions):
        if curr_end > end:
            end = curr_end
        if i+1 == len( regions ): 
            if len(new_regions) == 0:
                new_regions = [curr_start, curr_end]
                break
            if new_regions[-1][1] == end:
                break
            else:
                new_regions.append( [ curr_start, curr_end] )
                break
        if regions[i+1][0]-end <= 1:
            curr_end = max( regions[i+1][1], end ) 
        else:
            new_regions.append( [curr_start, curr_end] )
            curr_start = regions[i+1][0]
            curr_end = regions[i+1][1]
    if type(new_regions[0]) == int:
        return [new_regions]
    else:
        return new_regions

MIN_REGION_LEN = 200
EMPTY_BPK = 0
MIN_BPK = 1

def get_contigs_and_lens( rnaseq_reads, cage_reads ):
    chrm_lengths = {}
    for bam in chain( rnaseq_reads, cage_reads ):
        for ref_name, ref_len in zip(bam.references, bam.lengths):
            if clean_chr_name( ref_name ) not in chrm_lengths:
                chrm_lengths[clean_chr_name( ref_name )] = ref_len
            else:
                assert chrm_lengths[clean_chr_name(ref_name)] == ref_len, \
                    "Chromosome lengths do not match between bam files"
    
    return chrm_lengths

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

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(\
        description='Find exons from RNAseq, CAGE, and poly(A) assays.')

    parser.add_argument( 'rnaseq_reads',type=argparse.FileType('rb'),nargs='+',\
        help='BAM files containing mapped RNAseq reads ( must be indexed ).')
    
    parser.add_argument( '--cage-reads', type=file, default=[], nargs='*', \
        help='BAM files containing mapped cage reads.')
    
    parser.add_argument( '--polya-candidate-sites', type=file, nargs='*', \
        help='files with allowed polya sites.')
    
    parser.add_argument( '--out-filename', '-o', 
                         default="discovered_elements.bed",\
        help='Output file name. (default: discovered_elements.bed)')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true',\
        help='Whether or not to print status information.')
    parser.add_argument('--write-debug-data',default=False,action='store_true',\
        help='Whether or not to print out gff files containing intermediate exon assembly data.')
    parser.add_argument( '--threads', '-t', default=1, type=int,
        help='The number of threads to use.')
        
    args = parser.parse_args()

    global num_threads
    num_threads = args.threads
    
    global WRITE_DEBUG_DATA
    WRITE_DEBUG_DATA = args.write_debug_data
    
    global VERBOSE
    VERBOSE = args.verbose
        
    ofp = open( args.out_filename, "w" )
    
    return args.rnaseq_reads, args.cage_reads, args.polya_candidate_sites, ofp

class Bin( object ):
    def __init__( self, start, stop, left_label, right_label, 
                  bin_type=None, score=1000 ):
        self.start = start
        self.stop = stop
        assert stop - start > 0
        self.left_label = left_label
        self.right_label = right_label
        self.type = bin_type
        self.score = score
    
    def mean_cov( self, cov_array ):
        return cov_array[self.start:self.stop].mean()
    
    def reverse_strand(self, contig_len):
        return Bin(contig_len-self.stop, contig_len-self.start, 
                   self.right_label, self.left_label, self.type)
    
    def __repr__( self ):
        if self.type == None:
            return "%i-%i\t%s\t%s" % ( self.start, self.stop, self.left_label,
                                       self.right_label )

        return "%s:%i-%i" % ( self.type, self.start, self.stop )
    
    def __hash__( self ):
        return hash( (self.start, self.stop) )
    
    def __eq__( self, other ):
        return ( self.start == other.start and self.stop == other.stop )
    
    _bndry_color_mapping = {
        'CONTIG_BNDRY': '0,0,0',
        
        'POLYA': '255,255,0',

        'CAGE_PEAK': '0,255,0',
        
        'D_JN': '173,255,47',
        'R_JN': '0,0,255',
        
        'ESTART': '0,0,0',
        'ESTOP': '0,0,0'
    }
    
    def find_bndry_color( self, bndry ):
        return self._bndry_color_mapping[ bndry ]
    
    def _find_colors( self, strand ):
        if self.type != None:
            if self.type =='GENE':
                return '0,0,0'
            if self.type =='CAGE_PEAK':
                return '0,0,0'
            if self.type =='EXON':
                return '0,0,0'
            if self.type =='EXON_EXT':
                return '0,0,255'
            if self.type =='RETAINED_INTRON':
                return '255,255,0'
            if self.type =='TES_EXON':
                return '255,0,0'
            if self.type =='TSS_EXON':
                return '0,255,0'
            if self.type =='SE_GENE':
                return '255,255,0'

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
        if left_label == 'ESTART' and right_label  == 'ESTOP':
            return '159,153,87'
        
        return ( self.find_bndry_color(left_label), 
                 self.find_bndry_color(right_label) )

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
    
    def writeBed( self, ofp, contig_len, reverse_strand=True ):
        """
            chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
        """
        if reverse_strand and self.strand == '-':
            writetable_bins = self.reverse_strand( contig_len )
        else:
            writetable_bins = self
        
        for bin in writetable_bins:
            length = max( (bin.stop - bin.start)/4, 1)
            colors = bin._find_colors( self.strand )
            if isinstance( colors, str ):
                op = self._bed_template.format(
                    start=bin.start,stop=bin.stop+1,color=colors, 
                    name="%s_%s"%(bin.left_label, bin.right_label) )
                ofp.write( op )
            else:
                op = self._bed_template.format(
                    start=bin.start,stop=(bin.start + length),color=colors[0], 
                    name=bin.left_label)
                ofp.write( op )
                
                op = self._bed_template.format(
                    start=(bin.stop-length),stop=bin.stop,color=colors[1],
                    name=bin.right_label)
                ofp.write( op )
        
        return

    def writeGff( self, ofp, contig_len, filter=None ):
        """
            chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
        """
        if self.strand == '-':
            writetable_bins = self.reverse_strand( contig_len )
        else:
            writetable_bins = self
        
        for bin in writetable_bins:
            if filter != None and bin.type != filter:
                continue
            region = GenomicInterval(self.chrm, self.strand, 
                                     bin.start, bin.stop)
            grp_id = "%s_%s_%i_%i" % region
            ofp.write( create_gff_line(region, grp_id) + "\n" )
        
        return

def find_gene_boundaries((chrm, strand, contig_len), rnaseq_reads, polya_sites):
    def check_for_low_coverage_region( l_split, r_split ):
        region_len = r_split - l_split
        cnt, max_cnt = 0, region_len*EMPTY_BPK
        for rd in rnaseq_reads[0].iter_reads( chrm, strand, l_split, r_split ):
            cnt += 1
            if cnt > max_cnt: 
                return False
        return True
    
    def find_splits_in_region(rnaseq_reads, chrm, strand, start, stop):
        n_splits = int( (stop-start)/MIN_REGION_LEN )
        segment_len = int((stop-start)/(n_splits+1))
        return [ start + i*segment_len for i in xrange(1,n_splits+1) ]
    
    def find_initial_segmentation( chrm, strand, rnaseq_reads, polya_sites ):
        locs = {0: 'ESTOP', contig_len-1:'ESTART'}
        for polya in polya_sites:
            locs[ polya ] = "POLYA"

        for (start, stop), cnt in junctions:
            locs[start-1] = "D_JN"
            locs[stop+1] = "R_JN"

        boundaries = sorted( locs.keys() )
        for start, stop in zip( boundaries[:-1], boundaries[1:] ):
            region_len = stop - start
            if region_len > 3*MIN_REGION_LEN:
                splits = find_splits_in_region( 
                    rnaseq_reads, chrm, strand, start, stop )
                for l_split, r_split in izip( splits[:-1], splits[1:] ):
                    if check_for_low_coverage_region( l_split, r_split ):
                        locs[l_split] = 'ESTART' 
                        locs[r_split] = 'ESTOP' 
        
        boundaries = sorted( locs.keys() )
        bins = Bins( chrm, strand, [] )
        for index, (start, stop) in enumerate(
                izip(boundaries[:-1], boundaries[1:])):
            bins.append( Bin(start, stop, locs[start], locs[stop]) )
        
        del locs, boundaries
        return bins

    def find_gene_segments( bins ):
        gene_starts_indices = []
        for i, bin in enumerate( bins ):
            if bin.left_label == 'POLYA':
                if bin.right_label in ('D_JN', 'ESTART', 'ESTOP'):
                    gene_starts_indices.append( i )
                elif bin.right_label in 'POLYA':
                    gene_starts_indices[-1] = gene_starts_indices[-1] + 1
        
        gene_bndry_bins = Bins( chrm, strand )
        for start_i, stop_i in \
                zip( gene_starts_indices[:-1], gene_starts_indices[1:] ):
            start = bins[start_i].start
            left_label = bins[start_i].left_label
            stop = bins[stop_i-1].stop
            right_label = bins[stop_i-1].right_label
            gene_bndry_bins.append(
                Bin(start, stop, left_label, right_label, 'GENE'))
        
        gene_bndry_bins[0].start = 1
        gene_bndry_bins[-1].stop = contig_len - 1

        return gene_bndry_bins
    
    def merge_gene_segments(junctions, gene_bndry_bins):
        # find the junctions that overlap multiple gene bins
        # we use the interval overlap join algorithm
        sorted_jns = sorted( junctions )
        start_i = 0
        merge_jns = []
        for bin_i, bin in enumerate( gene_bndry_bins ):
            # increment the start pointer
            new_start_i = start_i
            for (jn_start, jn_stop), cnt in sorted_jns[start_i:]:
                if jn_stop >= bin.start:
                    break
                new_start_i += 1

            start_i = new_start_i

            # find matching junctions
            for (jn_start, jn_stop), cnt in sorted_jns[start_i:]:
                if jn_start > bin.stop:
                    break
                if jn_stop > bin.stop:
                    merge_jns.append( ( (jn_start, jn_stop), bin_i ) )

        genes_graph = Graph(len( gene_bndry_bins ))
        
        for jn, bin_i in merge_jns:
            # find the bin index that the jn merges into
            for end_bin_i, bin in enumerate( gene_bndry_bins[bin_i:] ):
                if jn[1]+1 < bin.start:
                    break
            
            assert bin_i + end_bin_i < len( gene_bndry_bins )
            for i in xrange( bin_i+1, bin_i+end_bin_i ):
                genes_graph.add_edge(bin_i, i)

        conn_nodes = genes_graph.clusters()
        new_bins = Bins( chrm, strand )
        for g in genes_graph.clusters():
            start_bin = gene_bndry_bins[ g[0] ]
            stop_bin = gene_bndry_bins[ g[-1] ]
            new_bins.append( 
                Bin( start_bin.start, stop_bin.stop, 
                     start_bin.left_label, stop_bin.right_label, "GENE")
            )
        
        return new_bins
    
    # find all of the junctions
    junctions = extract_junctions_in_contig( rnaseq_reads[0], chrm, strand )
    if VERBOSE: print "Finished extracting junctions for %s %s" % (chrm, strand)
    
    # find segment boundaries
    initial_segmentation = find_initial_segmentation( 
        chrm, strand, rnaseq_reads, polya_sites )
    if VERBOSE: print "Finished initial segmentation for %s %s" % (chrm, strand)
    
    gene_segments = find_gene_segments( initial_segmentation )
    if VERBOSE: print "Finished initial gene segments for %s %s" %(chrm, strand)
    
    if 0 == len( gene_segments  ):
        return gene_segments

    merged_gene_segments = merge_gene_segments( junctions, gene_segments )
    
    return merged_gene_segments
    

    gene_segments_graph = Graph(len(initial_segmentation)-1)
    gene_segments_graph.add_edges( all_edges )
        
        # flatten the segments
    gene_segments = []
    for g in gene_segments_graph.components():
        start = sorted_bnds[min(g)]
        stop = sorted_bnds[max(g)+1]
        gene_segments.append( (start, stop) )
        
    gene_segments.sort()
    return flatten( gene_segments )

    return initial_segmentation

    if 0 == len( gene_bndry_bins  ):
        return gene_bndry_bins

    # cluster the segments
    clustered_segments = cluster_segments( initial_segmentation )
    
    polya_sites = numpy.array( sorted( polya_sites ) )
    # merge the segments
    merged_segments = []
    for segment_i, (start, stop) in enumerate(clustered_segments[1:]):
        polyas = [x for x in polya_sites if x > stop - 1000 and x < stop + 1000]
        if len(polyas) == 0: continue
        merged_segments.append( (start, max(polyas)) )
    
    # build the gene bins, and write them out to the elements file
    genes = Bins( chrm, strand, [] )
    for start, stop in merged_segments:
        genes.append( Bin(start, stop, "ESTART", "POLYA", "GENE" ) )
    
    return genes

def find_exons_in_contig( (chrm, strand, contig_len), ofp,
                          rnaseq_reads, cage_reads, polya_sites):
    gene_bndry_bins = find_gene_boundaries( 
       (chrm, strand, contig_len), rnaseq_reads, polya_sites[(chrm, strand)])
    
    gene_bndry_bins.writeBed( ofp, contig_len, reverse_strand=False )
    
    return

def main():
    rnaseq_bams, cage_bams, polya_candidate_sites_fps, ofp \
        = parse_arguments()
    
    rnaseq_reads = [ RNAseqReads(fp.name).init(reverse_read_strand=False) 
                     for fp in rnaseq_bams ]
    
    cage_reads = [ CAGEReads(fp.name).init(reverse_read_strand=False) 
                   for fp in cage_bams ]

    if VERBOSE: print >> sys.stderr,  'Loading candidate polyA sites'
    polya_sites = find_polya_sites([x.name for x in polya_candidate_sites_fps])
    for fp in polya_candidate_sites_fps: fp.close()
    if VERBOSE: print >> sys.stderr, 'Finished loading candidate polyA sites'
    
    contig_lens = get_contigs_and_lens( rnaseq_reads, cage_reads )
    for contig, contig_len in contig_lens.iteritems():
        if contig != '4': continue
        for strand in '+':
            find_exons_in_contig( (contig, strand, contig_len), ofp,
                                  rnaseq_reads, cage_reads, polya_sites)
    
if __name__ == '__main__':
    main()



















#
