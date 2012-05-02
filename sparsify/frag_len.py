#!/usr/bin/python

import numpy
from numpy import mean, std, append
import operator
import sys
from collections import defaultdict
import copy

import os

VERBOSE = False
MIN_FLS_FOR_FL_DIST = 100

def get_fl_dist_base_from_fnames( reads_fname, genes_fname, cluster_read_groups=True ):
    """Get the fl dist base filename. 
        
    This is used for the cached fldist and fldist analysis objects.
    """
    def get_hash_from_fname( fname ):
        fsize = os.path.getsize( fname )
        hash_key = os.path.basename(fname) + "_" + str( fsize  )
        return hash_key
        
    reads_hash_key = get_hash_from_fname( reads_fname )
    exons_hash_key = get_hash_from_fname( genes_fname )    
    return ".fldist_" + reads_hash_key + "_" + exons_hash_key \
        + "_%i" % cluster_read_groups

    
def find_diagnostic_plot_fname( reads_fname, genes_fname, cluster_read_groups=True ):
    fname = get_fl_dist_base_from_fnames( reads_fname, genes_fname, cluster_read_groups )
    return fname + ".pdf"

def generate_gtf( reads, genes , out, rmsk=None ):
    """Generate a gtf file with num_exons_to_generate single exon genes with highest read depth

    The exons will not overlap a repeat region if a repeat masker 
    bed file is provided.
    """
    num_exons_to_generate = 250
    
    exon_read_depths = []
    for gene_mo in genes.values():
        # skip genes with more than 1 exon or short exons
        if len( gene_mo.exon_bndrys ) != 1:
            continue
        if (gene_mo.exon_bndrys[0][1] - gene_mo.exon_bndrys[0][0]) < 700:
            continue
        
        # calculate the length of the exon that we are looking at
        exon_len = gene_mo.exon_bndrys[0][1] - gene_mo.exon_bndrys[0][0]
        
        # iterate through read pairs to count
        rd_pr_count = 0
        for read1, read2 in reads.iter_paired_reads( gene_mo.boundaries ):
            # skip all secondary alignments
            if read1.is_secondary or read2.is_secondary:
                continue
            rd_pr_count += 1
            
        exon_read_depths.append( ( rd_pr_count/exon_len , gene_mo) )

    #sort exon_read_depths by read_depth
    def get_read_depth( obj ): return obj[0]
    exon_read_depths.sort(key=get_read_depth)

    exons2print = []
    # use repeat masker file to remove exons in low mapability regions
    # if repeat masker bed file is provided 
    if rmsk != None:
        import subprocess, shlex
        # loop through exons sorted by read_depth
        for exon_read_depth in exon_read_depths:
            # write single line exon bed file for each exon to be tested
            exon_fh = open( '.single_exon_for_high_read_depth_gtf_generator.txt' ,'w' )
            gene_mo = exon_read_depth[1]
            bed_line = "chr" + gene_mo.chromosome + "\t" + \
                str(gene_mo.exon_bndrys[0][0]) + "\t" + str(gene_mo.exon_bndrys[0][1])
            exon_fh.write( bed_line )
            exon_fh.close()
            # run bedtools intersect to see if exon is in repeat region
            cmd = shlex.split('intersectBed -u -b ' + rmsk + \
                                  ' -a .single_exon_for_high_read_depth_gtf_generator.txt')
            P = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            exon_in_rpt = False
            for line in P.stdout:
                exon_in_rpt = True
            # if the exon does not overlap a repeat region
            if not exon_in_rpt:
                exons2print.append( exon_read_depth )
            # when we have found num_exons_to_generate exons to print out break loop
            if len( exons2print ) >= num_exons_to_generate: 
                break
    else:
        exons2print = exon_read_depths[:num_exons_to_generate]

    def write_gtf_line( exon ):
        # return the following information in tab-delimited format as in a gtf file
        # (seqname, source, feature, start, end, score, strand, frame, group)
        line = [ exon.chromosome, 'frag_len.py_generated', 'exon', \
                     str(exon.bp_bndry(0)[0]), str(exon.bp_bndry(0)[1]), \
                     ".", exon.strand, '.', 'gene_id "' + exon.name + \
                     '"; transcript_id "' + exon.name + \
                     '"; exon_number "1"\n']
        return '\t'.join( line )

    # write out gtf lines corresponding to exons with highest read depth
    fid = open( out , 'w' )
    for exon in zip(*exons2print)[1]:
        fid.write( write_gtf_line( exon ) )
    fid.close()

    return

class FlDist( object ):
    """Store a fragment length dist.

    """
    """fl_min is the minimum allowed fragment length
       fl_max is the maximum allowed fragment length
       fl_density is a vector of length fl_max - fl_min + 1,
          where each entry is the fraction of fragments of that
          length. ie, the vector should sum to 1.

    """
    def __init__( self, fl_min, fl_max, fl_density, stats=None ):
        assert fl_min <= fl_max
        try:
            assert type( fl_min ) == int
            assert type( fl_max ) == int
            assert len( fl_density ) == ( fl_max - fl_min + 1 )
        except:
            print fl_min, fl_max
            print fl_density
            raise
        
        self.fl_min = fl_min
        self.fl_max = fl_max
        self.fl_density = fl_density
        self.fl_density_cumsum = fl_density.cumsum()
        # cumsum weighted by the fragment length, just a caching optimization
        self.fl_density_weighted_cumsum = \
            self.fl_density*numpy.arange( fl_min, fl_max+1 )
        self.fl_density_weighted_cumsum = self.fl_density_weighted_cumsum.cumsum()

        # build and set the hash value
        self._hash_value = hash( ( self.fl_min, self.fl_max, \
                                   tuple( self.fl_density ) ) )
        self.stats = stats
    
    def __eq__( self, other ):
        if self.fl_min == other.fl_min \
            and self.fl_max == other.fl_max \
            and self.fl_density == other.fl_density:
            return True
        else:
            return False
        assert False
    
    def plot( self ):
        import matplotlib
        import matplotlib.pyplot as plt
        plt.scatter( range(self.fl_min, self.fl_max+1), self.fl_density )
        plt.show()
    
    def __hash__( self ):
        return self._hash_value
            
    
def build_uniform_density( fl_min, fl_max ):
    length = fl_max - fl_min + 1
    return FlDist( fl_min, fl_max,  numpy.ones( length )/length )

def build_normal_density( fl_min, fl_max, mean, sd ):
    length = fl_max - fl_min + 1
    from scipy.stats import norm
    density = norm(loc=mean, scale=sd)
    density_array = numpy.array( [ density.pdf( index ) \
                                   for index in xrange( fl_min, fl_max+1 ) ] )
    density_array = density_array/density_array.sum()
    return FlDist( fl_min, fl_max, density_array )

def build_multi_modal_density( mean_sd_freq ):
    """This function is for creating an accurate fragment length distribution from several normal densities (i.e. simulated Illumina and 454)
    Note: This function has not been tested!!!
    """
    # set fl_min and fl_max and create scipy normal densities for each
    fl_min = mean_sd_freq[0][0]
    fl_max = mean_sd_freq[0][0]
    from scipy.stats import norm
    norms = []
    for stats in mean_sd_freq:
        norms.append( norm( loc=stats[0], scale=stats[1] ) )
        if density[0] - (4 * density[1]) < fl_min:
            fl_min = density[0] - (4 * density[1])
        if density[0] + (4 * density[1]) < fl_max:
            fl_max = density[0] + (4 * density[1])

    density_array = numpy.zeros( fl_max - fl_min + 1 )
    for index in xrange( fl_min, fl_max + 1):
        for norms_index in len(mean_sd_freq):
            density_array[ index ] += \
                norms[ norms_index ].pdf( index ) * mean_sd_freq[2]

    return FlDist( fl_min, fl_max, density_array )

def group_fragments_by_readgroup( fragments ):
    grouped_fragments = {}
    for ( read_group, strand, exon_len ), fl in fragments:
        try:
            grouped_fragments[ read_group ].append( fl )
        except KeyError:
            grouped_fragments[ read_group ] = [ fl, ]
    
    # make a numpy array out of each entry, for easier access to the cumdist, etc.
    for read_group in grouped_fragments.keys():
        grouped_fragments[ read_group ] \
            = numpy.array( grouped_fragments[ read_group ] )
    
    return grouped_fragments

def group_fragments_by_strand( fragments ):
    grouped_fragments = {}
    for ( read_group, strand, exon_len ), fl in fragments:
        try:
            grouped_fragments[ strand ].append( fl )
        except KeyError:
            grouped_fragments[ strand ] = [ fl, ]
    
    # make a numpy array out of each entry, for easier access to the cumdist, etc.
    for strand in grouped_fragments.keys():
        grouped_fragments[ strand ] \
            = numpy.array( grouped_fragments[ strand ] )
    
    return grouped_fragments

def group_fragments_by_exonlen( fragments ):
    # sort fragments by exon_len, randomizing for identical exon lengths
    import random
    def get_exon_len( obj ): return obj[0][2] + random.random()
    fragments.sort(key=get_exon_len)
    sorted_frags = zip(*fragments)[1]

    grouped_fragments = []
    for i in xrange(10):
        grouped_fragments.append( numpy.array( sorted_frags[ \
                    int(i*len(sorted_frags)/10) : int((i+1)*len(sorted_frags)/10) ] ) )

    return grouped_fragments

def analyze_fl_dists( fragments, out_filename='diagnostic_plots'):
    """Analyze the fragment lengths. 

    Produce diagnostic plots and statistics for fragment length distributions 
    separated by read group, exon length and strand.
    plots are output into a pdf with filename out_filename
    """
    # open pdf for showing multiple diagnostic plots
    import matplotlib
    matplotlib.use('PDF')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(out_filename)

    # Compare fragment lengths to the read group from which they came
    read_grouped_frag_lengths = group_fragments_by_readgroup( fragments )
    grouped_frag_lengths_for_bxplt = []
    for read_group, rg_frag_lengths in read_grouped_frag_lengths.iteritems():
        # do not plot read_groups which do not have enough data for significant
        # plots to be created adn which were not clustered for the same reason
        if len( rg_frag_lengths ) < MIN_FLS_FOR_FL_DIST:
            continue
        fl_dist = build_robust_fl_dist_with_stats( rg_frag_lengths )
        grouped_frag_lengths_for_bxplt.append( rg_frag_lengths[ ( \
                    rg_frag_lengths > fl_dist.fl_min) \
                    & (rg_frag_lengths < fl_dist.fl_max) ] )
    # plot separate box plots for each read group with trimmed data
    plt.boxplot( grouped_frag_lengths_for_bxplt )
    plt.xlabel('Read Group')
    plt.ylabel('Fragment Size')
    plt.suptitle('Fragement Size Distribution Separated by Read Group')
    pp.savefig()
    plt.clf()

    # Compare fragment lengths to the exons from which they came
    exon_len_grouped_frag_lengths = group_fragments_by_exonlen( fragments )
    grouped_frag_lengths_for_bxplt = []
    #print 'Exon length binned statistics(mean,sd,skew,se_m,se_sd,se_sk,low,up):'
    for exon_binned_frag_lengths in exon_len_grouped_frag_lengths:
        fl_dist = build_robust_fl_dist_with_stats( exon_binned_frag_lengths )
        grouped_frag_lengths_for_bxplt.append( exon_binned_frag_lengths[ \
                (exon_binned_frag_lengths > fl_dist.fl_min) \
                    & (exon_binned_frag_lengths < fl_dist.fl_max) ] )
    # plot separate box plots for each exon length with trimmed data
    plt.boxplot( grouped_frag_lengths_for_bxplt )
    plt.xlabel('Exon Length Bin (Short to Long)')
    plt.ylabel('Fragment Size')
    plt.suptitle('Fragement Size Distribution Separated by Exon Length')
    pp.savefig()
    plt.clf()

    # Compare fl_dists separated by strand
    strand_grouped_frag_lengths = group_fragments_by_strand( fragments )
    strand_fl_dists = {}
    for strand , strand_grouped_frag_lengths in strand_grouped_frag_lengths.iteritems():
        #print strand , 'strand statistics(mean,sd,skew,se_m,se_sd,se_sk,low,up):'
        strand_fl_dists[ strand ] = build_robust_fl_dist_with_stats( strand_grouped_frag_lengths )
    # plot scatter plots on the same plot for both forward and reverse strand.
    plt.scatter( range(strand_fl_dists['+'].fl_min, strand_fl_dists['+'].fl_max+1), \
                     strand_fl_dists['+'].fl_density , s=10, c='b' )
    plt.scatter( range(strand_fl_dists['-'].fl_min, strand_fl_dists['-'].fl_max+1), \
                     strand_fl_dists['-'].fl_density , s=10, c='r' )
    plt.xlabel('Fragment Size')
    plt.ylabel('Proportion of Fragments')
    plt.suptitle('Fragement Size Distribution Separated by Strand')
    pp.savefig()
    plt.clf()

    pp.close()

    return

def find_fragments( reads, genes ):
    fragment_sizes = []
    for gene_mo in genes.values():
        # skip genes with more than 1 exon or short exons
        if len( gene_mo.exon_bndrys ) != 1:
            continue
        if (gene_mo.exon_bndrys[0][1] - gene_mo.exon_bndrys[0][0]) < 700:
            continue
                
        # calculate the length of the exon
        exon_len = gene_mo.exon_bndrys[0][1] - gene_mo.exon_bndrys[0][0]
        
        # iterate through read pairs in exon to get fragment lengths
        for read1, read2 in reads.iter_paired_reads( gene_mo.boundaries ):
            # skip all secondary alignments
            if read1.is_secondary or read2.is_secondary:
                continue

            # get read group from tags
            # put all frags w/o a read group into mean read group
            try:
                read_group = [ val for key, val in read1.tags if key == 'RG' ][0]
            except IndexError:
                read_group = 'mean'
                
            strand = '-' if read1.is_reverse else '+'
            
            if not read1.is_reverse:
                frag_len = read2.aend - read1.pos
            else:
                frag_len = read1.aend - read2.pos
    
            key = ( read_group, strand, exon_len )
            fragment_sizes.append( ( key, frag_len ) )
        
    return fragment_sizes

def build_robust_fl_dist_with_stats( fragment_lengths ):
    """Trim outliers from a numpy array of fragment lengths.

    First, we estimate the mean and sd on the trimmed data. Then we use the initial
    estimate to truncate at +/- NUM_SDS SD's
    
    Calculate statistics for use in read_group distribution clustering
    """
    NUM_SDS = 4
    
    sorted_fls = numpy.sort( fragment_lengths )
    
    ## find the initial estiamte of the mean and sd
    # if the list is long, trim the top and bottom 5%
    if len( sorted_fls ) > 40:
        fifteen_percent_cnt = int( len( sorted_fls )*.15 )
        guess = sorted_fls[ fifteen_percent_cnt:-fifteen_percent_cnt  ]
        mn = numpy.mean( guess )
        sd = numpy.std( guess )
        lower_bnd = max( 0, mn - NUM_SDS*sd )
        upper_bnd = mn + NUM_SDS*sd
    # if the list is too short, just include everything
    else:
        lower_bnd = sorted_fls[0]
        upper_bnd = sorted_fls[-1]

    new_fls = sorted_fls[ (sorted_fls > lower_bnd) & (sorted_fls < upper_bnd) ]
    if len( new_fls ) == 0:
        print sorted_fls
        print lower_bnd
        print upper_bnd
        raise
    lower_bnd = new_fls[0]
    upper_bnd = new_fls[-1]

    # calculate mean, standard deviation and skew as well as associated SEs
    from scipy.stats import skew, moment, sem
    from math import sqrt
    mn = numpy.mean( new_fls )
    sd = numpy.std( new_fls )
    skw = skew( new_fls )
    n = len( new_fls )
    # calculate SE for mean, sd and skew
    se_m = sd / sqrt(n)
    se_sd = sqrt( (2 / float(n**2 - n)) * moment( new_fls, 4) )
    se_sk = sqrt( float( 6*n*(n-1) )/float( (n-2) * (n+1) * (n-3) ) )
    #assert skw - se_sk < 0, 'skew: %f skew_standard_error: %f' % (skew, se_sk)
    stats = ( (mn, sd, skw) , (se_m, se_sd, se_sk) )

    # build the empirical distribution
    emp_dist = numpy.zeros( upper_bnd - lower_bnd + 1 )
    for fl in new_fls:
        emp_dist[ fl - lower_bnd ] += 1
    emp_dist = emp_dist/emp_dist.sum()

    # build the fl dist object
    fl_dist = FlDist( int(lower_bnd), int(upper_bnd), emp_dist, stats )
    
    return fl_dist

def cluster_rdgrps( stats ):
    """

    stats is a dictioanry of readgroups, keyed by rg name, with values
    that take tuples of summary statistics.
    """
    def are_diff( val1, val2, tol ):
        """We say two things are different if they're off by > 10%
        
        """
        if abs( val1*0.10 ) > tol: tol = abs( val1 )*0.10
        if abs(val1 - val2) < tol: return False
        else: return True

    def grp_sorted_data( data ):
        """Return a mapping of the keys in data ( the first item ) 
           to 'groups' which are ints from 0 to num clusters.
        """
        grps = {}
        grp_key = 0
        prev_key, prev_val, prev_tol = data[0]
        for key, val, tol in data:
            if are_diff( val, prev_val, tol ):
                grp_key += 1
            grps[ key ] = grp_key
        return grps
    
    def cluster_by_summary_stat( stat_index ):
        def sort_key( obj ): return float( obj[1][0][stat_index] )
        data = sorted( stats.iteritems(), key=sort_key )
        
        def iter_relevant_sorted_data():
            for i in data:
                yield i[0], i[1][0][stat_index], 3*i[1][1][stat_index]
            return
        
        relevant_sorted_data = list( iter_relevant_sorted_data() )
        return grp_sorted_data( relevant_sorted_data )

    num_cluster_grps = 0
    def iter_combined_clusters( cluster_mappings ):
        """Iterate trough groups which are the same in *every* cluster.

        """
        cluster_group_index = {}
        num_cluster_grps = 0
        for rg in stats.keys():
            cluster_group = []
            for cmap in cluster_mappings:
                cluster_group.append( cmap[ rg ] )
            cluster_group = tuple( cluster_group )
            # change the cluster groups ( currently in tuples )
            # to unique integers, in no particular order
            try:
                cluster_group = cluster_group_index[ cluster_group ]
            except KeyError:
                cluster_group_index[ cluster_group ] = num_cluster_grps
                cluster_group = num_cluster_grps
                num_cluster_grps += 1
            yield rg, cluster_group
    
    mean_cl = cluster_by_summary_stat( 0 )
    sd_cl = cluster_by_summary_stat( 1 )
    skew_cl = cluster_by_summary_stat( 2 )
    
    clusters = sorted( iter_combined_clusters( ( mean_cl, sd_cl, skew_cl ) ), \
                       key=operator.itemgetter(1) )
    
    if VERBOSE:
        print 'Read group clustering statistics'
        print 'read_group\tmean\t\t\tSD\t\t\tskew\t\t\tcluster_group'
        for read_group, cluster in clusters:
            # mean cluster will not have stats and will be skipped
            try:
                rg_stats = stats[ read_group ]
            except KeyError:
                print 'Should be mean and last cluster group:', \
                    read_group, str(cluster)
                continue
            print '%s\t%.5f+/-%.5f\t%.5f+/-%.5f\t%.5f+/-%.5f\t%s' % (read_group, \
                                 rg_stats[0][0], 0.10*rg_stats[0][0], \
                                 rg_stats[0][1], 0.10*rg_stats[0][1], \
                                 rg_stats[0][2], 4*rg_stats[1][2], \
                                 cluster )
    
    return dict( clusters )

def merge_clustered_fl_dists( clustered_read_groups, grouped_fragments ):
    """After we have determined the clustering, estimate the fl dists
       for the clustered read groups.

    Clustered read groups is a dictionary of cluster group names,
    keyed by read group. 
    
    Grouped fragments is a dictionary of lists of fragments, keyed by 
    read groups. 
    """
    # build a mapping from cluster groups to read groups
    cluster_grp_mapping = defaultdict(list)
    for rg, cg in clustered_read_groups.iteritems():
        cluster_grp_mapping[ cg ].append( rg )

    # initialize cluster groups mapping to fl_dists
    clustered_fl_dists = {}

    for cg, rgs in cluster_grp_mapping.iteritems():
        # initialize cluster group
        fls = numpy.array([])
        # combine the fragments from each read group
        for rg in rgs:
            fls = numpy.append( fls, grouped_fragments[ rg ] )
        # estimate the robust fl distribution
        clustered_fl_dists[ cg ] = build_robust_fl_dist_with_stats( fls )
    
    return clustered_fl_dists

def estimate_fl_dists( reads, genes, cluster_read_groups=True ):
    fragments = find_fragments( reads, genes )
    if len( fragments ) == 0:
        err_str = "There are no reads for this data file in the " \
            + "high read depth exons."
        raise ValueError, err_str

    # distributions of individual read_groups is not currently used
    # only clustered distributions are used for downstream analysis
    fl_dists = {}
    fl_dist_stats = {}
    grouped_fragments = group_fragments_by_readgroup( fragments )
    for read_group, fragment_lengths in grouped_fragments.iteritems():
        if len( fragment_lengths ) < MIN_FLS_FOR_FL_DIST:
            continue
        fl_dist = build_robust_fl_dist_with_stats( fragment_lengths )
        fl_dists[ read_group ] = fl_dist
        fl_dist_stats[ read_group ] = fl_dist.stats

    if cluster_read_groups:
        # no clustering occurs with large n using chained clustering
        # would need bootstrapping method to create valid clusters
        # get a dict of read_groups mapping to cluster group
        clustered_read_groups = cluster_rdgrps( fl_dist_stats )
        # get a dict of cluster groups mapping to fl_dists of merged read_groups
        fl_dists = merge_clustered_fl_dists( \
            clustered_read_groups, grouped_fragments)
    else:
        clustered_read_groups = dict( (rg, rg) for rg in grouped_fragments.keys() )
    
    return fl_dists, clustered_read_groups, fragments

def get_fl_dists( reads, single_exon_genes, \
                  cluster_read_groups=True, analyze_fl_dist=True ):
    """Wrap estimate_fl_dists, but include a caching layer. 
    
    """
    import os
    import pickle
    
    if VERBOSE:
        print "Loading FL Dists."
        
    reads_path = os.path.dirname( reads.filename )
    
    try:
        cached_obj_fname = get_fl_dist_base_from_fnames( \
                reads.filename, single_exon_genes.filename, cluster_read_groups ) + ".obj"
        
        # first, try to open the cached object in the reads directory
        try:
            with open( os.path.join( reads_path, cached_obj_fname ) ) as fp:
                fl_dists, read_group_mappings = pickle.load( fp )
        # if we can't try from the current directory
        except IOError:            
            with open( cached_obj_fname ) as fp:
                fl_dists, read_group_mappings = pickle.load( fp )
    except IOError:
        print "Estimating FL Dists. This may take some time..."
        fl_dists, read_group_mappings, fragments = \
            estimate_fl_dists( reads, single_exon_genes, cluster_read_groups )

        # first, try to open a file inb the reads directory. If this fails
        # ( ie, we dont have write proivleges or are out of space ) then 
        # fall back to the local directory.
        try:
            ofp = open( os.path.join( reads_path, cached_obj_fname ), "w" )
        except IOError:
            ofp = open( cached_obj_fname, "w" )
        
        pickle.dump( ( fl_dists, read_group_mappings ), ofp )
        ofp.close()
        
        if analyze_fl_dist:
            plot_filename = get_fl_dist_base_from_fnames( \
                reads.filename, single_exon_genes.filename ) + ".pdf"
            # try to put the plot in the reads directory
            try: 
                plot_full_fname = os.path.join( reads_path, plot_filename )
                
                with open( plot_full_fname, "w" ) as fp: pass
                analyze_fl_dists( fragments , plot_full_fname )
            # if we cant open a file there, then put the plot in the local directory
            except IOError:
                analyze_fl_dists( fragments, plot_filename )
        
    return fl_dists, read_group_mappings

def build_objs( gtf_fp, bam_fn ):
    import slide
    genes = slide.GeneBoundaries( gtf_fp )
    gtf_fp.close()

    # load the bam file
    reads = slide.Reads( bam_fn , "rb" )
    # make sure that it is indexed by trying to get a read
    try:
        reads.fetch("X", 0, 1)
    except ValueError:
        raise RuntimeError, "The bam file must be indexed."
    
    for cnt, read in enumerate( reads.fetch() ):
        if read.is_paired:
            break
        if cnt > 100:
            raise RuntimeError, "SLIDE is not compatible with single end reads"
        

    return genes, reads

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Estimate fragment length (fl) distributions by getting fragments from high read count single exon genes.' )
    parser.add_argument( 'gtf', type=file, \
                             help='GTF file for fl dist extimation (or generation)' )
    parser.add_argument( 'bam', help='bam file to get fl dist from')
    parser.add_argument( '--analyze', '-a', action="store_true", default=False, \
                          help='produce analysis graphs for fl dists' )
    parser.add_argument( '--generate', '-g', action="store_true", default=False, \
                          help='generate gtf file containing single exon genes with the highest depth of coverage' )
    parser.add_argument( '--out', '-o', default='high_read_depth_exons.gtf', \
                          help='output file for generation of gtf [default : %(default)s]' )
    parser.add_argument( '--repeat', '-r', \
                          help='optional file containing repeat masker bed file to remove generated genes' )
    args = parser.parse_args()

    return args.gtf, args.bam, args.analyze, args.generate, args.out, args.repeat

if __name__ == "__main__":
    gtf_fp, bam_fn, analyze, generate, out, repeat = parse_arguments()
    genes, reads = build_objs( gtf_fp, bam_fn )

    VERBOSE = True
    
    if generate:
        generate_gtf( reads, genes, out, repeat )
    else:
        get_fl_dists( reads, genes, cluster_read_groups=True, analyze_fl_dist=analyze )
