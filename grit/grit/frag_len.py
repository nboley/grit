#!/usr/bin/python

# Copyright (c) 2011-2012 Nathan Boley

import numpy
from numpy import mean, std, append
import operator
import sys
from collections import defaultdict
import copy
import pickle
import random

import os

VERBOSE = False
MIN_FLS_FOR_FL_DIST = 100

# we stop collecting fragments after this many. We also make sure that we 
# have observed at least MAX_NUM_FRAGMENTS/10 exons
MAX_NUM_FRAGMENTS = 25000
MAX_NUM_FRAGMENTS_PER_EXON = 1000
MIN_EXON_LENGTH = 1000

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
            and (self.fl_density == other.fl_density).all():
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
    
    def __repr__( self ):
        return "<FL Dist: min=%.2f max=%.2f>"%(self.fl_min, self.fl_max)
    
def load_fl_dists( fnames ):
    all_fl_dists = {}
    all_read_group_mappings = []
    
    # load all of the individual fl dists
    for fname in fnames:
        with open( fname ) as fp:
            fl_dists = pickle.load( fp )
        for key, fl_dist in fl_dists.iteritems():
            all_fl_dists[ key ] = fl_dist
    
    # clustered_read_groups = cluster_rdgrps( all_fl_dists )
    clustered_read_groups = dict( (k,k) for k in all_fl_dists )
    
    return all_fl_dists, clustered_read_groups
    
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
    """Create a fragment length distribution from a mixture 
           of normal densities (i.e. simulated Illumina and 454)
    
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

def analyze_fl_dists( fragments, out_filename='diagnostic_plots.pdf'):
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
    for strand, strand_grouped_frag_lengths in \
            strand_grouped_frag_lengths.iteritems():
        strand_fl_dists[ strand ] = build_robust_fl_dist_with_stats( \
            strand_grouped_frag_lengths )
    
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

def find_fragments_in_exon( reads, exon ):
    fragment_sizes = []
    
    # calculate the length of the exon
    exon_len = exon.stop - exon.start + 1

    # iterate through read pairs in exon to get fragment lengths
    for cnt, (read1, read2) in enumerate(reads.iter_paired_reads( *exon )):
        # skip all secondary alignments
        if read1.is_secondary or read2.is_secondary:
            continue
        
        if len( read1.cigar ) > 1 or len( read2.cigar ) > 1:
            continue
        
        # get read group from tags
        # put all frags w/o a read group into mean read group
        try:
            read_group = [ val for key, val in read1.tags if key == 'RG' ][0]
        except IndexError:
            read_group = 'mean'
        read_group = 'mean'
        
        strand = '-' if read1.is_reverse else '+'

        frag_start = min( read1.pos, read2.pos, read1.aend, read2.aend )
        frag_stop = max( read1.pos, read2.pos, read1.aend, read2.aend )
        frag_len = frag_stop - frag_start + 1
        
        key = ( read_group, strand, exon_len )
        fragment_sizes.append( ( key, frag_len ) )
        
        if cnt > MAX_NUM_FRAGMENTS_PER_EXON:
            break
        
    return fragment_sizes

def find_fragments( reads, exons ):
    fragment_sizes = []
    for exon_i, exon in enumerate(exons):
        if (exon.stop - exon.start + 1) < MIN_EXON_LENGTH:
            continue
        
        fragment_sizes.extend( find_fragments_in_exon( reads, exon ) )
        if len( fragment_sizes ) > MAX_NUM_FRAGMENTS:
            break
    
    return fragment_sizes

def build_robust_fl_dist_with_stats( fragment_lengths ):
    """Trim outliers from a numpy array of fragment lengths.

    First, we estimate the mean and sd on the trimmed data. Then we use the 
    initial estimate to truncate at +/- NUM_SDS SD's
    
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
        print >> sys.stderr, "WARNING: There aren't any fragments after filtering."
        new_fls = sorted_fls
    
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

def cluster_rdgrps( fl_dists ):
    """

    stats is a dictioanry of readgroups, keyed by rg name, with values
    that take tuples of summary statistics.
    """
    stats = dict( (key, fl_dist.stats) for key, fl_dist in fl_dists.iteritems())

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

def estimate_fl_dists( reads, exons ):
    fragments = find_fragments( reads, exons )
    if len( fragments ) == 0:
        err_str = "There are no reads for this data file in the " \
            + "high read depth exons."
        raise ValueError, err_str
    
    # distributions of individual read_groups is not currently used
    # only clustered distributions are used for downstream analysis
    fl_dists = {}
    grouped_fragments = group_fragments_by_readgroup( fragments )
    for read_group, fragment_lengths in grouped_fragments.iteritems():
        if len( fragment_lengths ) < MIN_FLS_FOR_FL_DIST:
            continue
        fl_dist = build_robust_fl_dist_with_stats( fragment_lengths )
        fl_dists[ read_group ] = fl_dist
    
    return fl_dists, fragments

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description=\
        "Estimate fragment length (fl) distributions by getting fragments " +  \
        "from high read count single exon genes." )
    parser.add_argument( 'gff', type=file, \
        help='GFF file for fl dist extimation (or generation)' )
    parser.add_argument( 'bam', help='bam file to get fl dist from')
    parser.add_argument( 'outfname', help='output filename' )

    parser.add_argument( '--analyze', '-a', action="store_true", default=False,\
        help='produce analysis graphs for fl dists' )
    
    args = parser.parse_args()
    
    return args.gff, args.bam, args.outfname, args.analyze

def main():
    gff_fp, bam_fn, ofname, analyze = parse_arguments()
    
    from files.bed import parse_bed_line
    from files.reads import Reads
    
    exons = []
    for line in gff_fp:
        exon = parse_bed_line(line)
        if None == exon: continue
        exons.append( exon )
    
    # randomly shuffle the exons so that they aren't biased 
    random.shuffle(exons)
    gff_fp.close()
    
    # load the bam file
    reads = RNASeqReads( bam_fn , "rb" )
    
    VERBOSE = True
    fl_dists, fragments = estimate_fl_dists( reads, exons )
    
    if analyze:
        try:
            analyze_fl_dists( fragments, ofname + ".pdf" )
        except:
            pass
    
    with open( ofname, "w" ) as ofp:
        pickle.dump( fl_dists, ofp )
    
    return
    
    

if __name__ == "__main__":
    main()
