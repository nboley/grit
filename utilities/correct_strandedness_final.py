import sys, os, copy, numpy, time, pickle
#from numpy import polyfit, poly1d
from scipy import stats, signal

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../file_types/fast_gtf_parser/" ) )
from gtf import load_gtf

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "./fast_wiggle_parser/" ) )

from wiggle import Wiggle

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../file_types/" ) )
from genomic_intervals import GenomicInterval
from gtf_file import iter_gff_lines


def list_samples( samp_fn ):
    '''
    obtain and organize the list of samples on which to fit the forest
    '''
    all_samples = {}
    fid = open(samp_fn)
    for fn in fid:
        short_name = '.'.join(fn.strip().split('/')[-1].split('.')[1:-2])
        very_short_name = '.'.join(short_name.split('.')[:-1])
        if not all_samples.has_key( very_short_name ):
            all_samples[very_short_name] = { short_name : [] }
        if not all_samples[very_short_name].has_key(short_name):
            all_samples[very_short_name][short_name] = []
        all_samples[very_short_name][short_name].append(fn.strip())
    for samp in all_samples.iterkeys():
        for rd in all_samples[samp].iterkeys():
            assert len(all_samples[samp][rd]) == 2
    return all_samples
        

def parse_wiggle(sample):
    '''
    sample should look like:
        { sample.rd1 : [+,-], sample.rd2 : [+,-] }

    returns a dict of Wiggle objects with exactly two entries, rd1 and rd2

    wiggles, out_fname_prefix, chrm_sizes_fp, track_name_prefix, filter_region \
        = parse_arguments()
    '''
    chrm_sizes_fp = open('/media/scratch/genomes/drosophila/dm3.chrom.sizes')
    wiggle_dict = {}
    for short_name in sample.keys():
        rd = short_name.split('.')[-1]
        wiggle_dict[rd] = Wiggle( chrm_sizes_fp )
        for fn in sample[short_name]:
            fid = open(fn)
            wiggle_dict[rd].load_data_from_fp( fid )
    return wiggle_dict



def fit_linear( RNAseq ):
    ''' This function conducts an L1 fit of a linear model of antisense.
    The vast majority of this code is concerned with finding the appropriate
    genomic coordinates to serve as input to this fitting.
    '''

    def get_range_nonsym( myp, mym, L, low_perc, up_perc ):
        ind_p = numpy.nonzero(myp)[0]
        ind_m = numpy.nonzero(mym)[0]

        lower_p = stats.scoreatpercentile(myp[ind_p],low_perc)
        upper_p = stats.scoreatpercentile(myp[ind_p],up_perc)
        lower_m = stats.scoreatpercentile(mym[ind_m],low_perc)
        upper_m = stats.scoreatpercentile(mym[ind_m],up_perc)

        ind_lower_p = set( numpy.nonzero( myp <= lower_p )[0] )
        ind_upper_p = set( numpy.nonzero( myp >= upper_p )[0] )
        ind_lower_m = set( numpy.nonzero( mym <= lower_m )[0] )
        ind_upper_m = set( numpy.nonzero( mym >= upper_m )[0] )

        myind_p = ind_lower_m.intersection(ind_upper_p)
        myind_m = ind_lower_p.intersection(ind_upper_m)

        myind_p = numpy.asarray( list( myind_p.intersection(ind_m) ), dtype=int )
        myind_m = numpy.asarray( list( myind_m.intersection(ind_p) ), dtype=int )
        
        return myind_p, myind_m


    def compute_obj_L1( top, bottom, alpha ):
        return abs(bottom - alpha*top).mean()

    ### I have found that 10% and 90% does a pretty good job
    # this seems to result in a noise-driven, rather than an overlapping-gene-
    # driven threshold.
    low_perc = 10 
    up_perc = 90  

    PLUS_top = []
    MINUS_bottom = []
    MINUS_top = []
    PLUS_bottom = []
    for (chrm,strand) in RNAseq['rd1'].iterkeys():
        if strand == '-': continue
        T1 = time.time()
        P1 = RNAseq['rd1'][(chrm,strand)] 
        M2 = RNAseq['rd2'][(chrm,'-')]
        L = len(P1)
        inds_p, inds_m = get_range_nonsym( P1, M2, L, low_perc, up_perc )
        try:
            print >>sys.stderr, len(inds_p), len(inds_m), max(P1[inds_m]), max(M2[inds_p]), min(P1[inds_p]), min(M2[inds_m])
        except:
            print >>sys.stderr, len(inds_p), len(inds_m)
        PLUS_top.extend(list(P1[inds_p]))
        PLUS_bottom.extend(list(P1[inds_m]))
        MINUS_top.extend(list(M2[inds_m]))
        MINUS_bottom.extend(list(M2[inds_p]))
        print >>sys.stderr, time.time()-T1, chrm 
    TOP = PLUS_top
    TOP.extend(MINUS_top)
    BOT = PLUS_bottom
    BOT.extend(MINUS_bottom)
    TOP = numpy.asarray(TOP)
    BOT = numpy.asarray(BOT)
    T1 = time.time()
    ans = []
    X = list( numpy.linspace(0,0.01,1000) )
    X.pop(0)
    for alpha in X:
        ans.append( compute_obj_L1( TOP, BOT, alpha ) )
    print >>sys.stderr, time.time()-T1
    ans = numpy.asarray(ans)
    #print ','.join( map(str,ans) ) 
    return X[ans.argmin()]


def main():
    # get all the RNA-seq wiggle file names, organize by sample and by read #
    all_samples = list_samples( sys.argv[1] )
    assert len(all_samples.keys()) == 1
    name = all_samples.keys()[0]
    # get the RNA-seq wiggle
    RNAseq = parse_wiggle(all_samples[name])
    alpha = fit_linear(RNAseq)
    print alpha 


main()
    



