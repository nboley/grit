import sys, os, copy, numpy, time, pickle
from sklearn.ensemble import RandomForestClassifier
from bx.intervals.intersection import Intersecter, Interval
from itertools import izip

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


################################################################################
# The file names that need to be loaded in:
genome_fname = '/media/scratch/DATA/genomes/hg19/ucsc_all.fa'
annotation_fname = '/home/ben/human_polyA/UCSC_known.gtf'
polyA_reads_fname = '/home/ben/human_polyA/epilepsy_polyA.gff'
cDNA_tes_fname = '/home/ben/human_polyA/polyA_DB.gtf'
wiggle_list = sys.argv[1] # A file like:
# /media/scratch/final_transcriptome_v2/read_cov_bedgraphs/rnaseq_cov.AdMatedF_Ecl_1day_Heads.minus.bedGraph
# /media/scratch/final_transcriptome_v2/read_cov_bedgraphs/rnaseq_cov.AdMatedF_Ecl_1day_Heads.plus.bedGraph
# /media/scratch/final_transcriptome_v2/read_cov_bedgraphs/rnaseq_cov.AdMatedF_Ecl_20days_Heads.minus.bedGraph
# /media/scratch/final_transcriptome_v2/read_cov_bedgraphs/rnaseq_cov.AdMatedF_Ecl_20days_Heads.plus.bedGraph
#
################################################################################


################################################################################


def reverse_strand( seq ):
    flip = {'a' : 't', 'c' : 'g', 'g' : 'c', 't' : 'a', 'n' : 'c'}
    return ''.join([flip[base] for base in seq[::-1]])


# The upstream and downstream motifs
#
# Retelska et al. BMC Genomics 2006 7:176   doi:10.1186/1471-2164-7-176
use = '''92.1 2.26 1.36 4.26
74.72 0.54 4.57 20.14
1.76 1.02 3.25 93.94
98.43 0.15 1.36 0.04
96.67 2.49 0.28 0.55
99.46 0.18 0.18 0.17'''.split('\n')
#import pdb; pdb.set_trace()
mRNA_LUSE = [ map(float,u.split(' ')) for u in use ]
LUSE = [ mRNA[::-1] for mRNA in mRNA_LUSE[::-1] ] 


# Retelska et al. BMC Genomics 2006 7:176   doi:10.1186/1471-2164-7-176 and meme
mRNA_list = ['ataaa', 'attaaa', 'agtaaa', 'tataaa', 'aataaa', 'aattaaa', 'aagtaaa', 'atataaa']
word_list = []
for word in mRNA_list:
    word_list.append( reverse_strand( word ) ) 



# T G/T G/T T/G G/T G/T C/T
# T T T T G T T 
# Retelska et al. BMC Genomics 2006 7:176   doi:10.1186/1471-2164-7-176
dse = '''8.72 6.62 10.52 74.13
1.72 18.64 37.31 42.3
4.94 20.65 9.25 65.15
1.52 68.43 14.13 15.89
8.66 0.15 0.00 91.16
0.11 7.63 59.4 32.85
9.08 20.42 22.58 47.9'''.split('\n')

mRNA_LDSE = [ map(float,d.split(' ')) for d in dse ]
LDSE = [ mRNA[::-1] for mRNA in mRNA_LDSE[::-1] ] 


meme_use_cDNA = '''0.718000 0.068000 0.046000 0.168000 0
0.732000 0.074000 0.040000 0.154000 0
0.000000 0.104000 0.020000 0.876000 0
0.762000 0.074000 0.158000 0.006000 0
0.710000 0.160000 0.030000 0.100000 0
0.848000 0.034000 0.104000 0.014000 0
0.262000 0.222000 0.216000 0.300000 0'''.split('\n')

mRNA_MUSE = [ map(float,u.split(' ')) for u in meme_use_cDNA ]
MUSE = [ mRNA[::-1] for mRNA in mRNA_MUSE[::-1] ] 
################################################################################



def list_samples( samp_fn ):
    '''
    obtain and organize the list of samples on which to fit the forest
    '''
    all_samples = {}
    fid = open(samp_fn)
    for fn in fid:
        short_name = '.'.join(fn.strip().split('/')[-1].split('.')[:-2])
        very_short_name = short_name.split('.')[0]
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
    chrm_sizes_fp = open('/media/scratch/DATA/genomes/hg19/hg19.chrom.sizes')
    wiggle_dict = {}
    for short_name in sample.keys():
        rd = short_name.split('.')[-1]
        wiggle_dict[rd] = Wiggle( chrm_sizes_fp )
        for fn in sample[short_name]:
            fid = open(fn)
            wiggle_dict[rd].load_data_from_fp( fid )
    return wiggle_dict

def parse_fasta( fn ):
    '''
    load a fasta file into a dictionary pointing to sinlge strings, one for
    each chromosome
    '''
    genome = dict()
    fid = open(fn)
    chrm = ''
    for line in fid:
        data = line.strip()
        if data.startswith('>'):
            chrm = data[1:]
        else:
            if not genome.has_key(chrm):
                genome[chrm] = []
                print >>sys.stderr, chrm
            genome[chrm].append(data.lower())
    for chrm in genome.keys():
        genome[chrm] = ''.join(genome[chrm])
    fid.close()
    return genome


def polyA_gff_2_dict( fn ):
    '''
    load a polyA gff file into a dictionary object
    chr3L	Read	CDS	17393446	17393502	.	-	.	@HWI-ST382_0049:1:1:4511:58676#CGATGT/1
    '''
    polyA = dict()
    fid = open(fn)
    for line in fid:
        data = line.strip().split('\t')
        chrm = data[0]  
        strand = data[6]
        if strand == '-':
            p_site = int(data[3])
        else:
            assert strand == '+'
            p_site = int(data[4])
        p_site -= 1 # subtract 1 so that this indexes into the fasta dict object
        if not polyA.has_key( (chrm,strand) ):
            polyA[(chrm,strand)] = dict()
        if not polyA[(chrm,strand)].has_key( p_site ):
            polyA[(chrm,strand)][ p_site ] = 0
        polyA[(chrm,strand)][ p_site ] += 1
    fid.close()
    return polyA

def polyA_dict_2_intersecter( polyA ):
    '''
    load a polyA gff file into a dictionary/intersecter object
    NOTE: bx.python is (open,open) in its interval searchers,
    e.g. T.find( 1,10 ) will not return true if T contains (1,1) or (10,10)
    '''
    polyA_I = dict()
    for (chrm, strand) in polyA.keys():
        if not polyA_I.has_key((chrm,strand)):
            polyA_I[(chrm,strand)] = Intersecter()
        for p_site in polyA[(chrm,strand)]:
            polyA_I[(chrm,strand)].add( p_site,p_site, [ p_site, polyA[(chrm,strand)][p_site] ] ) 
    return polyA_I
            

def get_elements_from_gene( gene, get_tss=True, get_jns=True, \
                                get_tes=True, get_exons=False ):
    tss_exons = set()
    tes_exons = set()
    introns = set()
    exons = set()
    
    chrm, strand = gene[1], gene[2]
    transcripts = gene[-1]
    
    for trans in transcripts:
        bndries = trans[1]

        fp_region = GenomicInterval(chrm, strand, bndries[0], bndries[1])
        tp_region = GenomicInterval(chrm, strand, bndries[-2], bndries[-1])
        if strand == '+':
            if get_tss:
                tss_exons.add( fp_region )
            if get_tes:
                tes_exons.add( tp_region )
        else:
            if strand != '-':
                print >> sys.stderr, "BADBADBAD", strand
                continue
            assert strand == '-'
            if get_tss:
                tss_exons.add( tp_region )
            if get_tes:
                tes_exons.add( fp_region )
        
        if get_jns:
            for start, stop in izip( bndries[1:-2:2], bndries[2:-1:2] ):
                # add and subtract 1 to ge tthe inclusive intron boundaries,
                # rather than the exon boundaries
                if start >= stop:
                    continue
                introns.add( GenomicInterval(chrm, strand, start+1, stop-1) )

        if get_exons:
            for start, stop in izip( bndries[::2], bndries[1::2] ):
                exons.add( GenomicInterval(chrm, strand, start, stop) )
    
    return tss_exons, introns, tes_exons, exons

def get_element_sets( genes, get_tss=True, get_jns=True, \
                          get_tes=True, get_exons=True ):
    tss_exons = set()
    introns = set()
    tes_exons = set()
    exons = set()
    for gene in genes:
        i_tss_exons, i_introns, i_tes_exons, i_exons = \
            get_elements_from_gene( gene, get_tss, get_jns, get_tes, get_exons )
        
        tss_exons.update( i_tss_exons )
        introns.update( i_introns  )
        tes_exons.update( i_tes_exons )
        exons.update( i_exons )

    return sorted(tss_exons), sorted(introns), sorted( exons ), sorted(tes_exons)


                      
def gtf_2_intersecters_and_dicts( gtf_fname ):
    '''
    parse a gtf file into two intersecters: CDSs and introns
    use the fast_gtf parser to get introns, and get CDSs via brute force
    (I realize this is incredibly stupid, but don't want to muck up intron boundaries)
    '''
    # get the Intron intersecter and interval objects
    def GenomicInterval_2_intersecter_and_dict( GI ):
        II = dict()
        ID = dict()
        for intron in GI:
            chrm = "chr" + intron.chr
            strand = intron.strand
            if not II.has_key( ( chrm, strand ) ):
                II[ (chrm,strand) ] = Intersecter()
            II[ (chrm,strand) ].add( intron.start,intron.stop, [ intron.start,intron.stop ] )
            if not ID.has_key( ( chrm, intron.strand ) ):
                ID[ (chrm,strand) ] = []
            ID[ (chrm,strand) ].append( [intron.start,intron.stop] )  
        return II, ID

    # get the CDS intersecters and interval objects
    def gtf_CDSs_2_intersecter_and_dict( gtf_fn ):
        fid = open(gtf_fname)
        CDS_I = dict()    
        CDS_D = dict()    
        for line in fid:
            data = line.strip().split('\t')
            if not data[2] == 'CDS':
                continue
            chrm = data[0]
            strand = data[6]
            start = int(data[3])
            end = int(data[4])
            if not CDS_I.has_key( (chrm,strand) ):
                CDS_I[ (chrm,strand) ] = Intersecter()
            CDS_I[ (chrm,strand) ].add( start,end, [start,end] )
            if not CDS_D.has_key( (chrm,strand) ):
                CDS_D[ (chrm,strand) ] = []
            CDS_D[ (chrm,strand) ].append( [start,end] )
        return CDS_I, CDS_D
            

    # load the genes and build sorted, unique lists
    genes = load_gtf( gtf_fname )
    tss_exons, introns, exons, tes_exons = get_element_sets( \
        genes, True, True, True, True )

    # generate all the intersecters and intervals for the annotation
    Introns_Sect, Introns_Dict = GenomicInterval_2_intersecter_and_dict( introns )
    Exons_Sect, Exons_Dict = GenomicInterval_2_intersecter_and_dict( exons )
    CDSs_Sect, CDSs_Dict = gtf_CDSs_2_intersecter_and_dict( gtf_fname ) 
    return Introns_Sect, Introns_Dict, Exons_Sect, Exons_Dict, CDSs_Sect, CDSs_Dict



def purify_introns( Introns_Dict, Exons_Sect ):
    '''
    select only introns that don't contain exons, these are likely to be "pure"
    introns that lack any real polyA sites.
    '''
    pure = dict()
    for (chrm, strand) in Introns_Dict.keys():
        for intron in Introns_Dict[(chrm, strand)]:
            if not Exons_Sect[(chrm, strand)].find( intron[0], intron[1] ):
                if not pure.has_key((chrm, strand)):
                    pure[(chrm, strand)] = Intersecter()
                pure[(chrm, strand)].add( intron[0], intron[1], intron )
    return pure


def get_overlapping_elements( tes_dict, elements_I, w ):
    '''
    Find all elements (tes's) overlapping another element type
    '''
    start = w
    end = w+1
    over = dict()
    for (chrm,strand) in tes_dict.keys():
        if not elements_I.has_key((chrm,strand)):
            print >>sys.stderr, "warning, element_intersecter does not contain the chrm: ", chrm
            continue
        for tes in tes_dict[(chrm,strand)].keys():
            H = elements_I[(chrm,strand)].find(tes-start,tes+end)
            if H:
                if not over.has_key( (chrm,strand) ):
                    over[ (chrm,strand) ] = dict()
                over[ (chrm,strand) ][tes] = copy.deepcopy( tes_dict[ (chrm,strand) ][tes] )
    return over


def remove_overlapping_elements( tes_dict, elements_I, w ):
    '''
    Remove all elements (tes's) overlapping another element type
    '''
    start = w
    end = w+1
    over = dict()
    for (chrm,strand) in tes_dict.keys():
        if not elements_I.has_key((chrm,strand)):
            print >>sys.stderr, "warning, element_intersecter does not contain the chrm: ", chrm
            continue
        for tes in tes_dict[(chrm,strand)].keys():
            H = elements_I[(chrm,strand)].find(tes-start,tes+end)
            if not H:
                if not over.has_key( (chrm,strand) ):
                    over[ (chrm,strand) ] = dict()
                over[ (chrm,strand) ][tes] = copy.deepcopy( tes_dict[ (chrm,strand) ][tes] )
    return over



def extract_genome_sequence( genome, tes_dict, w ):
    '''
    Return an array of sequences each of size 2*w + 1 
    '''
    seqs = []
    start = w
    end = w+1
    for (chrm,strand) in tes_dict.keys():
        if not genome.has_key(chrm):
            print >>sys.stderr, "warning, genome sequence does not contain the chrm: ", chrm
            continue
        for tes in tes_dict[(chrm,strand)].keys():
            seq = genome[chrm][tes-start:tes+end]
            if strand == "-":
                seqs.append([[chrm,strand,tes,tes_dict[(chrm,strand)][tes]], reverse_strand(seq)])
            else:
                assert strand == "+"
                seqs.append([[chrm,strand,tes,tes_dict[(chrm,strand)][tes]], seq])
    return seqs



def find_indexes_of_word( word, seq ):
    '''
    find each location where a given word occurs in a sequcence
    '''
    all_indexes = set()
    curr_index = seq.find(word)
    L = len(seq)
    while curr_index >= 0:
        all_indexes.add(curr_index)
        curr_index = seq.find(word, curr_index+1, L)
    return all_indexes

def seq_2_index( seq ):
    '''
    Turn a DNA sequence into indicies 0-4 for speedy motif searches
    '''
    ind = {'a' : 0, 'c' : 1, 'g' : 2, 't' : 3, 'n' : 4}
    sind = [ ind[s] for s in seq[1] ]
    return sind


def search_for_motif( seq_ind, motif ):
    '''
    Do a reasonably speedy motif search of a sequence, return the vector of 
    scores
    '''
    L = len(motif)
    ls = len(seq_ind)
    scores = []
    for i in xrange( 0, ls - L ):
        curr_score = 0
        for j,m in enumerate(motif):
            curr_score += m[seq_ind[i+j]]
        scores.append(curr_score)
    return numpy.asarray(scores)



def extract_covariates_from_seqs( seqs, w, polyA_density_curr, RNA_density, RNA_header ):
    '''
    All the heavy lifting is done here, a massive function to get all the covariates
    that turn out to be important.  
    '''

    # initilize point names:
    all_points = []


    # here is the massive header of covariates that will be generated in this function
    # note that RNA-seq and polyA (local density) covariates will be appended.
    header = ['name','read_count','triplet_ID-1','triplet_ID_center','triplet_ID+1', 'reads_within_10bp','reads_within_20bp','reads_within_50bp',
        'count_ATAAA_20_40', 'dist_ATAAA_20', 'dist_ATAAA_40', 'count_ATAAA_0_20', 'loc_ATAAA_0_20', 
        'count_ATTAAA_20_40', 'dist_ATTAAA_20', 'dist_ATTAAA_40', 'count_ATTAAA_0_20', 'loc_ATTAAA_0_20',  
        'count_AGTAAA_20_40', 'dist_AGTAAA_20', 'dist_AGTAAA_40', 'count_AGTAAA_0_20', 'loc_AGTAAA_0_20', 
        'count_TATAAA_20_40', 'dist_TATAAA_20', 'dist_TATAAA_40', 'count_TATAAA_0_20', 'loc_TATAAA_0_20',
        'count_AATAAA_20_40', 'dist_AATAAA_20', 'dist_AATAAA_40', 'count_AATAAA_0_20', 'loc_AATAAA_0_20',
        'count_AATTAAA_20_40', 'dist_AATTAAA_20', 'dist_AATTAAA_40', 'count_AATTAAA_0_20', 'loc_AATTAAA_0_20',
        'count_AAGTAAA_20_40', 'dist_AAGTAAA_20', 'dist_AAGTAAA_40', 'count_AAGTAAA_0_20', 'loc_AAGTAAA_0_20',
        'count_ATATAAA_20_40', 'dist_ATATAAA_20', 'dist_ATATAAA_40', 'count_ATATAAA_0_20', 'loc_ATATAAA_0_20',
        'word_total_USE_counts_20_40','word_total_USE_counts_0_20',
        'mx_score_U_20_40', 'loc_mx_U_20_40', 'mx_score_mU_20_40', 'loc_mx_mU_20_40', 
        'mx_score_D_55_80', 'loc_mx_D_55_80', 'mx_score_D_80_100', 'loc_mx_D_80_100', 
        'sum_D_55_80', 'sum_D_80_100']

    # extend the header to include the RNA-seq covariates.
    header.extend(RNA_header)

    # Initialize the "Big X", the set of covariates, the predictor matrix
    Big_X = []

    # turn all the sequences into the indices 0-4
    seqs_inds = []
    for seq in seqs:
        seqs_inds.append( seq_2_index( seq ) )

    delete_this = []
    for ind,seq in enumerate(seqs):
        Big_X.append([])
        if len(seq[1]) < 101:
            key_code = '_'.join(map(str,seq_index[:-1]))
            local_density = polyA_density_curr[key_code]
            delete_this.append(ind)
            print >>sys.stderr, w, seq, local_density
            continue
            
        seq_ind = seqs_inds[ind]
        seq_index = seq[0]
        chrm = seq_index[0]
        sequence_name = '_'.join(map(str,seq_index))
        all_points.append(sequence_name)

        ##### Add a covariate ##################################################
        # get the local read count #############################################
        Big_X[ind].extend( [seq_index[-1]] )
        
        seq = seq[1] # don't need the positional information any more.

        ##### Add a covariate ##################################################
        # encode the letter triplet at the polyA site itself ###################
        Big_X[ind].extend( seq_ind[40:60] ) ####################################

        ##### Add a covariate ##################################################
        # get the local density ################################################
        key_code = '_'.join(map(str,seq_index[:-1])) ###########################
        local_density = polyA_density_curr[key_code] ###########################
        Big_X[ind].extend( local_density )

        # search for all the words. NOTE: These are currently all upstream elements
        total_count = 0
        total_count_0_20 = 0
        word_cov = []
        for word in word_list:
            all_occur = numpy.asarray(list(find_indexes_of_word( word, seq[19:40] )))
            # select covariates from the word occurrence locations
            # number occur in 20-40, location nearest 20, location nearest 40
            number = len(all_occur)
            nearest_20 = 0
            nearest_40 = 0
            if number > 0:
                nearest_20 = min( all_occur )
                nearest_40 = max( all_occur )
            occur_first_20 = numpy.asarray(list(find_indexes_of_word( word, seq[0:20] )))
            number_first_20 = len( occur_first_20 )
            total_count_0_20 += number_first_20
            mx_20 = 0
            if number_first_20 > 0:
                mx_20 = max( occur_first_20 )
            total_count += number
            word_cov.append( [number, nearest_20, nearest_40, number_first_20, mx_20] )

            ##### Add a covariate ##############################################
            # Add the position and counts of word occurances ###################
            Big_X[ind].extend( [number, nearest_20, nearest_40, number_first_20, mx_20] )

        ##### Add a covariate ##################################################
        # Add the total counts of word occurences ##############################
        Big_X[ind].extend( [total_count,total_count_0_20] )

            
        # do all the motif searchs
        U_20_40 = search_for_motif( seq_ind[19:40], LUSE )
        mU_20_40 = search_for_motif( seq_ind[19:40], MUSE )
        D_55_80 = search_for_motif( seq_ind[54:80], LDSE )
        D_80_100 = search_for_motif( seq_ind[79:100], LDSE )

        # get the motif-based covariates
        mx_U_20_40 = U_20_40.max()
        loc_mx_U_20_40 = U_20_40.argmax()
        mx_mU_20_40 = mU_20_40.max()
        loc_mx_mU_20_40 = mU_20_40.argmax()
        mx_D_55_80 = D_55_80.max()
        loc_mx_D_55_80 = D_55_80.argmax()
        mx_D_80_100 = D_80_100.max()
        loc_mx_D_80_100 = D_80_100.argmax()
        D_sum_55_80 = D_55_80.sum() # similarity of nucleotide frequencies
        D_sum_80_100 = D_80_100.sum() # similarity of nucleotide frequencies


        ##### Add a covariate ##################################################
        # add all the motif-related covariates #################################
        Big_X[ind].extend( [mx_U_20_40, loc_mx_U_20_40, mx_mU_20_40, loc_mx_mU_20_40, mx_D_55_80, loc_mx_D_55_80, mx_D_80_100, loc_mx_D_80_100, D_sum_55_80, D_sum_80_100] )


        # get the RNA densities 
        if chrm.startswith('chr'):
            key_code = [seq_index[0][3:]]
            key_code.extend(seq_index[1:-1])
            key_code = '_'.join(map(str,key_code))
        try:
            local_RNA = RNA_density[key_code]
        except:
            import pdb; pdb.set_trace()

        ##### Add a covariate ##################################################
        # add all the RNA-seq related covariates ###############################
        Big_X[ind].extend( local_RNA )

        Big_X[ind] = numpy.asarray(Big_X[ind])

        #if not len(Big_X[ind]) == len(header)-1:
        #    import pdb; pdb.set_trace()
    if len(delete_this) > 0:
        del Big_X[numpy.asarray(delete_this)]
    return numpy.asarray(Big_X), header, all_points




def get_local_read_density(polyA_reads_D, polyA_reads_I):
    '''
    get the local polyA read density.
    '''
    seq_dict = dict()
    for (chrm,strand) in polyA_reads_D.keys():
        for pos in polyA_reads_D[(chrm,strand)].keys():
            # e.g. 'key' will end up looking like: chr2L_+_2030538
            key = '_'.join([chrm,strand,str(pos)])
            seq_dict[key] = [ len( polyA_reads_I[(chrm,strand)].find(pos-10,pos+11) ),
                len( polyA_reads_I[(chrm,strand)].find(pos-20,pos+21) ),
                len( polyA_reads_I[(chrm,strand)].find(pos-50,pos+51) ) ]
    return seq_dict


def get_RNAseq_densities( all_samples, polyA ):
    '''
    get the local RNA-seq read densities 
    '''
    dense = dict()
    header = []
    wiggles = dict()
    for sample in all_samples:
        header.extend( [ sample + '_up_10_rd1', 
                         sample + 'down_10_rd1', 
                         sample + '_up_50_rd1', 
                         sample + 'down_50_rd1', 
                         sample + '_up_100_rd1', 
                         sample + 'down_100_rd1', 
                         sample + '_up_down_rat_10_rd1',
                         sample + '_up_down_rat_50_rd1', 
                         sample + '_up_down_rat_100_rd1' ] )
        header.extend( [ sample + '_up_10_rd2', 
                         sample + 'down_10_rd2', 
                         sample + '_up_50_rd2', 
                         sample + 'down_50_rd2', 
                         sample + '_up_100_rd2', 
                         sample + 'down_100_rd2', 
                         sample + '_up_down_rat_10_rd2', 
                         sample + '_up_down_rat_50_rd2', 
                         sample + '_up_down_rat_100_rd2' ] )
        header.extend( [ sample + '_up_down_rat_10_rd1_rd2', 
                         sample + '_up_down_rat_50_rd1_rd2', 
                         sample + '_up_down_rat_100_rd1_rd2' ] )
    import pdb; pdb.set_trace()
    for sample in all_samples.iterkeys():
        t1 = time.time()
        print >>sys.stderr, "Now loading wiggle for: " + sample
        wiggle = parse_wiggle(all_samples[sample])
        print >>sys.stderr, "To load, it took : " + str( time.time()-t1 )
        t1 = time.time()
        for (chrm,strand) in polyA.keys():
            for pos in polyA[(chrm,strand)].keys():
                if chrm.startswith('chr'):
                    chrm = chrm[3:]             
                key = '_'.join([chrm,strand,str(pos)])
                if not dense.has_key(key):
                    dense[key] = []

                # now there is an rd1 and and rd2 entry in wiggles
                upstream_10_rd1 = wiggle['rd1'][(chrm,strand)][pos-10:pos].sum()
                downstream_10_rd1 = wiggle['rd1'][(chrm,strand)][pos:pos+10].sum()

                upstream_50_rd1 = wiggle['rd1'][(chrm,strand)][pos-50:pos].sum()
                downstream_50_rd1 = wiggle['rd1'][(chrm,strand)][pos:pos+50].sum()

                upstream_100_rd1 = wiggle['rd1'][(chrm,strand)][pos-100:pos].sum()
                downstream_100_rd1 = wiggle['rd1'][(chrm,strand)][pos:pos+100].sum()

                upstream_10_rd2 = wiggle['rd2'][(chrm,strand)][pos-10:pos].sum()
                downstream_10_rd2 = wiggle['rd2'][(chrm,strand)][pos:pos+10].sum()

                upstream_50_rd2 = wiggle['rd2'][(chrm,strand)][pos-50:pos].sum()
                downstream_50_rd2 = wiggle['rd2'][(chrm,strand)][pos:pos+50].sum()

                upstream_100_rd2 = wiggle['rd2'][(chrm,strand)][pos-100:pos].sum()
                downstream_100_rd2 = wiggle['rd2'][(chrm,strand)][pos:pos+100].sum()

                if strand == '+':
                    dense[key].extend([
                            upstream_10_rd1, downstream_10_rd1,
                            upstream_50_rd1, downstream_50_rd1,
                            upstream_100_rd1, downstream_100_rd1, 
                            upstream_10_rd1/max(downstream_10_rd1,1), 
                            upstream_50_rd1/max(downstream_50_rd1,1), 
                            upstream_100_rd1/max(downstream_100_rd1,1),
                            upstream_10_rd2, downstream_10_rd2,
                            upstream_50_rd2, downstream_50_rd2,
                            upstream_100_rd2, downstream_100_rd2, 
                            upstream_10_rd2/max(downstream_10_rd2,1), 
                            upstream_50_rd2/max(downstream_50_rd2,1), 
                            upstream_100_rd2/max(downstream_100_rd2,1),
                            upstream_10_rd1/max(downstream_10_rd2,1), 
                            upstream_50_rd1/max(downstream_50_rd2,1), 
                            upstream_100_rd1/max(downstream_100_rd2,1)])
                else:
                    dense[key].extend([
                            downstream_10_rd1, upstream_10_rd1,
                            downstream_50_rd1, upstream_50_rd1,
                            downstream_100_rd1, upstream_100_rd1, 
                            downstream_10_rd1/max(upstream_10_rd1,1), 
                            downstream_50_rd1/max(upstream_50_rd1,1),
                            downstream_100_rd1/max(upstream_100_rd1,1),
                            downstream_10_rd2, upstream_10_rd2,
                            downstream_50_rd2, upstream_50_rd2,
                            downstream_100_rd2, upstream_100_rd2, 
                            downstream_10_rd2/max(upstream_10_rd2,1), 
                            downstream_50_rd2/max(upstream_50_rd2,1), 
                            downstream_100_rd2/max(upstream_100_rd2,1),
                            downstream_10_rd1/max(upstream_10_rd2,1), 
                            downstream_50_rd1/max(upstream_50_rd2,1), 
                            downstream_100_rd1/max(upstream_100_rd2,1)])
        print >>sys.stderr, "To process, it took : " + str( time.time()-t1 ) 
    return dense, header


def print_bed_from_D( D ):
    for (chrm,strand) in D.keys():
        for pos in D[(chrm,strand)].keys():
            print '\t'.join( map(str,[chrm, pos, pos+1, strand, D[(chrm,strand)][pos]]) )


def print_fasta_from_seq( seq, out_fn, ind_start, ind_end ):
    fid = open(out_fn,'w')
    for line in seq:
        print >>fid, '>' + '_'.join(map(str,line[0]))
        print >>fid, line[1][ind_start:ind_end]
    return


def fit_forests( X_pos, X_neg_set, total_sets, size_train, size_test ):
    '''
    X_pos -- a numpy matrix.

    X_neg_set -- an array of numpy matrices corresponding the various negative
    control datasets that will be used to fit the forest, e.g. CDSs and Introns.
    
    total_sets -- the number of RFs to fit. An ensemble of ensembles is used 
    because the training data is not as larger or diverse as one might like.
    This is good for making sure that each classifier is not overfit while at 
    the same time utilizing all of the data for training.

    size_train -- a vector of training set sizes.  The first entry is the number
    of positive examples that will be drawn from X_pos, and remainder are for 
    the negs and should be in the same order as X_neg_set.

    size_test -- same as size_train but for the test set. 

    clf = RandomForestClassifier(n_estimators=10)
    sklearn.ensemble.RandomForestClassifier(n_estimators=10, 
    criterion='gini', max_depth=None, min_samples_split=1, min_samples_leaf=1, 
    min_density=0.1, max_features='auto', bootstrap=True, 
    compute_importances=False, 
    oob_score=False, n_jobs=1, random_state=None, verbose=0)
    '''

    # compute the sizes of the input datasets
    Lp = len(X_pos)
    Ln = []
    for X_neg in X_neg_set:
        Ln.append( len(X_neg) )

    # initilize the ensemble of ensembles
    Forests = []
    # store the error information
    Errs = []

    # build the labels for the training and test data
    train_labels = numpy.zeros(sum(size_train))
    train_labels[:size_train[0]] += 1
    test_labels = numpy.zeros(sum(size_test))
    test_labels[:size_test[0]] += 1
    for j in xrange( 0, total_sets ):
        t1 = time.time()
        # do the randomization to select the training and test sets
        pos_perm = numpy.random.permutation( Lp )
        neg_perm_set = []
        for L in Ln:
            neg_perm_set.append( numpy.random.permutation( L ) )

        # build the training data
        X_train = list( X_pos[ pos_perm[:size_train[0]] ] )
        for i,X_neg in enumerate(X_neg_set):
            X_neg_train = list( X_neg[ neg_perm_set[i][:size_train[i+1]] ] )
            X_train.extend(X_neg_train)
        X_train = numpy.asarray(X_train)

        # build the test data
        top = size_train[0]+size_test[0]
        X_test = list( X_pos[ pos_perm[size_train[0]:top] ] )
        for i,X_neg in enumerate(X_neg_set):
            top = size_test[i+1]+size_train[i+1]
            X_neg_test = list( X_neg[ neg_perm_set[i][size_train[i+1]:top] ] )
            X_test.extend(X_neg_test)
        X_test = numpy.asarray(X_test)

        # initilize and fit the forest
        Forests.append(RandomForestClassifier(n_estimators=100,max_features=80))
        Forests[j].fit( X_train, train_labels )
        
        # test the forest on the held-out test data
        test = Forests[j].predict( X_test )
        FN = sum(test[:size_test[0]]==0)/float(size_test[0])
        FP = sum(test[size_test[0]:]==1)/float(sum(size_test[1:]))
        Errs.append([FN,FP])
        print >>sys.stderr, FN, FP, time.time()-t1
    import pdb; pdb.set_trace()
    return Forests, Errs

   


def main():
    # load in the polyA reads
    polyA_reads_D = polyA_gff_2_dict( polyA_reads_fname )
    polyA_reads_I = polyA_dict_2_intersecter( polyA_reads_D )
    #import pdb; pdb.set_trace()
   

    # get all the RNA-seq wiggle file names, organize by sample and by read #
    all_samples = list_samples( wiggle_list )
    
    # set the size of the window we will extract
    window = 50
    
    # get the RNA-seq read densities
    RNA_dense, RNA_header = get_RNAseq_densities( all_samples, polyA_reads_D )

    # get local read density
    polyA_density = get_local_read_density(polyA_reads_D, polyA_reads_I)

    # load in the reference GTF
    Introns_Sect, Introns_Dict, Exons_Sect, Exons_Dict, CDSs_Sect, CDSs_Dict = (
        gtf_2_intersecters_and_dicts( annotation_fname )

    # load in the cDNA polyA ends
    cDNA_polyA_D = polyA_gff_2_dict( cDNA_tes_fname )
    cDNA_polyA_I = polyA_dict_2_intersecter( cDNA_polyA_D )
    #cDNA_density = get_local_read_density(cDNA_polyA_D, cDNA_polyA_I)

    # purify cDNAs to remove those that overlap CDSs:
    cDNA_polyA_noCDS_D = remove_overlapping_elements( 
            cDNA_polyA_D, CDSs_Sect, window )

    # get a set of "positive", polyA reads that we believe
    polyA_reads_cDNA_ends_D = get_overlapping_elements( 
            polyA_reads_D, cDNA_polyA_I, window )

    # purify "positives" to remove those that overlap CDSs:
    polyA_reads_cDNA_noCDS_D = remove_overlapping_elements( 
            polyA_reads_cDNA_ends_D, CDSs_Sect, window )

    # load in the reference genome indexed by chrm
    FA = parse_fasta( genome_fname )

    # extract genome sequences around cDNA polyA ends that don't overlap CDSs
    cDNA_polyA_noCDS_seqs = extract_genome_sequence( 
            FA, cDNA_polyA_noCDS_D, window )

    # extract genome sequences around cDNA polyA ends that don't overlap CDSs
    polyA_reads_cDNA_noCDS_seqs = extract_genome_sequence( 
            FA, polyA_reads_cDNA_noCDS_D, window )
    X_polyA_cDNA, header,point_names_polyA_cDNA = extract_covariates_from_seqs( 
            polyA_reads_cDNA_noCDS_seqs, 50, polyA_density,RNA_dense,RNA_header)

    # 1.a) get a set of introns that overlap no exons, TESs in these 
    # should be largely rubbish
    pure_introns_I = purify_introns( Introns_Dict, Exons_Sect )

    # 1.b) find the polyA reads that fall in these "pure" introns
    polyA_intronic_reads_D = get_overlapping_elements( 
            polyA_reads_D, pure_introns_I, 0 )
    
    # find the polyA reads that fall in CDSs
    polyA_CDS_reads_D = get_overlapping_elements( polyA_reads_D, CDSs_Sect, 0 )

    # extract sequences corresponding to negatives
    polyA_CDS_reads_seqs = extract_genome_sequence( 
            FA, polyA_CDS_reads_D, window )
    X_polyA_CDS, header, point_names_CDS = extract_covariates_from_seqs( 
            polyA_CDS_reads_seqs, 50, polyA_density, RNA_dense, RNA_header )
    polyA_intronic_reads_seqs = extract_genome_sequence( 
            FA, polyA_intronic_reads_D, window )
    X_polyA_intronic,header,point_names_intronic = extract_covariates_from_seqs(
            polyA_intronic_reads_seqs, 50, polyA_density, RNA_dense, RNA_header)

    import pdb; pdb.set_trace()

    # fit the forests:
    Forests, Errs = fit_forests( X_polyA_cDNA, [X_polyA_CDS, X_polyA_intronic], 
                                 3, [2000, 2000, 800], [2000, 2000, 800] )

    # get all polyA seqs
    polyA_reads_seqs = extract_genome_sequence( FA, polyA_reads_D, window )
    X_polyA_all, header, point_names_all = extract_covariates_from_seqs(
            polyA_reads_seqs, 100, polyA_density, RNA_dense, RNA_header )

    # do all the predictions for each forest
    preds = []
    L = len(Forests)
    fl = 1
    for forest in Forests:
        preds.append( forest.predict(X_polyA_all) )
    all_preds = numpy.zeros(len(preds[0]))
    for i in xrange(0,len(preds[0])):
        curr_pred = 0
        for j in xrange(0,L):
            curr_pred += preds[j][i]
        if curr_pred > fl:
            all_preds[i] = 1


    # collect all the polyA ends that pass prediction
    every_site = {}
    for i,p in enumerate(all_preds):
        if p == 1:
            every_site[point_names_all[i]] = X_polyA_all[i][0]
    # add on all the polyA ends from cDNAs
    #for (chrm,strand) in cDNA_polyA_noCDS_D.iterkeys():
    #    for pos in cDNA_polyA_noCDS_D[(chrm,strand)]:
    #        key_code = '_'.join([chrm,strand,str(pos)])
    #        if not every_site.has_key(key_code):
    #            every_site[key_code] = 10.5


    # print out bedGraphs of all clean polyA ends
    posfid = open('clean_454_polyA_sites_above_10.plus.bedGraph','w')
    minfid = open('clean_454_polyA_sites_above_10.minus.bedGraph','w')
    for key in every_site.iterkeys():
        data = key.split('_')
        chrm = data[0]
        strand = data[1]
        pos = data[2]
        score = str(every_site[key])
        if strand == '+':
            print >>posfid, '\t'.join( [chrm, pos, pos, score] )
        else:
            print >>minfid, '\t'.join( [chrm, pos, pos, score] )
        
    posfid.close()
    minfid.close()

    # Pickle Forest using protocol 0.
    #pkl_forest_fid = open('pickled_forest_antiCDS_thin.pkl', 'wb')
    #pickle.dump([Forests, Errs, header, all_samples, all_preds, every_site], pkl_forest_fid)
    #pkl_forest_fid.close()


    # Pickel the covariates in case we want to try to re-fit without having to bloody reload everything.
    #pkl_cov_fid = open('pickled_cov_antiCDS_thin.pkl', 'wb')
    #pickle.dump([X_polyA_cDNA, X_polyA_CDS, X_polyA_intronic, X_polyA_all, header, [5, [1000, 1500, 800], [1000, 1500, 800]]], pkl_cov_fid)
    #pkl_cov_fid.close()

    import pdb; pdb.set_trace()
    # gets 85% of FlyBase r5.45 3' ends, and maintains 13,865 intergenic polyA sites 

main()


