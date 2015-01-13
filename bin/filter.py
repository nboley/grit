import os, sys
import re

from collections import defaultdict, namedtuple

pat = re.compile('transcript_id "(.*?)";')
pat_gene = re.compile('gene_id "(.*?)";')

ExpressionTrackingLine = namedtuple( 
    'ExpressionTrackingLine', 
    ["tracking_id", "gene_id", 
     "coverage", "FPKM", "FPKM_lo", "FPKM_hi", "status"])


def load_expression_tracking_data(fp):
    rv = {}
    for line_num, line in enumerate(fp):
        # skip the header
        if line_num == 0: continue
        data = line.split()
        for i, val in enumerate(data[2:6]):
            if val == '-': val = None
            else: val = float(val)
            data[i+2]  = val
        data = ExpressionTrackingLine(*data)
        rv[data.tracking_id] = data
    
    return rv

def find_transcript_stats(expression_fps):
    rv = defaultdict(lambda: defaultdict(list))
    for fp in expression_fps:
        for data in load_expression_tracking_data(fp).itervalues():
            rv[data.gene_id][data.tracking_id].append(data)
    
    return rv

def find_valid_transcripts_in_gene( gene_expression_data,
                                    min_fpkm_lb, min_fpkm_ub,
                                    intrasample_max_fpkm_ratio,
                                    intersample_max_fpkm_ratio ):
    max_fpkm_lb_across_samples = -1.0
    max_fpkm_lb_in_sample = defaultdict(lambda: -1.0)
    for t_id, t_data in gene_expression_data.iteritems():
        for i, sample_data in enumerate(t_data):
            max_fpkm_lb_across_samples = max(
                max_fpkm_lb_across_samples, sample_data.FPKM_lo )
            max_fpkm_lb_in_sample[i] = max(
                max_fpkm_lb_in_sample[i], sample_data.FPKM_lo )

    good_transcripts = set()
    for t_id, t_data in gene_expression_data.iteritems():
        # skip transcripts with lower bounds all below the threshold
        if all(t.FPKM_lo < min_fpkm_lb for t in t_data): continue
        
        for i, sample_data in enumerate(t_data):
            min_max_fpkm = max(
                min_fpkm_ub, 
                max_fpkm_lb_in_sample[i]/intrasample_max_fpkm_ratio,
                max_fpkm_lb_across_samples/intersample_max_fpkm_ratio)
            if sample_data.FPKM_hi >= min_max_fpkm: 
                good_transcripts.add(t_id)
    return good_transcripts

def find_valid_transcripts( expression_data,
                            min_fpkm_lb, min_fpkm_ub,
                            intrasample_max_fpkm_ratio,
                            intersample_max_fpkm_ratio ):
    valid_transcripts = set()
    for gene_id, transcript_data in expression_data.iteritems():
        valid_transcripts.update(
            find_valid_transcripts_in_gene(
                transcript_data,
                min_fpkm_lb, min_fpkm_ub,
                intrasample_max_fpkm_ratio,
                intersample_max_fpkm_ratio))
    return valid_transcripts

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(
        description = 'Produce an expression gtf from a gtf and expression file.' )
    parser.add_argument(
        'gtf', type=file,
        help='GTF file to filter.' )
    parser.add_argument(
        'expression_tracking_files', type=file, nargs='+',
        help='Expression tracking files to build transcripts from.' )
    parser.add_argument(
        '--min-fpkm-ub', default=0.0, type=float,
        help='Filter transcripts with fpkms upper bounds below this value.')
    parser.add_argument(
        '--min-fpkm-lb', default=0.0, type=float,
        help='Filter transcripts with fpkms upper bounds below this value.')
    parser.add_argument(
        '--intrasample-max-fpkm-ratio', type=float, default=1e6,
        help='For each gene cluster and sample, filter transcripts whose fpkm upper bound is this many times lower than the maximum fpkm lower bound.')
    parser.add_argument(
        '--intersample-max-fpkm-ratio', type=float, default=1e6,
        help='For each gene cluster and between all samples, filter transcripts whose fpkm upper bound is this many times lower than the maximum fpkm lower bound.')

 
    args = parser.parse_args()   
    return ( args.gtf, args.expression_tracking_files, 
             args.min_fpkm_lb, args.min_fpkm_ub,
             args.intrasample_max_fpkm_ratio, 
             args.intersample_max_fpkm_ratio )

def build_track_line(track_line, 
                     (min_fpkm_lb, min_fpkm_ub,
                      intrasample_max_fpkm_ratio, 
                      intersample_max_fpkm_ratio)):
    data = track_line.split()
    name_field = [(i, x) for i, x in enumerate(data) if x.startswith("name")]
    assert len(name_field) <= 1
    # if there's no name field, return
    if len(name_field) == 0: return track_line

    filter_str = ".filtered"
    name = "=".join(name_field[0][1].split("=")[1:])
    if name[-1] == '"':
        new_name = name[:-1] + filter_str + '"'
    else:
        new_name = name + filter_str
    return track_line.replace(name, new_name).strip()

def main():
    ( gtf_fp, expression_fps, 
      min_fpkm_lb, min_fpkm_ub, 
      intrasample_max_fpkm_ration,
      intersample_max_fpkm_ratio
      ) = parse_arguments()

    transcript_stats = find_transcript_stats(expression_fps)
    
    valid_transcripts = find_valid_transcripts(
        transcript_stats,       
        min_fpkm_lb, min_fpkm_ub, 
        intrasample_max_fpkm_ration,
        intersample_max_fpkm_ratio )
    
    for line in gtf_fp:
        if line.startswith( "track" ): 
            print build_track_line(
                line, (
                    min_fpkm_lb, min_fpkm_ub, 
                    intrasample_max_fpkm_ration,
                    intersample_max_fpkm_ratio ) )
            continue
        
        t_id = re.findall( pat, line)[0]
        if t_id not in valid_transcripts: continue
        gene_id = re.findall( pat_gene, line)[0]
        fpkm = transcript_stats[gene_id][t_id][0].FPKM
        print line.strip() + ' FPKM "%s";' % fpkm

if __name__ == '__main__':
    main()
