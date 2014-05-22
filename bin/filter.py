import os, sys
import re

pat = re.compile('transcript_id "(.*?)";')

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(
        description = 'Produce an expression gtf from a gtf and expression file.' )
    parser.add_argument(
        'gtf', type=file,
        help='GTF file to filter.' )
    parser.add_argument(
        'expression_tracking', type=file,
        help='Expression tracking file to build transcripts from.' )
    parser.add_argument(
        '--FPKM', type=float, default=-1,
        help='Filter transceripts with an FPKM below this')
    parser.add_argument(
        '--conf-lo', type=float, default=-1,
        help='Filter transcripts with a conf_lo below this')
    parser.add_argument(
        '--conf-hi', type=float, default=-1,
        help='Filter transcripts with a conf_hi below this')
 
    args = parser.parse_args()   
    return ( args.gtf, args.expression_tracking, 
             args.FPKM, args.conf_lo, args.conf_hi )

def load_expression_data(fp):
    rv = {}
    for line_num, line in enumerate(fp):
        # skip the header
        if line_num == 0: continue
        data = line.split()
        rv[data[0]] = data
    return rv

def main():
    ( gtf_fp, expression_fp, fpkm_thresh, conf_lo_thresh, conf_hi_thresh
      ) = parse_arguments()
    
    expression_data = load_expression_data(expression_fp)
    for line in gtf_fp:
        if line.startswith( "track" ): 
            print line,
            continue
        
        t_id = re.findall( pat, line)[0]
        data = expression_data[t_id]
        
        fpkm = data[3]
        if fpkm == '-' or float(fpkm) < fpkm_thresh: continue
        conf_lo = float(data[4])
        if conf_lo == '-' or conf_lo < conf_lo_thresh: continue
        conf_hi = float(data[5])
        if conf_hi == '-' or conf_hi < conf_hi_thresh: continue
        
        print line,

if __name__ == '__main__':
    main()
