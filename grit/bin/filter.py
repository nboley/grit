import os, sys
import re

def parse_arguments():
    import argparse
    
    parser = argparse.ArgumentParser(
        description = 'Filter out GRIT transcripts.' )
    parser.add_argument(
        'gtf', type=file,
        help='GTF file to filter.' )
    parser.add_argument(
        '--FPKM', type=float, default=-1,
        help='Filter transceripts with an FPKM below this')
    parser.add_argument(
        '--conf-lo', type=float, default=-1,
        help='Filter transceripts with a conf_lo below this')
 
    args = parser.parse_args()   
    return args.gtf, args.FPKM, args.conf_lo

def main():
    gtf_fp, fpkm_thresh, conf_lo_thresh = parse_arguments()
    
    for line in gtf_fp:
        if line.startswith( "track" ): 
            print line,
            continue
        
        fpkm = float(re.findall( '.*?FPKM "(.*?)";', line )[0])
        if fpkm < fpkm_thresh: continue
        conf_lo = float(re.findall( '.*?conf_lo "(.*?)";', line )[0])
        if conf_lo < conf_lo_thresh: continue
        
        print line,

if __name__ == '__main__':
    main()
