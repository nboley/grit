import os,  sys

sys.path.append( os.path.join(os.path.dirname(__file__), "..", "file_types") )
from gtf_file import parse_gff_line, GenomicInterval

def get_region( region_str ):
    chrm, strand, locs = region_str.strip().split(':')
    if chrm.startswith("chr"): chrm = chrm[3:]
    locs.replace( ",", "" )
    start, stop = map( int, locs.split("-") )
    return GenomicInterval( chrm, strand, start, stop )

def main():
    filter_region = get_region( sys.argv[1] )
    with open( sys.argv[2]) as fp:
        for line in fp:
            data = parse_gff_line( line )
            if data == None: continue
            region = data[0]
            if region.chr != filter_region.chr: continue
            if filter_region.strand != '.' \
                    and filter_region.strand != region.strand: 
                continue
            if region.stop < filter_region.start \
                    or region.start > filter_region.stop: 
                continue
            print line.strip()
    
    return

if __name__ == '__main__':
    main()
