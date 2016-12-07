"""
Copyright (c) 2011-2015 Nathan Boley

This file is part of GRIT.

GRIT is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GRIT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GRIT.  If not, see <http://www.gnu.org/licenses/>.
"""

import os, sys

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), ".." ) )
from grit.files.reads import clean_chr_name, fix_chrm_name_for_ucsc
from grit.files.bed import GenomicInterval, parse_bed_line


def parse_arguments():
    allowed_assays = ['cage', 'rampage', 'rnaseq', 'polya']
    
    import argparse
    parser = argparse.ArgumentParser(
        description='Extract elements from a bed or gtf file.')
    parser.add_argument( 'region', 
                         help='A region string of the form contig:strand:start-stop')
    parser.add_argument( 'input_file', type=file,
                         help='An bed file to extract elements from')
    args = parser.parse_args()
    
    chrm, strand, poss = args.region.strip().split(":")
    start, stop = [int(x) for x in poss.replace(",", "").split('-')]
    
    return ( GenomicInterval(clean_chr_name(chrm), strand, start, stop), 
            args.input_file )
        
def build_track_line(track_line, region):
    data = track_line.split()
    name_field = [(i, x) for i, x in enumerate(data) if x.startswith("name")]
    assert len(name_field) <= 1
    # if there's no name field, return
    if len(name_field) == 0: return track_line

    region_str = "_region_" + "_".join(str(x) for x in region)
    name = "=".join(name_field[0][1].split("=")[1:])
    if name[-1] == '"':
        new_name = name[:-1] + region_str + '"'
    else:
        new_name = name + region_str
    return track_line.replace(name, new_name)

def main():
    region, fp = parse_arguments()
    for line in fp:
        if line.startswith('track'): 
            print(build_track_line(line, region), end=' ')
            continue
        element_region = parse_bed_line(line)
        if region.chr != element_region.chr: continue
        if region.strand != element_region.strand: continue
        if element_region.start > region.stop: continue
        if element_region.stop < region.start: continue
        print(line, end=' ')

if __name__ == '__main__':
    main()
