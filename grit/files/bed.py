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

from collections import namedtuple
GenomicInterval = namedtuple('GenomicInterval', ['chr', 'strand', 'start', 'stop'])

from reads import clean_chr_name

def create_bed_line( chrm, strand, start, stop, 
                     name='.', score=1000, color='00,00,00',
                     use_thick_lines=True, blocks=[]):
    data = [ '%s' % chrm,
             "%s" % start,
             "%s" % stop,
             "%s" % name,
             "%s" % score,
             strand,
             "%s" % start,
             "%s" % stop,
             color  ]
    
    if len(blocks) > 0:
        n_blocks = len(blocks)
        b_starts, b_sizes = [], []
        for b_start, b_stop in blocks:
            assert start <= b_start <= b_stop  <= stop
            b_starts.append( b_start - start )
            b_sizes.append( b_stop - b_start + 1 )
    else:
        if use_thick_lines:
            n_blocks = 1
            b_starts = [0,]
            b_sizes = [stop-start,]
        else:
            n_blocks = 2
            b_starts = [0,stop-start-1]
            b_sizes = [1,1]
            
    data.append( "%i" % n_blocks )
    data.append( ",".join(map(str, b_sizes)) )
    data.append( ",".join(map(str, b_starts)) )
    
    return "\t".join( data )

def parse_bed_line(line):
    data = line.split()
    if len(data) < 6: return
    return GenomicInterval( clean_chr_name(data[0]), data[5], 
                            int(data[1]), int(data[2]) )
