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

from string import maketrans

intab = "ACGTNacgtn"
outtab = "TGCANtgcan"
seq_trans_table = maketrans(intab, outtab)

def reverse_comp_seq(seq):
    return seq[::-1].translate(seq_trans_table)

def iter_x_char_lines( data, x=80 ):
    pos = 0
    while pos < len(data):
        yield data[pos:pos+x]
        pos += x
    return
