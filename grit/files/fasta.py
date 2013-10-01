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
