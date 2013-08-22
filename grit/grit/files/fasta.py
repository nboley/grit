def iter_x_char_lines( data, x=80 ):
    pos = 0
    while pos < len(data):
        yield data[pos:pos+x]
        pos += x
    return
