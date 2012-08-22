import sys

def main():
    assert len( sys.argv ) == 2

    num_tracks = 0
    
    curr_chr = None
    prev_pos = None
    prev_val = None
    length = 0
    
    fname = sys.argv[1]
    fp = open( sys.argv[1] )

    assert fname.endswith( ".wig" )
    ofp = open( fname[:-4] + ".bedgraph", "w" )
    
    for line in fp:
        if line.startswith( "track" ):
            assert line.index("type=wiggle_0" ) > 0
            num_tracks += 1
            assert num_tracks <= 1
            continue
        
        if line.startswith( "variableStep" ):
            curr_chr = line.split("chrom=")[-1].split()[0]
            prev_pos = None
            prev_val = None
            length = 0
            print "New Chr: ", curr_chr
            continue
        
        pos, value = map( float, line.split() )
        if prev_pos == None:
            assert prev_val == None
            prev_pos = int( pos )
            prev_val = value
            length = 1
            continue
        elif pos == prev_pos + 1 and value == prev_val:
            length += 1
        else:
            ofp.write( "%s\t%i\t%i\t%e\n" % (curr_chr, prev_pos, prev_pos+length, prev_val) )
            prev_pos = int( pos )
            prev_val = value
            length = 1

    if prev_pos != None:
        ofp.write( "%s\t%i\t%i\t%e\n" % (curr_chr, prev_pos, prev_pos+length, prev_val) )

    ofp.close()
    fp.close()

main()
