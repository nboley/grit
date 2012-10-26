import os, sys


def main():
    ofname = ".".join(sys.argv[1].split(".")[:-1]) + ".bedgraph"
    print "writing to ", ofname
    ofp = open( ofname, "w" )
    with open( sys.argv[1] ) as fp:
        curr_strand = None
        curr_chrm = None
        curr_track = None
        for line in fp:
            # get the track name and strand
            if line.startswith("track"):
                curr_track = line.strip().split("=")[-1]
                if curr_track.endswith("plus") or curr_track.endswith("+"):
                    curr_strand = '+' 
                elif curr_track.endswith("minus") or curr_track.endswith("-"):
                    curr_strand = '-' 
                else:
                    raise ValueError, "Unrecognized strand for track name '%s'"\
                        % curr_track
                
                curr_strand = None
                curr_chrm = None
                
                print >> ofp, "track type=bedGraph name=%s" % curr_track
                continue
            
            if line.startswith( "variableStep" ):
                curr_chrm = line.strip().split("=")[-1]
                continue
            
            data = line.split()
            pos, value = int(data[0]), float(data[1])
            
            if value >= 1:
                print >> ofp, "%s\t%s\t%i\t%e" % ( 
                    curr_chrm, pos, pos+1, value )
    ofp.close()
    return

if __name__ == '__main__':
    main()
