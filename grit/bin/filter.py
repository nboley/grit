import os, sys
import re

def main():
    min_rpkm = float( sys.argv[1] )
    for line in open( sys.argv[2] ):
        if line.startswith( "track" ): continue
        frac = float(re.findall( '.*?frac "(.*?)";', line )[0])
        if frac >= min_rpkm:
            print line.strip()


if __name__ == '__main__':
    main()
