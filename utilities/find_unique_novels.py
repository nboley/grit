import sys

new_loci = set()
matched_loci = set()

for line in open( sys.argv[1] ):
    data = line.strip().split()
    if data[2] != '!': 
        matched_loci.add( int( data[3].split("_")[1] ) )
    else:
        new_loci.add( int( data[3].split("_")[1] ) )

novel_loci = sorted( new_loci - matched_loci )

for entry in novel_loci:
    print "cluster_%i" % entry
