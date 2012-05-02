import sys
from collections import defaultdict

all_keys = defaultdict( int )
for line in open( sys.argv[1] ):
    data = line.split()
    if len( data ) == 0 or data[0] == 'track': continue
    num1 = data[11][1:-2]
    num2 = data[13][1:-1]
    all_keys[ ( num1, num2 ) ] += 1

num_dups = 0
for key, val in all_keys.iteritems():
    if val > 1:
        print key, val
        num_dups += 1

if num_dups == 0:
    print "PASSED: TRANCRIPT ID / EXON NUMBERS are unique"
else:
    print "FAILED: TRANCRIPT ID / EXON NUMBERS are NOT unique"
