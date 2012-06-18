import sys
from collections import defaultdict

genes = defaultdict( lambda: defaultdict( int ) )

for line in open( sys.argv[1]):
    data = line.split()
    code = data[2]
    gene_name = data[0]
    genes[gene_name][code] += 1

code_sums = defaultdict( int )
for gene_name, codes in genes.iteritems():
    if codes['='] > 0:
        code_sums['='] += 1
    elif codes['c'] > 0:
        code_sums['c'] += 1
    elif codes['j'] > 0:
        code_sums['j'] += 1
    elif codes['u'] > 0:
        code_sums['u'] += 1
    else:
        print codes
        assert False

print code_sums
print sum( code_sums.values() )
