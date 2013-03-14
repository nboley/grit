from collections import namedtuple
GenomicInterval = namedtuple('GenomicInterval', ['chr', 'strand', 'start', 'stop'])

from reads import clean_chr_name

def create_bed_line( chrm, strand, start, stop, 
                     name='.', score=1000, color='00,00,00',
                     use_thick_lines=True):
    return "\t".join( [
                'chr%s' % chrm,
                "%s" % start,
                "%s" % stop,
                "%s" % name,
                "%s" % score,
                strand,
                "%s" % start,
                "%s" % (stop if use_thick_lines else start),
                color
            ])

def parse_bed_line(line):
    data = line.split()
    if len(data) < 6: return
    if data[3] != 'single_exon_gene': return
    return GenomicInterval( clean_chr_name(data[0]), data[5], 
                            int(data[1]), int(data[2]) )
