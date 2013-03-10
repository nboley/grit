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
