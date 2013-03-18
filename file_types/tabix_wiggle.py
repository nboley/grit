import os, sys
import pysam
import numpy
import gzip
import subprocess

from chrm_sizes import ChrmSizes

VERBOSE = True

def clean_chr_name( chrm ):
    if chrm.startswith( "chr" ):
        chrm = chrm[3:]
    # convert the dmel chrm to M, to be consistent with ucsc
    if chrm == 'dmel_mitochondrion_genome':
        chrm = "M"
    return chrm

def guess_strand_from_fname( fname ):
    if fname.lower().rfind( "plus" ) >= 0:
        return '+'
    elif fname.lower().rfind( "+" ) >= 0:
        return '+'
    elif fname.lower().rfind( "minus" ) >= 0:
        return '-'
    elif fname.lower().rfind( "-" ) >= 0:
        return '-'
    else:
        raise ValueError, "Couldn't infer strand from filename '%s'" % fname
    
    assert False

class TabixBackedArray(object):
    def __init__( self, files, chrm, strand=None, contig_len=None ):
        self._data_files = [ pysam.Tabixfile( file ) 
                             if type(file) == str else file 
                             for file in files ]
        self.chrm = clean_chr_name(chrm)
        self.strand = strand
        self.contig_len = contig_len
        return
    
    def __getitem__( self, item ):
        start, stop = item.start, item.stop
        if start == None: start = 0 
        if stop == None: stop = self.contig_len 
        rv = numpy.zeros( stop - start  )
        for data_file in self._data_files:
            print data_file.filename
            for line in data_file.fetch( 
                    'chr' + self.chrm, start, stop, parser=pysam.asTuple() ):
                rv[max(0, int(line[1]))-start:int(line[2])-start] += float(line[3])
        
        return rv
    
    def asarray( self, contig_len=None ):
        if contig_len == None:
            contig_len = self.contig_len
        return self[0:contig_len+1]
    
    def __len__(self):
        return self.contig_len
    
def build_tabix_index( fp ):
    # check to see if we've been passed the compressed version. If
    # it ends with .gz, assume that we have
    fname = fp.name + ".gz" if not fp.name.endswith('.gz') else fp.name
    
    # check to see if the compressed version exists. if not, create it
    if not os.path.exists(fname):
        if VERBOSE: print >> sys.stderr, "Compressing ", fp.name        
        cmd = "bgzip -c %s > %s" % (fp.name, fname)
        subprocess.check_call( cmd, shell=True )
    
    # check to see if the tabix index exists. If not, create it
    if not os.path.exists(fname + '.tbi'):
        # check for a header line
        fp.seek(0)
        nskip = 1 if fp.readline().startswith("track") else 0
        fp.seek(0)
        
        if VERBOSE: print >> sys.stderr, "Indexing", fname
        cmd = "tabix %s -p bed -S %i" % ( fname, nskip )
        subprocess.check_call( cmd, shell=True )
    
    return fname

class Wiggle( dict ):
    def __init__( self, chrm_sizes_fp, fps, strands=None ):
        self.chrm_sizes = ChrmSizes( chrm_sizes_fp.name )
        
        for i, fp in enumerate(fps):
            # find the strand
            strand = strands[i] if strands != None \
                else guess_strand_from_fname(fp.name)
            
            # compress the file if necessary
            # index the file, if necessary
            fname = build_tabix_index( fp )
            
            for contig, contig_size in self.chrm_sizes.iteritems():
                if not self.has_key((contig, strand)):
                    self[(contig, strand)] = []
                self[(contig, strand)].append( fname )
        
        for (contig, strand), vals in self.iteritems():
            self[(contig, strand)] = TabixBackedArray( 
                self[(contig, strand)], contig, strand, self.chrm_sizes[contig])
        
        return

if __name__ == '__main__':
    chrm_lens_fp = open(sys.argv[1])
    fps = [ open( x ) for x in sys.argv[2:] ]
    x = Wiggle( chrm_lens_fp, fps )
    import wiggle
    y = wiggle.Wiggle( chrm_lens_fp, fps )
    for key in y.keys():        
        print "Diff indices: ", ( x[key].asarray() - y[key] ).nonzero()
    
