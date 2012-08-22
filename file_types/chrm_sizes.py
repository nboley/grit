# Copyright (c) 2011-2012 Nathan Boley

import tempfile

class ChrmSizes( dict ):
    @staticmethod
    def _parse_genome_data_file_line( line ):
        chrm, size = line.split()
        if chrm.startswith( 'chr' ): 
            chrm = chrm[3:]
        return chrm, int( size )
    
    def build_file_str( self, append_chr=True ):
        res = []
        keys = sorted( self.keys() )
        for key in keys:
            res.append("%s\t%i" % ( key, self[key] ) )
        return "\n".join( res )
    
    
    def get_tmp_fname( self, with_chr=True ):
        return self.genome_data_fname
        
    def __init__( self, genome_data_fname ):
        # load the data file
        with open( genome_data_fname ) as genome_data_fp:
            for line in genome_data_fp:
                try:
                    chrm, size = self._parse_genome_data_file_line( line )
                except ValueError:
                    continue

                self[ chrm ] = size
        
        self.genome_data_fname = genome_data_fname
        
        return

def main():
    raise NotImplementedError

if __name__ == '__main__':
    main()

