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
    
    def _build_temp_files( self ):
        genome_lns_fp = tempfile.NamedTemporaryFile( 
            "w", prefix="genome_len_wo_chr", suffix='.tmp.txt' )
        genome_lns_fp.write( self.build_file_str( False ) )
        genome_lns_fp.flush()
        self.tmp_w_chrm_fp = genome_lns_fp
        
        genome_lns_fp_2 = tempfile.NamedTemporaryFile( 
            "w", prefix="genome_len_with_chr", suffix='.tmp.txt' )
        genome_lns_fp_2.write( self.build_file_str( True ) )
        genome_lns_fp_2.flush()
        self.tmp_wo_chrm_fp = genome_lns_fp_2
        
        return
    
    def get_tmp_fname( self, with_chr=True ):
        if with_chr:
            return self.tmp_w_chrm_fp.name
        return self.tmp_wo_chrm_fp.name
    
    def __del__( self ):
        try:
            self.tmp_w_chrm_fp.close()
        except OSError:
            pass
        
        try:
            self.tmp_wo_chrm_fp.close()
        except OSError:
            pass
    
    def __init__( self, genome_data_fp ):
        # load the data file
        for line in genome_data_fp:
            try:
                chrm, size = self._parse_genome_data_file_line( line )
            except ValueError:
                continue
            
            self[ chrm ] = size
        
        # build temporary files to pass as arguments to ecternal scripts
        self._build_temp_files()
        
        return

def main():
    raise NotImplementedError

if __name__ == '__main__':
    main()

