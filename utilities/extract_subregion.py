# Copyright (c) 2011-2012 Nathan Boley

import sys
import os
import subprocess
import pysam

EXTRACT_WIG_CMD = os.path.join( os.path.dirname( __file__ ), 
                                "extract_region_from_wiggle.py" )
EXTRACT_GFF_CMD = os.path.join( os.path.dirname( __file__ ), 
                                "extract_region_from_gff.py" )

class Region( object ):
    def __init__( self, chrm, start, stop ):
        self.chrm = chrm
        self.start = start
        self.stop = stop
        self.region_str = "%s_%i_%i" % ( chrm, start, stop )
        return
    
    def __str__( self ):
        return self.region_str
    
    def __repr__( self ):
        return self.region_str

def build_extract_bam_cmd( region, base_dir, 
                           sample_type, sample_id, fname, datatype=None):
    # try to detemrine if the chrm has a 'chr' prefix or not
    tmp_chr_name = None
    with pysam.Samfile( fname, "rb" ) as samfp:
        if samfp.gettid( "chr" + region.chrm ) != -1:
            tmp_chr_name = "chr" + region.chrm
        else:
            tmp_chr_name = region.chrm
            if samfp.gettid( region.chrm ) == -1:
                raise ValueError, "Neither chr%s nor %s exist within %s" \
                    % (region.chrm, region.chrm, fname)
    assert tmp_chr_name != None
    
    new_fname = os.path.join(base_dir, "%s_%s_%s.bam" % ( 
            sample_type, sample_id, region ))
    cmd1 = "samtools view -bh " + fname + " " \
         +  "%s:%i-%i" % (tmp_chr_name, region.start, region.stop) \
         +  " > " + new_fname
    cmd2 = "samtools index " + new_fname
    cmd = cmd1 + " && " + cmd2
    return cmd, new_fname

def build_extract_wig_cmd( region, base_dir, sample_type, sample_id, strand, 
                           fname, chrm_sizes_fname):
    """
    python extract_region_from_wiggle.py 
        chr4:+:1-500000 
        CAGE_AdMatedF_Ecl_1day_Heads_Trizol_Tissues.+.wig 
        /media/scratch/genomes/drosophila/dm3.chrom.sizes 
        tmp

    """
    sample_type_str = "merged" if sample_type == "*" else sample_type
    sample_id_str = "merged" if sample_id == "*" else sample_id
    
    new_fname_prefix = "%s_%s_%s" % (
        sample_type_str, sample_id_str, region )
    new_fname_prefix = os.path.join(base_dir, new_fname_prefix )
    
    region_str = "%s:%s:%i-%i" % ( region.chrm,strand,region.start,region.stop )
    
    cmd_template  = "python %s {0} {1} {2} {3}" % EXTRACT_WIG_CMD
    call = cmd_template.format( 
        region_str, fname, chrm_sizes_fname, new_fname_prefix )
    
    new_fname = new_fname_prefix + ".%s.wig" % strand
    return call, new_fname

def build_extract_g_f_cmd( region, base_dir, fname ):
    #new_fname = os.path.join( base_dir, ".".join(\
    #        os.path.basename(fname).split(".")[:-1]) \
    #        + ".%s_%i_%i.bam" % (  region.chrm, start, region.stop ) )
    new_fname = os.path.join(base_dir, os.path.basename(fname)
                             + str(region) + "." + fname.split(".")[-1] )

    region_str = "%s:%s:%i-%i" % ( region.chrm, '.', region.start, region.stop )
    
    cmd_template  = "python %s {0} {1} > {2}" % EXTRACT_GFF_CMD
    call = cmd_template.format( region_str, fname, new_fname )
    
    return call, new_fname

def get_filetype_from_datatype( datatype ):
    if datatype.lower().endswith( "bam" ):
        return 'bam'
    elif datatype.lower().endswith( "wig" ):
        return "wig"
    elif datatype.lower().endswith( "bedgraph" ):
        return "wig"
    elif datatype.lower().endswith( "gff" ) \
            or datatype.lower().endswith( "gtf" ):
        return "gff"
    else:
        return "UNKNOWN"
    assert False

def get_cmds_from_input_file( fp, region, base_dir ):
    chrm_sizes_fname = None
    new_lines = []
    cmds = []
    for line_num, line in enumerate(fp):
        line = line.strip()
        
        # skip commented out lines
        if line.startswith( "#" ):
            new_lines.append( line.strip() )
            continue
        
        if line == "": continue
        
        try:
            datatype, sample_type, sample_id, strand, fname = line.split()
        except:
            print line_num, line
            raise
        # deal with the chrm sizes specially
        if datatype == 'chrm_sizes':
            new_lines.append( line.strip() )
            chrm_sizes_fname = fname
            continue
        
        filetype = get_filetype_from_datatype( datatype )
        if filetype == 'UNKNOWN':
            new_lines.append( line.strip() )
            continue
        else:
            cmd, op_fname = None, None
            if filetype == 'bam':
                cmd, op_fname = build_extract_bam_cmd(
                    region, base_dir, sample_type, sample_id, fname, datatype)
            elif filetype == 'wig':
                cmd, op_fname = build_extract_wig_cmd(
                    region, base_dir, sample_type, sample_id, 
                    strand, fname, chrm_sizes_fname)
            elif filetype in ('gff', 'gtf'):
                cmd, op_fname = build_extract_g_f_cmd( region, base_dir, fname )
            else:
                raise ValueError, "Unrecognized line '%s'" % line.strip()
            
            cmds.append( cmd )
            new_lines.append( "\t".join(
                    (datatype, sample_type, sample_id, strand, op_fname)))

    ofp = open( os.path.join( base_dir, "test_elements.txt"), "w" )
    ofp.write( "\n".join( new_lines ) )
    ofp.close()
    
    ps = []
    for cmd in cmds:
        ps.append( subprocess.Popen( cmd, shell=True ) )
    
    for p in ps:
        print p
        p.wait()
    
    return

def parse_arguments():
    import argparse
    desc = 'Extract data files for a subregion from a grit command file.'
    parser = argparse.ArgumentParser(description=desc)    
    parser.add_argument(
        "--control-file", "-c", type=file, required=True,
        help='The GRIT control file containing the files to be parsed'  )

    parser.add_argument(
        "--region-str", "-r", type=str, required=True,
        help='The region that you want to extract ( ie chr4:100-10000 )' )
     
    parser.add_argument(
        '--output-directory', '-o', type=str, required=True,
        help='The directory to write the data files to.')
   
    args = parser.parse_args()
    
    chrm, poss = args.region_str.split(":")
    chrm = chrm[3:] if chrm.startswith( 'chr' ) else chrm
    start, stop = map( int, poss.split( "-" ) )
    
    return args.control_file, \
        Region( chrm, start, stop ), \
        os.path.abspath( args.output_directory )

if __name__ == "__main__":
    control_fp, region, base_dir = parse_arguments()
    get_cmds_from_input_file( control_fp, region, base_dir )    
