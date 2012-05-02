import sys, os
import shutil

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "../file_types/" ) )
from gtf_file import parse_gff_line

OP_DIR_NAME = "browser_track_data"

def write_filtered_gff( ifname, plus_ofname, minus_ofname, track_name ):
    ifp = open( ifname )
    plus_ofp = open( plus_ofname, "w" )
    minus_ofp = open( minus_ofname, "w" )
    
    plus_ofp.write("track name=plus_%s\n" % track_name )
    minus_ofp.write("track name=minus_%s\n" % track_name)
    for line in ifp:
        try:
            if line.split()[6] == '+':
                plus_ofp.write(line)
            else:
                minus_ofp.write(line)
        except IndexError:
            continue
    
    ifp.close()
    plus_ofp.close()
    minus_ofp.close()
    
    return

def copy_exons():
    # extract stranded exon versions, and copy 
    ifname = os.path.join( "./exons/", "discovered_exons.gff" )
    plus_ofname = os.path.join( OP_DIR_NAME, "plus_exons.gff" )
    minus_ofname = os.path.join( OP_DIR_NAME, "minus_exons.gff" )
    track_name = "disc_exons"
    write_filtered_gff( ifname, plus_ofname, minus_ofname, track_name )
    
    # copy the clustered exons
    shutil.copy( os.path.join("./exons/", "clustered_exons.gff"), OP_DIR_NAME )

    return

def copy_jns():
    ifname = os.path.join( "./junctions/", "merged.filtered.jns.gff" )
    plus_ofname = os.path.join( OP_DIR_NAME, "jns.plus.gff" )
    minus_ofname = os.path.join( OP_DIR_NAME, "jns.minus.gff" )
    track_name = "jns"
    write_filtered_gff( ifname, plus_ofname, minus_ofname, track_name )
    
    return

def copy_tsss():
    ifname = os.path.join( "./tss_exons/", "discovered_cage_tss_exons.gff" )
    plus_ofname = os.path.join( OP_DIR_NAME, "cage_exons.plus.gff" )
    minus_ofname = os.path.join( OP_DIR_NAME, "cage_exons.minus.gff" )
    track_name = "cage_exons"
    write_filtered_gff( ifname, plus_ofname, minus_ofname, track_name )
    
    return

def copy_tess():
    ifname = os.path.join( "./tes_exons/", "discovered_polya_tes_exons.gff" )
    plus_ofname = os.path.join( OP_DIR_NAME, "tes_exons.plus.gff" )
    minus_ofname = os.path.join( OP_DIR_NAME, "tes_exons.minus.gff" )
    track_name = "tes_exons"
    write_filtered_gff( ifname, plus_ofname, minus_ofname, track_name )
    
    return

def main():
    os.chdir( sys.argv[1] )
    try:
        os.mkdir( OP_DIR_NAME )
    except OSError:
        pass
    
    # assume the signal files are already there
    copy_exons()
    copy_jns()
    copy_tsss()
    copy_tess()

if __name__ == "__main__":
    main()
