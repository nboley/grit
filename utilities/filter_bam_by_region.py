import sys
import os
import subprocess
"""
example command:
python filter_bam_by_region.py `ls /media/scratch/RNAseq/all_samples/ | grep -P '^Ad|^L3_|^WPP_'`
"""

region_chr = "2L"
start = 4605421
stop = 5481770
base_dir = "/media/scratch/chr2L_test_region_raw_data/unstranded_timepoints/"

def build_extract_bam_cmd( fname ):
    new_fname = os.path.join( base_dir, ".".join(\
            os.path.basename(fname).split(".")[:-1]) \
            + ".%s_%i_%i.bam" % (  region_chr, start, stop ) )
    cmd1 = "samtools view -bh " + old_fname + " " \
         +  "%s:%i-%i" % (region_chr, start, stop) +  " > " + new_fname
    cmd2 = "samtools index " + new_fname
    cmd = cmd1 + " && " + cmd2
    return cmd

ps = []
for fname in sys.argv[1:]:
    cmd = build_extract_bam_cmd( fname )
    ps.append( subprocess.Popen( cmd, shell=True ) )

for p in ps:
    print p
    p.wait()
