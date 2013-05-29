import sys, os
import pysam
import numpy
import shutil
import subprocess
import tempfile
import time
from profilehooks import profile

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "grit/file_types/" ) )
from reads import iter_coverage_regions_for_read, clean_chr_name, \
    read_pairs_are_on_same_strand, get_strand


from multiprocessing import Process
num_threads = 1
buffer_size = 5000000
bedGraphToBigWig_script = "/usr/local/bin/bedGraphToBigWig"


def update_buffer_array_from_read( buffer_array, buffer_offset, 
                                   read, reads, reverse_read_strand, 
                                   pairs_are_opp_strand, chrm_len ):
    """populate buffer with histogram of contiguous read regions

    """
    for chrm, strand, start, stop in iter_coverage_regions_for_read(
            read, reads, reverse_read_strand, pairs_are_opp_strand):
        if stop >= chrm_len:
            return
        else:
            buffer_array[(start-buffer_offset):(stop-buffer_offset) + 1] += 1

def populate_cvg_arrays( reads_fname, chrm,
                         reverse_read_strand, pairs_are_opp_strand, 
                         chrm_len, strand):
                     
    """Get all of the junctions represented in a reads object
    """
    reads = pysam.Samfile( reads_fname, "rb" )

    # check if the chromosome exists in the references
    ref_index = reads.references.index( chrm )
    if ref_index == -1:
        reads.close()
        return

    # try to get one read. If we cant, then skip this chromosome
    try:
        next( iter(reads.fetch(chrm)) )
    except StopIteration:
        reads.close()
        return

    #create file for output
    if strand:
        ofp = file(os.path.abspath("/tmp/{0}.{1}.bedGraph".format(clean_chr_name(chrm), strand)), "w")
    else:
        ofp = file(os.path.abspath("/tmp/{0}.bedGraph".format(clean_chr_name(chrm))), "w")
    # we make this two times the buffer size so that we can do a memmove rather than
    # a new alloc when we write the buffer to disk. Note that the buffer size must be
    # larger than the gap between reads for this to work. For a normal RNA experiment,
    # 1MB should be more than enough, so we set it to 5MB to be extra safe.
    buffer_array = numpy.zeros( buffer_size*2 )
    # stores after how many bases in the current contig ( chrm, strand ) the current
    # buffer starts.
    buffer_offset = 0
    for i, read in enumerate(reads.fetch(chrm)):
        if strand and get_strand(read, reverse_read_strand, pairs_are_opp_strand) != strand:
            continue

        if i is 0:
            buffer_offset = read.pos

        # skip reads that aren't mapped in pair
        if not read.is_proper_pair:
            continue


        if read.pos > buffer_offset + buffer_size:
            write_array_to_opstream( ofp, buffer_array[:buffer_size], clean_chr_name(chrm), buffer_offset)
            # move the unwritten portion to the start of the buffer,
            # and zero out the end
            buffer_array[:buffer_size] = buffer_array[buffer_size:]
            buffer_array[buffer_size:] = 0
            buffer_offset += buffer_size

        update_buffer_array_from_read( buffer_array, buffer_offset, read, reads, reverse_read_strand, pairs_are_opp_strand, chrm_len )
    #to make sure the rest of the buffer is stored on disk
    write_array_to_opstream(ofp, buffer_array[:buffer_size], clean_chr_name(chrm), buffer_offset)
    reads.close()
    return



def write_array_to_opstream(ofp, buffer, chrm, buffer_offset ):
    """write buffer to disk, buffer offset determines the index of the read
    """
    lines = []
    prev_pos = 0
    prev_val = buffer[0]
    for pos, val in enumerate(buffer[1:]):
        if val != prev_val:
            if prev_val > 1e-12:
                lines.append( "chr%s\t%i\t%i\t%.2f" % (
                        chrm, buffer_offset + prev_pos, buffer_offset + pos+1, prev_val) )
            prev_pos = pos+1
            prev_val = val
    
    if prev_val > 1e-12:
        lines.append( "chr%s\t%i\t%i\t%.2f" % (chrm, buffer_offset+ prev_pos, buffer_offset + pos+2, prev_val) )
    if lines:
        ofp.write( "\n".join( lines ) + "\n" )
    return


def generate_wiggle(reads_fname, op_prefix, reverse_read_strand, reads_opp_strand, stranded):
    reads = pysam.Samfile( reads_fname, "rb" )
    print reads_opp_strand
    if reads_opp_strand == None:

        pairs_are_opp_strand = not read_pairs_are_on_same_strand( reads, num_reads_to_check =500000 )
    else:
        pairs_are_opp_strand = reads_opp_strand
    #not all bam files follow the header format if they dont assume they are sorted and indexed
    #if "HD" in reads.header and  "SO" in reads.header["HD"]:
        #print reads.header["HD"]["SO"]
        #assert reads.header["HD"]["SO"] != "unsorted", "The input should be sorted"

    if stranded:
        rd1_fps = {
        '+': file("{0}.{1}.bedGraph".format( op_prefix, "plus"), "w"),
        '-': file("{0}.{1}.bedGraph".format( op_prefix, "minus"), "w")
        }
    else:

        rd1_fps = file("{0}.bedGraph".format(op_prefix), "w")


    chrm_sizes_file = tempfile.NamedTemporaryFile(delete=False)
    chrm_names = reads.references
    chrm_lengths = zip(chrm_names, reads.lengths)
    #write out the chromosomes and its corrosponding size to disk
    for chrm, chrm_length in chrm_lengths:
        chrm_sizes_file.write(chrm + "   " + str(chrm_length) +"\n")

    reads.close()


    all_args = []
    if stranded:
        for chrm, chrm_length in chrm_lengths:
            all_args.append(( reads_fname, chrm, reverse_read_strand, pairs_are_opp_strand, chrm_length, '+'))
            all_args.append(( reads_fname, chrm, reverse_read_strand, pairs_are_opp_strand, chrm_length, '-'))
    else:
        for chrm, chrm_length in chrm_lengths:
            all_args.append(( reads_fname, chrm, reverse_read_strand, pairs_are_opp_strand, chrm_length,None))

    if num_threads == 1:
        for args in all_args:
            populate_cvg_arrays(*args)
    else:
        ps = []
        #create a list of processes to be run based on num_threads argument
        for thread_id, args in zip( range(num_threads), all_args ):
            p = Process(target=populate_cvg_arrays, name=args[1], args=args)
            p.start()
            ps.append( p )

        del all_args[:num_threads]


        while len(all_args) > 0:
            time.sleep( 0.5 )
            for p_i, p in enumerate(ps):
                if len( all_args ) == 0:
                    break
                if not p.is_alive():
                    ps[p_i].join()
                    args = all_args.pop()
                    ps[p_i] = Process(target=populate_cvg_arrays,
                                      name=args[1], args=args)
                    ps[p_i].start()

        for p in ps:
            p.join()

    #spawn processes to populate the coverage array

    if stranded:
        for chrm in chrm_names:
            temp = file(os.path.abspath("/tmp/{0}.{1}.bedGraph".format(clean_chr_name(chrm), "+")), "r")
            shutil.copyfileobj(temp, rd1_fps['+'] )
            temp = file(os.path.abspath("/tmp/{0}.{1}.bedGraph".format(clean_chr_name(chrm), "-")), "r")
            shutil.copyfileobj(temp, rd1_fps['-'])
        chrm_sizes_file.close()
        rd1_fps['+'].close()
        rd1_fps['-'].close()

       # call the bed graph conversion program

        res1 = subprocess.call([os.path.abspath(bedGraphToBigWig_script), 
                                os.path.abspath("{0}.{1}.bedGraph".format( op_prefix, "plus")), 
                                os.path.abspath(chrm_sizes_file.name), 
                                os.path.abspath("{0}.{1}.bw".format( op_prefix, "plus"))])
        res2 = subprocess.call([os.path.abspath(bedGraphToBigWig_script), 
                                os.path.abspath("{0}.{1}.bedGraph".format( op_prefix, "minus")), 
                                os.path.abspath(chrm_sizes_file.name), 
                                os.path.abspath("{0}.{1}.bw".format( op_prefix, "minus"))])
        
        #clean up
        chrm_sizes_file.delete
        if res1 == 0 and res2 == 0:
            for chrm in chrm_names:
                os.system("rm " + os.path.abspath("/tmp/{0}.{1}.bedGraph".format(clean_chr_name(chrm), "+")))
                os.system("rm " + os.path.abspath("/tmp/{0}.{1}.bedGraph".format(clean_chr_name(chrm), "-")))

            os.system("rm " + os.path.abspath("{0}.*.bedGraph".format(op_prefix)))
    else:
        for chrm in chrm_names:
            temp = file(os.path.abspath("/tmp/{0}.bedGraph".format(clean_chr_name(chrm))), "r")
            shutil.copyfileobj(temp, rd1_fps)
        rd1_fps.close()
        chrm_sizes_file.close()

        # call the bed graph conversion program
        res1 = subprocess.call([os.path.abspath(bedGraphToBigWig_script), 
                                os.path.abspath("{0}.bedGraph".format( op_prefix)), 
                                os.path.abspath(chrm_sizes_file.name), 
                                os.path.abspath("{0}.bw".format( op_prefix))])

        #clean up
        chrm_sizes_file.delete
        #ensure you created the bw before cleaning up
        if res1 == 0:
            for chrm in chrm_names:
                os.system("rm " + os.path.abspath("/tmp/{0}.bedGraph".format(clean_chr_name(chrm))))

            os.system("rm " + os.path.abspath("{0}.bedGraph".format(op_prefix)))




def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(
        description='Get coverage bedgraphs from aligned reads.')
    parser.add_argument( '--mapped-reads-fname', required=True,
                         help='BAM or SAM file(s) containing the mapped reads.')
    parser.add_argument( '--out-fname-prefix', '-o', required=True, 
                         help='Output file(s) will be bigWig')
    
    parser.add_argument( '--verbose', '-v', default=False, action='store_true', 
                         help='Whether or not to print status information.')
    parser.add_argument( '--threads', '-t', default=10, type=int,
                         help='The number of threads to run.')

    parser.add_argument( '--reverse-read-strand', '-r', type=bool,
                         help='Whether or not to reverse the strand of the read. default: infer from junction reads')
    parser.add_argument('--stranded', type=bool,
                        help="true if reads are stranded")
    parser.add_argument('--reads-opp-strand', type=bool, 
                        help='True if reads are opp strand for each BAM file')
    parser.add_argument('--buffer-size', '-b', default=5000000, type=int,
                        help='The amount of memory(in base pairs) to use before flushing to disk')
    parser.add_argument('--force', '-f', default=False, action='store_true', 
                        help='force recreate of bigwig even if it exists')

    args = parser.parse_args()
    global buffer_size
    buffer_size = args.buffer_size
    global num_threads
    num_threads = args.threads

    if not args.reads_opp_strand:
        args.reads_opp_strand = None

    global VERBOSE
    VERBOSE = args.verbose

    global FORCE
    FORCE = args.force

    return args.mapped_reads_fname, args.out_fname_prefix, \
        args.reverse_read_strand, args.reads_opp_strand, not(args.not_stranded)

#python bam2wig.py --mapped-reads-fname ./50Mb/I008_50Mb_NA10831.bam -o /var/www/JBrowse-1.8.1/data/tracks/I008_50Mb_NA10831 -b 5000000 -t 20
def main():
    reads_fnames, op_prefixes, reverse_read_strand, reads_opp_strand, stranded  = parse_arguments()
    print stranded
    #BW already created skip
    if not FORCE and stranded and os.path.exists(os.path.abspath(op_prefixes + ".plus.bw")) and os.path.exists(os.path.abspath(op_prefixes + ".minus.bw")):
        print "skipped " + op_prefixes
        return
    if not FORCE and not stranded and os.path.exists(os.path.abspath(op_prefixes + ".bw")):
        print "skipped " + op_prefixes
        return
    generate_wiggle(reads_fnames, op_prefixes, reverse_read_strand, reads_opp_strand, stranded)



if __name__ == "__main__":
    main()

