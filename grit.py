"""TODO

Currently, if you add a sample the merged samples won't be re-run.
I should check to make sure there are no new elements for the merged
elements, and delete them if there are. 


Also, merged trasncripts now builds off of the merged samples, but 
this can be a little inefficient becaue the merge command has to wait on 
other merge commands. Cahnge this to merge on the individual samples.

build_transcripts
merge_transcripts

"""

import sys, os
import tempfile
from collections import namedtuple, defaultdict
from subprocess import Popen
import time
import random

import sqlite3

# import the fast chksum function
from zlib import crc32
import itertools
MAX_NUM_CKSUM_BLOCKS = 10
USE_TF_ELEMENTS = False

USE_MERGED_EXONS_FOR_TSS_EXONS = False
BUILD_SAMPLE_SPECIFIC_EXONS = True
USE_SAMPLE_SPECIFIC_CAGE = True and BUILD_SAMPLE_SPECIFIC_EXONS
USE_MERGED_JNS_FOR_EXONS = False
USE_MERGED_EXONS_FOR_TRANSCRIPTS = False

VERBOSE = True
STRESS_TEST_DEP_FINDER = False


EXTRACT_JNS_CMD = os.path.join( os.path.dirname( __file__ ), 
                              "./elements/junctions/", "extract_junctions.py" )
MERGE_AND_FILTER_JNS_CMD = os.path.join( os.path.dirname( __file__ ), 
                     "./elements/junctions/", "merge_and_filter_junctions.py" )
INTERSECT_JNS_CMD = os.path.join( os.path.dirname( __file__ ), 
                            "./elements/junctions/", "intersect_junctions.py" )

BUILD_READ_CVG_CMD = os.path.join( os.path.dirname( __file__ ), 
                             "./elements/", "build_read_coverage_bedgraph.py" )

FIND_EXONS_CMD = os.path.join( os.path.dirname( __file__ ), 
                                "./elements/exons/", "discover_exons.py" )

MERGE_EXONS_CMD = os.path.join( os.path.dirname( __file__ ), 
                                "./elements/exons/", "merge_exons.py" )

MERGE_WIGS_CMD = os.path.join( os.path.dirname( __file__ ), 
                              "./utilities/", "wiggle_merger.py" )
BUILD_GENELETS_CMD = os.path.join( os.path.dirname( __file__ ), 
                            "./elements/exons/", "build_genelets.py" )
DISC_CAGE_TSS_CMD = os.path.join( os.path.dirname( __file__ ), 
                              "./elements/distal_exons/", \
                                      "discover_tss_exons_from_cage_data.py" )
DISC_POLYA_TES_CMD = os.path.join( os.path.dirname( __file__ ), 
                              "./elements/distal_exons/", 
                                   "discover_tes_exons_from_polya_reads.py" )

BUILD_TRANSCRIPTS_CMD = os.path.join( os.path.dirname( __file__ ), 
                              "./elements/", "build_transcripts.py" )

MERGE_TRANSCRIPTS_CMD = os.path.join( os.path.dirname( __file__ ), 
                              "./analyze/", "merge_transcripts.py" )

TRANS_CMP_CMD = os.path.join( os.path.dirname( __file__ ), 
                              "./analyze/", "compare_annotations.py" )

EXTRACT_ELEMENTS_CMD = os.path.join( os.path.dirname( __file__ ), 
                         "./elements/", "extract_elements_from_transcripts.py" )




sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "./file_types/" ) )
from chrm_sizes import ChrmSizes

BaseDataTypes = [  "RNASeq", "CAGE", "PolyASeq", "Annotation"]
CHRM_LENS_DATA_TYPE = "chrm_sizes"
GENOME_DATA_TYPE = "genome"
RnaSeqDataTypes = [ "rnaseq_bam",
                    "rnaseq_rev_bam",
                    "rnaseq_cov_bedgraph",
                    "filtered_jns_gff",
                    "TF_jns_gff",
                    "jns_gff", 
                    "cage_fastq", 
                    "cage_cov_wig", 
                    "cage_tss_exons_gff", 
                    "exons_gff",
                    "TF_exons_gff",
                    "clustered_exons_gff",
                    "TF_clustered_exons_gff",
                    "polya_reads_gff",
                    "polya_tes_exons_gff",
                    "transcripts_gtf",
                    "TF_filt_element_transcripts_gtf",
                    "cdna_gtf",
                    "cdna_tes_gff",
                    "cdna_tss_gff",
                    "cdna_jns_gff",
                    "stats",
                    "TFE_filt_stats",
                    "annotation_gtf" ]

Element = namedtuple( "Element", ["data_type", "sample_type", \
                                  "sample_id", "strand", "fname", "status"] )
ElementType = namedtuple( "ElementType", \
                              ["data_type", "sample_type", \
                                   "sample_id", "strand"] )


def get_sample_type_name( sample_type ):
    if sample_type == '*':
        return "merged"
    else:
        return sample_type

get_sample_id_name = get_sample_type_name


def calc_chksum( fname, blocksize=8*1024*1024 ):
    """Calculate a checksum for fname, to verify that the file hasn't changed.

    """
    fp = open( fname )
    cksum = crc32( fp.read( blocksize ) )
    data = fp.read( blocksize )
    num_blocks = 1
    while num_blocks < MAX_NUM_CKSUM_BLOCKS and data != "":
        cksum = crc32( data, cksum )
        data = fp.read( blocksize )
        num_blocks += 1
    
    return cksum

class Elements( object ):    
    def _iter_elements_from_fp( self, elements_fp ):
        for line in elements_fp:
            # skip comment lines
            if line.startswith( "#" ): continue
            data = line.split()
            data.append( "finished" )
            yield Element( *data )
        
        return
    
    def add_element( self, element, dependencies ):
        c = self.conn.cursor()
        q = """INSERT INTO samples 
                  ( 'data_type', 'sample_type', 
                    'sample_id', 'strand', 
                    'fname', 'status' ) 
                  VALUES
                  ( ?, ?, ?, ?, ?, ?)
            """
        data = list(element)
        c.execute( q, data )
        
        q = """INSERT INTO dependencies
                  ( 'sample_fname', 'dependency_fname' )
                  VALUES
                  ( ?, ?)
            """
        data = [ element.fname, None]
        for dependency_fname in dependencies:
            data[1] = dependency_fname
            c.execute( q, data )
        
        c.close()
        self.conn.commit()
    
    def change_element_status( self, element, new_status ):
        q = """
        UPDATE samples
           SET status = '{1}' 
         WHERE data_type   = '{0.data_type}'
           AND sample_type = '{0.sample_type}'
           AND sample_id   = '{0.sample_id}'
           AND strand      = '{0.strand}';
        """.format( element, new_status )
        
        self.conn.execute( q )
        
        return
    
    def update_mod_time_and_cksum( self, fname ):
        chksum = calc_chksum( fname )
        last_mod_time = os.path.getmtime( fname )
        
        q = """
        UPDATE samples
           SET chksum = '{0}',
               last_mod_time = '{1}' 
         WHERE fname  = '{2}'
        """.format( chksum, last_mod_time, fname )
        self.conn.execute( q )
        self.conn.commit()

        return
    
    def get_dependencies( self, fname ):
        q = """
        SELECT dependency_fname, status
          FROM dependencies, samples
         WHERE dependencies.sample_fname = '{0}'
           AND dependencies.dependency_fname = samples.fname;
        """.format( fname )
        return list( self.conn.execute( q ) )

    def log_cmd( self, cmd_str ):
        q = """
        INSERT INTO commands_log VALUES ( ? );
        """
        self.conn.execute( q, ( cmd_str, ) )
        self.conn.commit()
        return
    
    def _init_db(self, new_elements):
        # turn on foreign keys
        self.conn.execute("PRAGMA foreign_keys = ON;" )
        self.conn.execute("PRAGMA journal_mode = MEMORY;")
        self.conn.execute("PRAGMA synchronous = OFF;")

        # create the commands log table
        self.conn.execute('''
        CREATE TABLE IF NOT EXISTS commands_log (
            commands        TEXT
        );
        ''')
        
        # create the file types table and insert them
        self.conn.execute('''
        CREATE TABLE IF NOT EXISTS data_types (
            data_type        TEXT,
            PRIMARY KEY( data_type )
        );
        ''')
        c = self.conn.cursor()
        for x in RnaSeqDataTypes:
            try:
                q = "INSERT INTO data_types VALUES ('{0}');".format(x)
                c.execute( q )
            except sqlite3.IntegrityError, inst:
                if str(inst) != 'column data_type is not unique': 
                    raise
        c.close()
        
        # create the sample table
        self.conn.execute('''
        CREATE TABLE IF NOT EXISTS samples (
             data_type       TEXT,
             sample_type     TEXT,
             sample_id       TEXT, 
             strand          TEXT CHECK( strand IN ('+', '-', '.' ) ),
             fname           TEXT,
             last_mod_time   FLOAT DEFAULT NULL,
             chksum          INT DEFAULT NULL,
             
             -- store whether the element has been built yet or not
             status      TEXT DEFAULT 'queued', 
             FOREIGN KEY (data_type) REFERENCES data_types(data_type),
             UNIQUE( data_type, sample_type, sample_id, strand ),
             UNIQUE( fname )
        );''')

        self.conn.execute('''
        CREATE TABLE IF NOT EXISTS dependencies (
            sample_fname        TEXT,
            dependency_fname    TEXT,
            PRIMARY KEY( sample_fname, dependency_fname ),
            FOREIGN KEY (sample_fname) REFERENCES samples(fname)
            --FOREIGN KEY (dependency_fname) REFERENCES data_types(fname)
        );
        ''')
        
        self.conn.commit()
        
        # insert the samples
        for element in new_elements:
            try:
                self.add_element( element, [] )
                self.update_mod_time_and_cksum( element.fname )
            except sqlite3.IntegrityError, inst:
                # only ignore unique constraints. We expect these 
                # because we are often reinserting the same base elements.
                if not str(inst).endswith( 'not unique' ):
                    raise
        
        return
    
    def get_elements_from_db( self, data_type, 
                              sample_type=None, 
                              sample_id=None,
                              strand=None ):
        assert data_type in RnaSeqDataTypes
        q_str = """
        SELECT  data_type, sample_type, sample_id, strand, fname, status
           FROM samples
          WHERE data_type   = '{0.data_type}'
        """
        if sample_type != None:
            q_str += " "*10 + "\n  AND sample_type = '{0.sample_type}'"
        if sample_id != None:
            q_str += " "*10 + "\n  AND sample_id   = '{0.sample_id}'"
        if strand != None:
            q_str += " "*10 + "\n  AND strand      = '{0.strand}'"
        q_str += ";"

        element_type = ElementType(data_type, sample_type, sample_id, strand)
        res = self.conn.execute( q_str.format( element_type ) ).fetchall()
        return map( Element._make, res )

    def get_elements_from_db_from_element( self, e ):
        return self.get_elements_from_db( \
            e.data_type, e.sample_type, e.sample_id, e.strand )

    def get_distinct_element_types_and_ids( self, data_types, get_merged=True ):
        if isinstance( data_types, str ):
            data_types = [ data_types, ]
        
        res = self.get_elements_from_db( data_types[0] )
        sample_types_and_ids = \
            set( (e.sample_type, e.sample_id) for e in res )

        for data_type in data_types[1:]:
            # get all of the distinct element types from the db
            res = self.get_elements_from_db( data_type )
            sample_types_and_ids.intersection_update( \
                (e.sample_type, e.sample_id) for e in res )
            
        if get_merged:
            return sample_types_and_ids
        else:
            return [ (s,t) for s,t in sample_types_and_ids \
                         if s not in '*M' and t not in '*M' ]
    
    def remove_element( self, fname ):
        print "REMOVING", fname
        # first, find the elements that depend on this elements
        # and remove them. 
        deps_q = "SELECT sample_fname \n" \
               + "  FROM dependencies \n" \
               + " WHERE dependency_fname = '{0}';"
        for dependent_fname, in self.conn.execute( deps_q.format( fname ) ):
            self.remove_element( dependent_fname )
        
        # next, remove this element from the deps table
        q = "DELETE FROM dependencies WHERE sample_fname = '{0}';".format(fname)
        self.conn.execute( q )

        # remove this element from the samples table
        q = "DELETE FROM samples WHERE fname = '{0}';".format( fname )
        self.conn.execute( q )
        
        self.conn.commit()
        
        return
    
    def element_has_changed( self, fname ):
        # make sure the file actually exists
        if not os.path.exists( fname ):
            return True
        
        # get the element
        q = "SELECT last_mod_time, chksum FROM samples\n" \
          + " WHERE fname = '{0}';".format( fname )
        res = list( self.conn.execute( q ) )
        
        #if the length of res is 0, then the element has
        # been removed from the DB and we can return True
        if len( res ) == 0: return True
        # otherwise, unpack the result
        assert len( res ) == 1
        db_last_mod_time, db_chksum = res[0]

        last_mod_time = os.path.getmtime( fname )
        # if the fname mod times are the same, the element is
        # the same ( in ther same sec to avoid bad float cmps )
        if round(db_last_mod_time - last_mod_time, 2) == 0:
            return False
        
        # otherwise, the element could still be the same, so
        # we verify with a chksum
        print "CALCULATING CHKSUM FOR", fname
        chksum = calc_chksum( fname )
        if chksum == db_chksum:
            # if the checksums are the same, update the access time in the DB
            return False
        
        return True
          
    def verify_db_elements( self  ):
        q = "SELECT fname, status FROM samples;"
        # if we can't run this query, it means that samples
        # hasn't been created yet so this is pointless
        try:
            res = self.conn.execute( q )
        except:
            return
        
        for fname, status in res:
            # if the element has changed, remove it
            # and the elements that depend on it from the 
            # database
            if status != 'finished' \
                    or self.element_has_changed( fname ):
                self.remove_element( fname )
                
        return
    
    def __init__( self, elements_fp ):
        # parse the input file, and add all of the input types
        new_elements = []
        self.chr_sizes = None
        self.genome_fname = None
        for element in self._iter_elements_from_fp( elements_fp ):
            # if this is a chrm sizes file, treat it specially
            if element.data_type == CHRM_LENS_DATA_TYPE:
                if self.chr_sizes != None:
                    raise ValueError, "Can only pass a single " \
                        + "chrm file sizes file."
                
                with open( element.fname ) as fp:
                    self.chr_sizes = ChrmSizes( fp )
            elif element.data_type == GENOME_DATA_TYPE:
                if self.genome_fname != None:
                    raise ValueError, "Can only pass a single genome file."
                self.genome_fname = element.fname
            else:
                new_elements.append( element )
        
        # make sure that we have access to the chrm sizes
        if self.chr_sizes == None:
            raise ValueError, "You must include chrm file sizes as input."

        # make sure that we have access to the genome 
        if self.genome_fname == None:
            raise ValueError, "You must include genome file as input."
        
        # create a sql lite database to store the info in
        db_fname = os.path.abspath( elements_fp.name ) + ".db"
        self.conn = sqlite3.connect(db_fname)

        # verify the db elements to remove anything that has been removed.
        # This prevents a key error for inew insertions
        self.verify_db_elements()

        # init the database tables if need be
        self._init_db( new_elements )

        # verify the newly added elements
        self.verify_db_elements()
        
        return
    
    def __str__( self ):
        rv = []
        rv.append(  "============= FULL DB =============" )
        rv.append( "Samples DB:" )
        for e in self.conn.execute( "SELECT * FROM samples;" ):
            rv.append( str( e ) )

        rv.append( "\nDependencies DB:" )
        for e in self.conn.execute( "SELECT * FROM dependencies \
                                     ORDER BY sample_fname;" ):
            rv.append( str(( os.path.basename(e[0]), \
                             os.path.basename(e[1]) )) )

        rv.append( "\nLogged Commands:" )
        for e in self.conn.execute( "SELECT * FROM commands_log;" ):
            rv.append( e[0] )
        
        rv.append( "============= END FULL DB =============" )
        return "\n".join( rv )

class Resource( int ):
    """Store the system resources that a process is using. 
    
    For now, we simply assume that this is the number of threads. In the 
    future, we may want, for instance, IO resources. 
    """
    pass

class ProcessServer( object ):
    """Serve and run processes.

    """
    def _update_running_commands( self ):
        # remove any finished processes
        indices_to_remove = []
        for i, ( p, cmd, n_rec ) in enumerate( self.running_commands ):
            # if the process is finished
            ret_code = p.poll()
            if ret_code != None:
                if ret_code != 0:
                    print >> sys.stderr, "Cmd: ", cmd
                    raise ValueError, "Invalid return code."
                indices_to_remove.append( i )
                self.available_resources += n_rec

        running_commands = []
        finished_commands = []
        for i, ( p, cmd, n_rec ) in enumerate(self.running_commands):
            if i not in indices_to_remove:
                running_commands.append( ( p, cmd, n_rec ) ) 
            else:
                finished_commands.append( cmd )
                # add the new element(s) to the database
                for output_e in cmd.iter_output_element_types():
                    # update the output element types to be 
                    # 'completed' in the DB
                    res = self.elements.get_elements_from_db_from_element(
                        output_e )
                    self.elements.change_element_status( \
                        output_e, "finished" )
                
                for op_fname in cmd.output_fnames:
                    self.elements.update_mod_time_and_cksum( op_fname )
        
        self.running_commands = running_commands
        
        if len( indices_to_remove ) > 0:
            if VERBOSE:
                for cmd in finished_commands:
                    print "FINISHED BUILDING: ".ljust(30), cmd.output_fnames
        
        return
    
    def _start_new_processes( self ):
        def dependencies_satisfied( cmd ):
            for output_fn in cmd.output_fnames:
                deps = self.elements.get_dependencies( output_fn )
                if any( status != 'finished' for fn, status in deps ):
                    return False
            
            return True
        
        # add any processes that we have space for
        cant_start_list = []
        num_new_cmds = 0
        while len( self.pending_commands ) > 0 \
                and self.available_resources > 0:
            new_cmd, ( min_res, max_res ) = self.pending_commands.pop(0)
            # if we don't have the recources to start this, or all
            # of the dependencies havn't been finished, put it on the
            # back of the queue
            if min_res > self.available_resources \
                    or not dependencies_satisfied( new_cmd ):
                cant_start_list.append((new_cmd, (min_res, max_res)))
            else:
                num_new_cmds += 1
                
                num_threads = min( max_res, self.available_resources )
                self.available_resources -= num_threads
                
                self.elements.log_cmd( new_cmd.get_cmd( num_threads )+ "\n" )
                
                p = Popen( new_cmd.get_cmd( num_threads ), shell=True )
                self.running_commands.append( (p, new_cmd, num_threads) )
                
                # change the status of the element to running
                for output_e in new_cmd.output_element_types:
                    self.elements.change_element_status( \
                        output_e, "running" )
            

        self.pending_commands.extend( cant_start_list )
        
        if num_new_cmds > 0:
            if VERBOSE:
                for ps, cmd, resources in self.running_commands[-num_new_cmds:]:
                    print "STARTED BUILDING:".ljust(30), cmd.output_fnames
        
        return
    
    def process_queue( self ):
        if STRESS_TEST_DEP_FINDER:
            random.shuffle( self.pending_commands  )
        
        while True:
            self._update_running_commands()
            
            self._start_new_processes()
            
            # if there's nothing left to do, break
            if len( self.pending_commands ) == 0 \
                    and len( self.running_commands ) == 0:
                break
            
            if len( self.running_commands ) == 0 \
                    and len( self.pending_commands ) > 0:
                print self.elements
                raise Exception, "There are pending commands but no running " \
                    + "commands. This probably indicates a dependency problem."
            
            time.sleep(0.5)
        
        return
    
    def __init__( self, elements, available_resources = Resource(1) ):
        # store a reference to elements. We will use these
        # to store the location of finished objects. 
        self.elements = elements
        
        # store the commands that havn't been started yet
        self.pending_commands = []
        
        # store the commands that are running
        self.running_commands = []
        
        # store the available thread resources
        self.available_resources = available_resources
        
    def add_process( self, cmd, min_resources, max_resources=None ):
        if max_resources == None:
            max_resources = min_resources
        
        # check to see if we already have these elements. If we do, then
        # don't add them. 
        ress = []
        for e in cmd.output_element_types:
            res = self.elements.get_elements_from_db_from_element( e )
            assert len( res ) <= 1
            ress.append( res )
        
        # if we got an element for every file, then we don't need 
        # to runthis command
        if all( len(res) == 1 for res in ress ):
            return
        
        # otherwise, remove the elements from the database that
        # are already there
        curr_ofnames = set()
        for res in ress:
            curr_ofnames.update( i.fname for i in res )
        
        for ofname in cmd.output_fnames:
            if ofname in curr_ofnames:
                self.elements.remove_element( ofname )
        
        # otherwise, add the element to the queue
        # first, add the elements to the database
        for e_t, ofname in zip( \
                cmd.output_element_types, cmd.output_fnames ):
            data = list( e_t )
            data.append( ofname )
            data.append( 'queued' )
            new_e = Element( *data )
            self.elements.add_element( new_e, cmd.dependencies )
        
        self.pending_commands.append( \
            ( cmd, (min_resources, max_resources) ) )
        
        return
    
class Cmd( object ):
    def __init__( self, cmd, output_element_types, 
                  output_fnames, dependency_fnames ):
        self.cmd = cmd
        self.output_element_types = output_element_types
        self.dependencies = dependency_fnames
        
        assert all(  isinstance( output_element_type, ElementType ) \
                     for output_element_type in output_element_types )
        self.output_fnames = output_fnames
        return
    
    def iter_output_element_types( self ):
        for t, f in zip( self.output_element_types, self.output_fnames ):
            yield t
        return
        
    def get_cmd( self, num_threads=1 ):
        return self.cmd.format( threads=num_threads )
    
    def __str__( self ):
        return self.cmd
    
    def print_all(self):
        data = (self.cmd, self.output_element_types, \
                    self.output_fnames, self.dependencies )
        return "\n".join( map( str, data ) )
    
def build_all_cov_wiggles( elements, process_server, derived_output_dir ):
    def build_wiggles_from_rnaseq_bam_fname(bam_fname, genome_lens, rev_strand):
        """Build the stranded bed graphs from bam_fname 
        python build_read_coverage_gtf.py 
               --chrm-sizes=/media/scratch/drosophila/dm3.chrom.sizes 
               -o test 
               --mapped-reads-fname L3_Carcass.574_BS311.bam
        
        """
        # find the output fname
        basename = os.path.basename( bam_fname )

        genome_lens_fname = genome_lens.get_tmp_fname( with_chr = False )

        dependencies = []
        
        # build the commands without strand
        cmd  = "python {0} ".format( BUILD_READ_CVG_CMD )
        cmd += "--chrm-sizes {0} ".format( \
            genome_lens.get_tmp_fname( with_chr = False ) )
        cmd += "--out-fname-prefix {0} ".format( \
            os.path.join( derived_output_dir, basename ) )
        cmd += "--mapped-reads-fname {0} ".format( bam_fname )
        if rev_strand:
            cmd += "--reverse-read-strand ".format( bam_fname )

        
        
        dependencies.extend( ( bam_fname, ) )
        output_fnames = []
        for strand in ( 'plus', 'minus' ):
            output_fnames.append( os.path.join( \
                    derived_output_dir, basename + ".%s.bedGraph" % strand ) )
        
        return cmd, output_fnames, dependencies
    
    # create the output directory
    try:
        os.makedirs( derived_output_dir )
    except OSError:
        # assume that the directory already exists, and continue
        pass
    
    def build_build_read_cov_cmds( res, cov_wig_fnames, rev_strand=False ):\
        # build all of the sample specific coverage files
        cmds = []
        for sample in res:
            cmd_str, output_fnames, dependencies = \
                    build_wiggles_from_rnaseq_bam_fname( \
                        sample.fname, elements.chr_sizes, rev_strand=rev_strand)

            cov_wig_fnames[ sample.sample_type ].extend( output_fnames )

            et1 = ElementType( \
               "rnaseq_cov_bedgraph", sample.sample_type, sample.sample_id, "+")
            et2 = ElementType( \
               "rnaseq_cov_bedgraph", sample.sample_type, sample.sample_id, "-")

            cmd = Cmd( cmd_str, [et1, et2], output_fnames, dependencies )
            process_server.add_process( cmd, Resource(1) )

    # first get all of the read bams
    cov_wig_fnames = defaultdict( list )
    res = elements.get_elements_from_db( "rnaseq_bam" )
    build_build_read_cov_cmds( res, cov_wig_fnames, rev_strand=False )
    res = elements.get_elements_from_db( "rnaseq_rev_bam" )
    build_build_read_cov_cmds( res, cov_wig_fnames, rev_strand=True )
    
    # add the wiggle merger process
    def add_merge_cmd( sample_type, cov_wig_fnames ):
        if sample_type == "*": 
            base_name = "rnaseq_cov"
        else:
            base_name =  "rnaseq_cov." + sample_type
                    
        #### Get the subprocess command info together
            
        # choose with chromosome because the input wiggles are always have a chr
        # ( to make them easier to visualize in UCSC ), and so we want to choose
        # a chrm sizes with chr as well
        chrm_sizes_fname = elements.chr_sizes.get_tmp_fname( with_chr = True )
        
        # get the output filename without strand
        out_fname_prefix = os.path.join( derived_output_dir, base_name )
        out_fname_template = out_fname_prefix + ".{strand}.bedGraph"
        
        # get the track name prefix
        track_name_prefix = base_name
        
        cmd_str_template = "python %s {0} --chrm-sizes {1} " % MERGE_WIGS_CMD
        cmd_str_template += "--out-fname-prefix {2} --track-name-prefix {3}"
        
        
        
        ##### get the output data types
        output_fnames = [ out_fname_template.format(strand="plus"), \
                          out_fname_template.format(strand="minus")    ]


        output_element_types = [ 
            ElementType( "rnaseq_cov_bedgraph", sample_type, "*", "+" ),  \
            ElementType( "rnaseq_cov_bedgraph", sample_type, "*", "-" )   \
        ]
        
        
        dependency_fnames = cov_wig_fnames
        
        cmd_str = cmd_str_template.format( \
            " ".join(cov_wig_fnames), chrm_sizes_fname, \
                out_fname_prefix, track_name_prefix )
        
        cmd = Cmd( cmd_str, output_element_types, \
                   output_fnames, dependency_fnames )
        
        process_server.add_process( cmd, Resource(1) )
        
    
    # build the merge sample type wiggles    
    for sample_type, fnames in cov_wig_fnames.iteritems():
        add_merge_cmd( sample_type, fnames )
    
    # add the merge everything sample
    add_merge_cmd( "*", list( itertools.chain( *cov_wig_fnames.values() ) ) )
    
    return 

def extract_all_junctions( elements, pserver, output_prefix ):
    try:
        os.makedirs( output_prefix )
    except OSError:
        # assume that the directory already exists, and continue
        pass
    
    def build_extract_junctions_cmd( bam_element, rev_strand=False ):
        bam_fn = bam_element.fname
        
        if bam_element.sample_id == "*" and bam_element.sample_type == "*":
            base_out_fname = "merged.jns.gff"
        else:
            sample_name = "merged_reps" if bam_element.sample_id == "*" \
                else bam_element.sample_id
            base_out_fname = bam_element.sample_type + "." \
                + sample_name + ".jns.gff"
        output_fname = os.path.join( output_prefix, base_out_fname )
        
        rev_strand_str = '--reverse-strand' if rev_strand else ""
        
        call_template = "python %s {0} --fasta {1} {2} > {3} " % EXTRACT_JNS_CMD
        call = call_template.format( \
            bam_fn, elements.genome_fname, rev_strand_str, output_fname ) 
        
        dependency_fnames = [ bam_fn, ]
        output_fnames = [ output_fname, ]
        output_element_type = ElementType( \
            "jns_gff", bam_element.sample_type, \
                bam_element.sample_id, bam_element.strand )
        output_element_types = [ output_element_type, ]
        
        return Cmd(call, output_element_types, output_fnames, dependency_fnames)
        
        
    # build all of the raw junctions
    res = elements.get_elements_from_db( "rnaseq_bam" )
    for e in res:
        cmd = build_extract_junctions_cmd( e, rev_strand = False )
        pserver.add_process( cmd, Resource(1) )

    res = elements.get_elements_from_db( "rnaseq_rev_bam" )
    for e in res:
        cmd = build_extract_junctions_cmd( e, rev_strand = True )
        pserver.add_process( cmd, Resource(1) )
    
    return

def build_high_quality_junctions( elements, pserver, output_prefix ):
    # check to see if there is a current merged junctions file. If there
    # is, then check to see what raw junctions it came from. If that list is 
    # the same as the list that is in the DB currently, then we don't need 
    # to do anything otherwise, remove it. 
    res = elements.get_elements_from_db( "jns_gff" )
    raw_jn_fnames = []
    for e in res:
        raw_jn_fnames.append( e.fname )
    
    merged_output_fname = os.path.join(output_prefix, "merged.filtered.jns.gff")
    
    # get all of the raw junctions fnames, and merge them
    call_template = "python %s --filter --fasta {0} {1} > {2} " \
        % MERGE_AND_FILTER_JNS_CMD
    call = call_template.format( \
        elements.genome_fname, " ".join( raw_jn_fnames ), merged_output_fname ) 
    
    dependency_fnames = raw_jn_fnames
    output_fnames = [ merged_output_fname, ]    
    output_element_types = [ ElementType("filtered_jns_gff", "*", "*", "."), ]
    
    cmd = Cmd( call, output_element_types, output_fnames, dependency_fnames )
    pserver.add_process( cmd, Resource(1) )
       
    return

def merge_sample_type_junctions( elements, pserver, output_prefix ):
    def merge_junctions( sample_type ):
        res = elements.get_elements_from_db( "jns_gff", sample_type )
        raw_jn_fnames = []
        for e in res:
            raw_jn_fnames.append( e.fname )

        merged_output_fname = os.path.join( \
            output_prefix, "%s.jns.gff" % sample_type )
        
        # get all of the raw junctions fnames, and merge them
        call_template = "python %s --fasta {0} {1} > {2} " \
            % MERGE_AND_FILTER_JNS_CMD
        call = call_template.format( \
            elements.genome_fname, " ".join( raw_jn_fnames ), \
            merged_output_fname ) 

        dependency_fnames = raw_jn_fnames
        output_fnames = [ merged_output_fname, ]    
        output_element_types = [ \
            ElementType("filtered_jns_gff", sample_type, "*", "."), ]

        cmd = Cmd(call, output_element_types, output_fnames, dependency_fnames)
        pserver.add_process( cmd, Resource(1) )
    
    e_types = elements.get_distinct_element_types_and_ids( \
        "jns_gff", get_merged=False )

    distinct_sample_types = sorted(set( _[0] for _ in e_types ))
    for sample_type in distinct_sample_types:
        merge_junctions( sample_type )
    
    return

def intersect_all_junctions( elements, pserver, output_prefix ):
    try:
        os.makedirs( output_prefix )
    except OSError:
        # assume that the directory already exists, and continue
        pass
    
    def build_filter_junctions_cmd( jn_element, merged_filtered_jn_fn ):
        jn_fn = jn_element.fname
        
        sample_name = "merged_reps" if jn_element.sample_id == "*" \
            else jn_element.sample_id
        base_out_fname = jn_element.sample_type + "." \
            + sample_name + ".filtered.jns.gff"
        output_fname = os.path.join( output_prefix, base_out_fname )
        
        call_template = "python %s {0} {1} > {2} " % INTERSECT_JNS_CMD
        call = call_template.format(merged_filtered_jn_fn, jn_fn, output_fname)
        
        dependency_fnames = [ jn_fn, merged_filtered_jn_fn ]
        output_fnames = [ output_fname, ]
        output_element_type = ElementType( \
            "filtered_jns_gff", jn_element.sample_type, \
                jn_element.sample_id, jn_element.strand )
        output_element_types = [ output_element_type, ]
        
        return Cmd(call, output_element_types, output_fnames, dependency_fnames)
    
    
    # build all of the raw junctions
    res = elements.get_elements_from_db( "jns_gff" )
    merged_filtered_jn_fn = list( elements.get_elements_from_db( \
        "filtered_jns_gff", "*", "*" ) )[0].fname
    for e in res:
        cmd = build_filter_junctions_cmd( e, merged_filtered_jn_fn )
        pserver.add_process( cmd, Resource(1) )

    return



def build_all_exon_files( elements, pserver, output_prefix ):
    try:
        os.makedirs( output_prefix )
    except OSError:
        # assume that the directory already exists, and continue
        pass
        
    """
    python ../get_elements/find_exons.py 
    ../raw_elements/all_tissues.plus.wig 
    ../raw_elements/all_tissues.minus.wig 
    ../raw_elements/all_tissues.merged_junctions.gff 
    ../raw_elements/dm3.chrom.sizes 
    -o exons.all_tissues.gff 
    """
    def build_find_exons_cmds( sample_type, sample_id ):
        call_template = "python {0} {1} {2} {3} -o {4}"
        
        # find all of the merged rnaseq coverage files
        res = elements.get_elements_from_db( \
            "rnaseq_cov_bedgraph", sample_type, sample_id )
        bedgraph_fnames = [ i.fname for i in res ]
        if len( bedgraph_fnames ) != 2:
            raise ValueError, "Can't find the rnaseq coverage files " \
                              + "necessary to build exons from."

        # find the merged, filtered jns file
        if USE_MERGED_JNS_FOR_EXONS:
            res = elements.get_elements_from_db( \
                "filtered_jns_gff", "*" , "*", "." )
        else:
            res = elements.get_elements_from_db( \
                "filtered_jns_gff", sample_type , sample_id, "." )
            
        if len( res ) < 1:
            raise ValueError, "Can't find the jns gff file " \
                              + "necessary to build exons."
        assert len( res ) == 1
        jns_fname = res[0].fname

        chrm_sizes_fname = elements.chr_sizes.get_tmp_fname( with_chr = True )

        output_fname = os.path.join( output_prefix, \
            "discovered_exons.%s.%s.gff" % ( sample_type, sample_id ) )
        output_fname = output_fname.replace( ".*", "" )

        call = call_template.format( FIND_EXONS_CMD, " ".join(bedgraph_fnames),\
                            jns_fname, chrm_sizes_fname,   \
                            output_fname, )

        output_fnames = [ output_fname, ]
        output_element_types = [ElementType( \
                "exons_gff", sample_type, sample_id, "." ),]
        dependency_fnames = [ jns_fname ]
        dependency_fnames.extend( bedgraph_fnames )

        cmd = Cmd(call, output_element_types, output_fnames, dependency_fnames)
        pserver.add_process( cmd, Resource(1), Resource(1) )

    def build_merge_exons_cmds( sample_type=None ):
        # find all of the merged rnaseq coverage files
        res = elements.get_elements_from_db( \
            "exons_gff", sample_type )
        exon_fnames = [ i.fname for i in res \
                           if i.sample_id not in '*M' \
                            and i.sample_type not in '*M']
        
        op_base_fname = "discovered_exons.%s.gff" % sample_type \
            if sample_type != None else "discovered_exons_merged.gff"
        output_fname = os.path.join( output_prefix, op_base_fname )
        
        call = "python %s %s > %s" % ( \
            MERGE_EXONS_CMD, " ".join(exon_fnames), output_fname )

        if sample_type == None: 
            sample_type = 'M'
        output_element_types = [ElementType( \
                "exons_gff", sample_type, "M", "." ),]
        
        dependency_fnames = exon_fnames

        cmd = Cmd( call, output_element_types, \
                       [ output_fname, ], dependency_fnames )
        
        pserver.add_process( cmd, Resource(1), Resource(1) )
    
    
    # build the merged sample exons
    build_find_exons_cmds( "*", "*" )
    
    if BUILD_SAMPLE_SPECIFIC_EXONS:
        sample_types_and_ids = elements.get_distinct_element_types_and_ids( \
            [ "filtered_jns_gff", "rnaseq_cov_bedgraph" ], get_merged=False )
        
        # buidl the build eoxn cmds
        for sample_type, sample_id in sample_types_and_ids:
            build_find_exons_cmds( sample_type, sample_id )    

        # build the merge exons cmds
        for sample_type in set( t for t, i in sample_types_and_ids ):
            build_merge_exons_cmds( sample_type )    

        # build the merged exons
        build_merge_exons_cmds( )    

    
    return

def cluster_exons( elements, pserver, output_prefix, \
                       input_exon_type="exons_gff",  \
                       input_jn_type="filtered_jns_gff" ):
    """
    python build_genelets.py exons_gff junctions_gff --out_fname = OUT
    """
    if input_exon_type == 'exons_gff':
        assert input_jn_type == "filtered_jns_gff"
        output_element_type = "clustered_exons_gff"
        output_base_prefix = "clustered_exons"
    elif input_exon_type == 'TF_exons_gff':
        assert input_jn_type == "TF_jns_gff"
        output_element_type = "TF_clustered_exons_gff"
        output_base_prefix = "TF_filt_clustered_exons"
    else:
        assert False
    
    def cluster_exons( sample_type, sample_id=None ):
        # get the merged exons
        if sample_id == None:
            sample_id = 'M'
        res = elements.get_elements_from_db( \
            input_exon_type, sample_type, sample_id )
        assert len( res ) == 1
        merged_exons_fname = res[0].fname

        # get the merged junctions
        jn_sample_id = '*' if sample_id == 'M' else sample_id
        res = elements.get_elements_from_db( \
            input_jn_type, sample_type, jn_sample_id )
        assert len( res ) == 1
        merged_jns_fname = res[0].fname

        op_fname = os.path.join( output_prefix, \
            "%s.%s.%s.gff" % ( output_base_prefix, sample_type, sample_id ) )
        op_fname = op_fname.replace( ".*", "" )
        
        cmd_tem = "python {0} --exons {1} --junctions-gff {2} --out_fname {3}"
        cmd_str = cmd_tem.format( \
            BUILD_GENELETS_CMD, merged_exons_fname, merged_jns_fname, op_fname )

        op_element_types = [ \
            ElementType( output_element_type, sample_type, sample_id, "." ), ]
        op_fnames = [ op_fname,]
        
        dependencies = [merged_exons_fname, merged_jns_fname]
        cmd = Cmd( cmd_str, op_element_types, op_fnames, dependencies )
        pserver.add_process( cmd, Resource(1) )

    if BUILD_SAMPLE_SPECIFIC_EXONS:
        sample_types_and_ids = elements.get_distinct_element_types_and_ids( \
            [ input_exon_type, input_jn_type ], get_merged=False )    
        
        for sample_type, sample_id in sample_types_and_ids:
            cluster_exons( sample_type, sample_id )    
        
        sample_types = set( st for st, sid in sample_types_and_ids )
        for sample_type, sample_id in sample_types_and_ids:
            cluster_exons( sample_type )    
    
    cluster_exons( '*', '*' )
    
    return

def build_cage_tss_exons( elements, pserver, output_prefix, filter_with_rnaseq=False ):
    """
    python ../get_elements/get_tss_exons.py
        [--tss_reads_wigs [TSS_READS_WIGS [TSS_READS_WIGS ...]]]
        [--tss_reads_nc_wigs [TSS_READS_NC_WIGS [TSS_READS_NC_WIGS ...]]]
        [--clustered_exons_gff CLUSTERED_EXONS_GFF]
        [--chrm_sizes_fname CHRM_SIZES_FNAME]
        [--out_fname OUT_FNAME] [--verbose
        
        [--transcripts_gtfs [TRANSCRIPTS_GTFS [TRANSCRIPTS_GTFS ...]]]
        [--tss_gffs [TSS_GFFS [TSS_GFFS ...]]]

    """
    try:
        os.makedirs( output_prefix )
    except OSError:
        # assume that the directory already exists, and continue
        pass

    def build_mk_cage_cmd( sample_type ):
        # get the cage coverage wiggles
        ress = elements.get_elements_from_db( "cage_cov_wig", "*", "*"  )
        if len( ress ) != 2:
            raise ValueError, "Can't build CAGE TSSs without CAGE read cov files."
        cage_cov_fnames = [ res.fname for res in ress ]

        # get the nc wigs ( the merged RNA seq data )
        if filter_with_rnaseq:
            ress = elements.get_elements_from_db( \
                "rnaseq_cov_bedgraph", sample_type, "*" )
            cage_NC_cov_fnames = [ res.fname for res in ress ]
        else:
            cage_NC_cov_fnames = []

        # get the exons
        if USE_MERGED_EXONS_FOR_TSS_EXONS: 
            exons_sample_type = "*"
        else:
            exons_sample_type = sample_type        
        exons_sample_id = "M" if exons_sample_type != '*' else '*'
        
        res = elements.get_elements_from_db( \
            "clustered_exons_gff", exons_sample_type, exons_sample_id )
        assert len( res ) == 1
        clustered_merged_exons_fname = res[0].fname

        res = elements.get_elements_from_db( "filtered_jns_gff", sample_type, "*" )
        assert len( res ) == 1
        jns_fname = res[0].fname
        
        output_fname = os.path.join( \
            output_prefix, "discovered_cage_tss_exons.%s.gff" % sample_type)
        output_fname = output_fname.replace(".*", "")

        call  = "python %s " %  DISC_CAGE_TSS_CMD
        call += "--chrm-sizes {0} ".format( \
            elements.chr_sizes.get_tmp_fname( with_chr = False ) )
        call += "--clustered-exons {0} ".format( clustered_merged_exons_fname )

        call += "--cage-read-wigs {0} ".format( " ".join(cage_cov_fnames) )

        call += "--junctions-gff {0} ".format( jns_fname )

        # cage read nc wigs is not implemented
        #if len( cage_NC_cov_fnames ) > 0:
        #    call += "--cage-read-nc-wigs {0} ".format( " ".join(cage_NC_cov_fnames))
        #call += "--smooth-coverage-for-quantile "
        call += " > {0} ".format( output_fname ) 
        
        op_element_types = [ \
            ElementType( "cage_tss_exons_gff", sample_type, "*", "." ), ]
        op_fnames = [ output_fname,]
        dependencies = [ clustered_merged_exons_fname, ]
        dependencies.extend( cage_cov_fnames )
        dependencies.extend( cage_NC_cov_fnames )
        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        pserver.add_process( cmd, Resource(1) )
    
    ip_etypes = [ "clustered_exons_gff", "rnaseq_cov_bedgraph" ]
    element_types = sorted( set( i[0] for i in  \
            elements.get_distinct_element_types_and_ids( ip_etypes, False ) ) )

    if USE_SAMPLE_SPECIFIC_CAGE:
        for element_type in element_types:
            build_mk_cage_cmd( element_type )
    else:
        build_mk_cage_cmd( "*" )
    
    build_mk_cage_cmd( "*" )
    
    return


def build_polya_tes_exons( elements, pserver, \
                               output_prefix, filter_with_rnaseq=False ):
    """
    python ~/slide/trunk/slide/distal_exons/discover_tes_exons_from_polya_reads.py 
    --chrm-sizes /media/scratch/genomes/drosophila/dm3.chrom.sizes.nochr 
    --clustered-exons exons/clustered_exons.gff 
    --polya-read-gffs chr4_raw/chr4.polyas.gff 
    --polya-read-nc-wigs /media/scratch/NEW_RNA_ELEMENTS/read_cov_wigs/rnaseq_cov.plus.bedGraph 
                         /media/scratch/NEW_RNA_ELEMENTS/read_cov_wigs/rnaseq_cov.minus.bedGraph  
    --verbose 
         > test.tes_exons.gff
    """
    try:
        os.makedirs( output_prefix )
    except OSError:
        # assume that the directory already exists, and continue
        pass
    
    # get the polya reads
    ress = elements.get_elements_from_db( "polya_reads_gff", "*", "*"  )
    if len( ress ) == 0:
        raise ValueError, "Can't build POLYA TESs without POLYA reads file."
    polya_read_fnames = [ res.fname for res in ress ]
    
    if filter_with_rnaseq:
        # get the nc wigs ( the merged RNA seq data )
        ress = elements.get_elements_from_db( "rnaseq_cov_bedgraph", "*", "*" )
        tes_reads_nc_wig_fnames = [ res.fname for res in ress ]
    else:
        tes_reads_nc_wig_fnames = []
    
    # get the merged exons
    res = elements.get_elements_from_db( "clustered_exons_gff", "*", "*" )
    assert len( res ) == 1
    clustered_merged_exons_fname = res[0].fname

    res = elements.get_elements_from_db( "filtered_jns_gff", "*", "*" )
    assert len( res ) == 1
    jns_fname = res[0].fname
    
    output_fname = os.path.join(output_prefix, "discovered_polya_tes_exons.gff")
    
    call  = "python %s " %  DISC_POLYA_TES_CMD
    call += "--chrm-sizes {0} ".format( \
        elements.chr_sizes.get_tmp_fname( with_chr = False ) )
    call += "--clustered-exons {0} ".format( clustered_merged_exons_fname )

    call += "--polya-read-gffs {0} ".format( " ".join(polya_read_fnames) )
    if len( tes_reads_nc_wig_fnames ) > 0:
        call += "--polya-read-nc-wigs {0} ".format( \
            " ".join(tes_reads_nc_wig_fnames) )
    call += "--junctions-gff {0} ".format( jns_fname )
    call += " --out-fname {0} ".format( output_fname ) 
    
    op_element_types = [ \
        ElementType( "polya_tes_exons_gff", "*", "*", "." ), ]
    op_fnames = [ output_fname,]
    dependencies = [ clustered_merged_exons_fname, ]
    dependencies.extend( polya_read_fnames )
    dependencies.extend( tes_reads_nc_wig_fnames )
    cmd = Cmd( call, op_element_types, op_fnames, dependencies )
    pserver.add_process( cmd, Resource(1) )

    return

def build_transcripts( elements, pserver, output_prefix, use_TF_elements=False ):
    """
    build_transcripts.py [-h] --exons EXONS [EXONS ...] 
                              --tss TSS [TSS ...]
                              --tes TES [TES ...] 
                              --junctions JUNCTIONS [JUNCTIONS ...] 
                              [--threads THREADS]
                              [--out-fname OUT_FNAME] 
                              [--log-fname LOG_FNAME]
                              [--verbose]

    """
    try:
        os.makedirs( output_prefix )
    except OSError:
        # assume that the directory already exists, and continue
        pass

    if use_TF_elements:
        exon_input_type = 'TF_clustered_exons_gff'
        jn_input_type = 'TF_jns_gff'
        transcript_output_type = "TF_filt_element_transcripts_gtf"
    else:
        exon_input_type = 'clustered_exons_gff'
        jn_input_type = 'filtered_jns_gff'
        transcript_output_type = "transcripts_gtf"
        
    def build_transcripts_for_sample( sample_type, sample_id ):
        if sample_type == '*': return
        
        # get the exons
        if not USE_MERGED_EXONS_FOR_TRANSCRIPTS and BUILD_SAMPLE_SPECIFIC_EXONS:
            ress = elements.get_elements_from_db( \
                exon_input_type, sample_type, sample_id  )
        else:
            ress = elements.get_elements_from_db( \
                exon_input_type, "*", "*"  )
            
        if len( ress ) < 1: return
        assert len( ress ) == 1
        exons_fname = ress[0].fname

        # get tss exons
        cage_st = sample_type if USE_SAMPLE_SPECIFIC_CAGE else "*"
        ress = elements.get_elements_from_db("cage_tss_exons_gff", cage_st, "*" )
        assert len( ress ) == 1
        
        tss_exons_fnames = [ ress[0].fname, ]
        ress = elements.get_elements_from_db( "cdna_tss_gff" )
        tss_exons_fnames.extend( e.fname for e in ress )

        # get tes exons
        ress = elements.get_elements_from_db( "polya_tes_exons_gff", "*", "*"  )
        assert len( ress ) == 1
        tes_exons_fnames = [ ress[0].fname, ]
        ress = elements.get_elements_from_db( "cdna_tes_gff" )
        tes_exons_fnames.extend( e.fname for e in ress )

        # get junctions
        ress = elements.get_elements_from_db( \
            jn_input_type, sample_type, sample_id  )
        assert len( ress ) == 1
        jns_fname = ress[0].fname

        sample_type_name = "merged_input" if sample_type == "*" else sample_type
        sample_id_name = "merged_input" if sample_id == "*" else sample_id
        suffix = ".transcripts.gtf" if not use_TF_elements \
            else ".tfelements.transcripts.gtf"
        op_fname = os.path.join( output_prefix, sample_type_name + "." \
                                     + sample_id_name + suffix )
        
        log_fname = op_fname + ".log"
        
        call = "python %s --exons {0} --tss {1} --tes {2} --junctions {3} " \
                      + " --out-fname {4} --log-fname {5}"
        call = call % BUILD_TRANSCRIPTS_CMD
        
        call = call.format( exons_fname, " ".join( tss_exons_fnames ), \
                                " ".join( tes_exons_fnames ),  \
                                jns_fname, op_fname, log_fname  )
        call +=  " --threads={threads}"
        
        op_element_types = [ \
            ElementType( transcript_output_type, sample_type, sample_id, "." ), ]
        op_fnames = [ op_fname,]

        dependencies = [ exons_fname, jns_fname ]
        dependencies.extend( tss_exons_fnames )
        dependencies.extend( tes_exons_fnames )

        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        pserver.add_process( cmd, Resource(1), Resource(8) )
        
        return
    
    ress = elements.get_elements_from_db( jn_input_type )
    for e in ress[:]:
        build_transcripts_for_sample( e.sample_type, e.sample_id  )
    
    return

def merge_transcripts( elements, pserver, output_prefix, \
                           input_e_type='transcripts_gtf' ):
    try:
        os.makedirs( output_prefix )
    except OSError:
        # assume that the directory already exists, and continue
        pass
    
    def add_merge_transcripts_sample_cmd(  \
                sample_type, input_fnames, \
                op_e_type="transcripts_gtf" ):
        if op_e_type == "TF_filt_element_transcripts_gtf":
            suffix = ".TFfilt_transcripts.gtf"
        elif op_e_type == 'transcripts_gtf':
            suffix = ".transcripts.gtf"
        else:
            assert False
        
        if sample_type == "*": 
            raise ValueError, "Can't merge merged transcripts."

        sample_type_name = "merged" if sample_type == "M" else sample_type
        op_fname = os.path.join(output_prefix, \
                                    sample_type_name + suffix)
        
        call = "python %s %s > %s" \
            % ( MERGE_TRANSCRIPTS_CMD, " ".join( input_fnames ), op_fname )
        
        op_element_types = [ \
            ElementType( op_e_type, sample_type, "M", "." ), ]
        op_fnames = [ op_fname,]
        dependencies = input_fnames
        
        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        pserver.add_process( cmd, Resource(1) )
    
    # add the merge all command first, since it takes the 
    # longest to complete
    in_es = elements.get_elements_from_db( input_e_type )
    input_fnames = [ e.fname for e in in_es if e.sample_id not in "*M"]
    if len( input_fnames ) > 0:
        add_merge_transcripts_sample_cmd( "M", input_fnames )
    
    # get all of the sample types, and add them
    sample_types = set()
    for e in elements.get_elements_from_db( input_e_type ):
        if e.sample_type != "*": sample_types.add( e.sample_type )

    for sample_type in sample_types:
        # get the input elements filenames
        in_es = elements.get_elements_from_db( input_e_type, sample_type )
        input_fnames = [ e.fname for e in in_es ]
        add_merge_transcripts_sample_cmd( sample_type, input_fnames )

    return

    
def run_all_slide_compares(elements, pserver, output_prefix, build_maps=False):
    assert build_maps == False
    
    try:
        os.makedirs( output_prefix )
    except OSError:
        # assume that the directory already exists, and continue
        pass
        
    def run_slide_compare( transcript_element, annotation_fname, \
                               output_e_type='stats' ):
        sample_type_name = transcript_element.sample_type
        if transcript_element.sample_type == "*": 
            sample_type_name = "merged_input"
        if transcript_element.sample_type == "M":
            sample_type_name = "merged"

        sample_id_name = transcript_element.sample_id
        if transcript_element.sample_id == "*": 
            sample_id_name = "merged_input"
        if transcript_element.sample_type == "M":
            sample_id_name = "merged"

        output_files_prefix = os.path.join( \
            output_prefix, "%s.%s.%s" % ( \
                sample_type_name, sample_id_name, output_e_type ) )
        
        call = "python %s {0.fname} {1} --out-prefix {2}" % TRANS_CMP_CMD
        call = call.format( transcript_element, \
                                annotation_fname, output_files_prefix  )
        call +=  " --threads={threads}"
        
        op_element_types = [ \
            ElementType( output_e_type, transcript_element.sample_type, \
                         transcript_element.sample_id, "." ), ]
        
        op_fnames = [ output_files_prefix + ".stats",]
        dependencies = [ transcript_element.fname, annotation_fname ]
        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        pserver.add_process( cmd, Resource(1) )

        return
    
    annotation_elements = elements.get_elements_from_db( "annotation_gtf" )
    ann_fnames = [ e.fname for e in annotation_elements ]

    # run slide compare for every combination of transcript and annotation
    transcript_elements = elements.get_elements_from_db( "transcripts_gtf" )
    for transcript_element in transcript_elements:
        for ann_fname in ann_fnames:
            run_slide_compare( transcript_element, ann_fname )

    # run for transcript element filtered transcroipts
    transcript_elements = elements.get_elements_from_db( \
        "TF_filt_element_transcripts_gtf" )
    for transcript_element in transcript_elements:
        for ann_fname in ann_fnames:
            run_slide_compare( transcript_element, ann_fname, "TFE_filt_stats" )

    
    return


def extract_cdna_elements( elements, pserver, tss_op_prefix, tes_op_prefix ):
    def build_extract_cdnas_cmd( input_element ):
        input_fname = input_element.fname
        
        sample_type_name = get_sample_type_name( input_element.sample_type )
        sample_id_name = get_sample_id_name( input_element.sample_id )
        op_fname_prefix = "%s.%s" % (sample_type_name, sample_id_name )

        tss_op_fname = os.path.join( \
            tss_op_prefix, op_fname_prefix + ".cdna_tss_exons.gff" )
        tes_op_fname = os.path.join( \
            tes_op_prefix, op_fname_prefix + ".cdna_tes_exons.gff" )

        call = "python %s --tss-fname {0} --tes-fname {1} {2}" \
            % EXTRACT_ELEMENTS_CMD
        call = call.format( tss_op_fname, tes_op_fname, input_fname )

        op_element_types = [ \
            ElementType( "cdna_tss_gff", input_element.sample_type, \
                         input_element.sample_id, "." ),
            ElementType( "cdna_tes_gff", input_element.sample_type, \
                         input_element.sample_id, "." ) 
            ]
        
        op_fnames = [ tss_op_fname, tes_op_fname ]
        
        dependencies = [ input_fname, ]
        
        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        pserver.add_process( cmd, Resource(1) )
        
    
    # get all of the cdna elements
    for e in elements.get_elements_from_db( "cdna_gtf" ):
        build_extract_cdnas_cmd( e )
    
    return

def extract_transcript_filtered_elements( \
        elements, pserver, exons_op_prefix, jns_op_prefix ):
    def build_extract_TF_filtered_cmd( input_element ):
        input_fname = input_element.fname
        
        sample_type_name = get_sample_type_name( input_element.sample_type )
        sample_id_name = get_sample_id_name( input_element.sample_id )
        op_fname_prefix = "%s.%s" % (sample_type_name, sample_id_name )

        exons_op_fname = os.path.join( \
            exons_op_prefix, "trans_filt." + op_fname_prefix + ".gff" )
        jns_op_fname = os.path.join( \
            jns_op_prefix, "trans_filt." + op_fname_prefix + ".gff" )

        call = "python %s --exon-fname {0} --jn-fname {1} {2}" \
            % EXTRACT_ELEMENTS_CMD
        call = call.format( exons_op_fname, jns_op_fname, input_fname )

        op_element_types = [ \
            ElementType( "TF_exons_gff", input_element.sample_type, \
                         input_element.sample_id, "." ),
            ElementType( "TF_jns_gff", input_element.sample_type, \
                         input_element.sample_id, "." ) 
            ]
        
        op_fnames = [ exons_op_fname, jns_op_fname ]
        
        dependencies = [ input_fname, ]
        
        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        pserver.add_process( cmd, Resource(1) )
        
        return exons_op_fname, jns_op_fname
        
    
    # get all of the cdna elements
    for e in elements.get_elements_from_db( "transcripts_gtf" ):
        build_extract_TF_filtered_cmd( e )
    
    return

def main():
    elements_fname = sys.argv[1]
    num_threads = int( sys.argv[2] )
    base_dir = sys.argv[3]
    
    with open( elements_fname ) as fp:
        elements = Elements( fp )
    
    pserver = ProcessServer( elements, Resource( num_threads )  )
    
    extract_all_junctions( elements, pserver, base_dir + "junctions/" )
    merge_sample_type_junctions( elements, pserver, base_dir + "junctions/" )
    build_high_quality_junctions( elements, pserver, base_dir + "junctions/" )
    intersect_all_junctions( elements, pserver, base_dir + "junctions/" )
        
    build_all_cov_wiggles( elements, pserver, base_dir + "read_cov_bedgraphs" )
    
    build_all_exon_files( elements, pserver, base_dir + "exons" )
    
    cluster_exons( elements, pserver, base_dir + "exons", \
                       "exons_gff", "filtered_jns_gff" )

    build_polya_tes_exons( elements, pserver, base_dir + "tes_exons" )
    build_cage_tss_exons( elements, pserver, base_dir + "tss_exons" )
    extract_cdna_elements( \
        elements, pserver, base_dir + "tss_exons", base_dir + "tes_exons" )
    
    build_transcripts( elements, pserver, base_dir + "transcripts" )

    if USE_TF_ELEMENTS:
        cluster_exons( elements, pserver, base_dir + "exons", \
                       "TF_exons_gff", "TF_jns_gff" )

        extract_transcript_filtered_elements( \
            elements, pserver, base_dir + "exons", base_dir + "junctions" )

        build_transcripts( elements, pserver, base_dir + "transcripts", True )

    merge_transcripts( elements, pserver, base_dir + "transcripts" )

    run_all_slide_compares( elements, pserver, base_dir + "stats" )

    pserver.process_queue()
    
    return


main()

