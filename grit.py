# Copyright (c) 2011-2012 Nathan Boley

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

USE_POLYA_SIGNAL = False

BUILD_SAMPLE_SPECIFIC_EXONS = True
USE_SAMPLE_SPECIFIC_CAGE = True and BUILD_SAMPLE_SPECIFIC_EXONS
USE_MERGED_RNASEQ_FOR_EXONS = True and BUILD_SAMPLE_SPECIFIC_EXONS

USE_MERGED_JNS_FOR_EXONS = False
USE_MERGED_CAGE_SIGNAL_FOR_EXONS = False
USE_MERGED_POLYA_SIGNAL_FOR_EXONS = True

USE_MERGED_EXONS_FOR_TRANSCRIPTS = False
USE_HEURISTIC_FILTERING = False

# whether or not to include the merged_input transripts
# in the final merged file
USE_MERGED_INPUT=True

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
                                "./elements/", "find_exons.py" )

MERGE_EXONS_CMD = os.path.join( os.path.dirname( __file__ ), 
                                "./elements/exons/", "merge_exons.py" )

MERGE_BEDGRAPHS_CMD = os.path.join( os.path.dirname( __file__ ), 
                              "./utilities/", "merge_bedgraphs.py" )
 
BUILD_TRANSCRIPTS_CMD = os.path.join( os.path.dirname( __file__ ), 
                              "./elements/", "build_transcripts.py" )

MERGE_TRANSCRIPTS_CMD = os.path.join( os.path.dirname( __file__ ), 
                              "./elements/transcripts/", "merge_transcripts.py")

TRANS_CMP_CMD = os.path.join( os.path.dirname( __file__ ), 
                              "./analyze/", "compare_annotations.py" )

EXTRACT_ELEMENTS_CMD = os.path.join( os.path.dirname( __file__ ), 
                         "./elements/", "extract_elements_from_transcripts.py" )

BUILD_FL_DIST_CMD = os.path.join( os.path.dirname( __file__ ), 
                        "./sparsify/", "frag_len.py" )

SPARSIFY_TRANSCRIPTS_CMD = os.path.join( 
    os.path.dirname( __file__ ), "./sparsify/", "sparsify_transcripts.py" )

FIND_ORFS_CMD = os.path.join( 
    os.path.dirname( __file__ ), "./proteomics/", "ORF_finder.py" )

SPLIT_MERGED_TRANS_CMD = os.path.join( 
    os.path.dirname( __file__ ), "./elements/transcripts/", 
    "split_merged_transcripts.py" )

HEURISTIC_FILTER_TRANS_CMD = os.path.join( 
    os.path.dirname( __file__ ), "./sparsify/", "filter_transcript_set.py" )

RENAME_TRANS_CMD = os.path.join( 
    os.path.dirname( __file__ ), "./utilities/", 
    "group_transcripts_by_reference.py" )

EST_EXON_EXP_CMD = os.path.join( 
    os.path.dirname( __file__ ), "expression", 
    "estimate_exon_expression.py" )

AGGREGATE_EXP_CMD = os.path.join( 
    os.path.dirname( __file__ ), "expression", 
    "make_csv_from_gtfs.py" )

sys.path.insert( 0, os.path.join( os.path.dirname( __file__ ), \
                                      "./file_types/" ) )
from chrm_sizes import ChrmSizes

BaseDataTypes = [  "RNASeq", "CAGE", "PolyASeq", "Annotation"]
CHRM_LENS_DATA_TYPE = "chrm_sizes"
GENOME_DATA_TYPE = "genome"
RnaSeqDataTypes = [ "rnaseq_polyaplus_bam",
                    "rnaseq_polyaplus_rev_bam",
                    "rnaseq_polyaplus_unstranded_bam",
                    "rnaseq_total_bam",
                    "annotation_gtf",
                    "gene_merges_white_list",
                    "rnaseq_cov_pair1_bedgraph",
                    "rnaseq_cov_pair2_bedgraph",
                    "filtered_jns_gff",
                    "jns_gff", 
                    "cage_fastq", 
                    "cage_cov_wig", 
                    "cage_tss_exons_gff", 
                    "exons_gff",
                    "polya_reads_gff",
                    "polya_tes_exons_gff",
                    "single_exon_genes_gff",
                    "cand_polya_sites",
                    "transcripts_gtf",
                    "sparse_transcripts_gtf",
                    "cdna_gtf",
                    "cdna_tes_gff",
                    "cdna_tss_gff",
                    "cdna_jns_gff",
                    "fl_dist",
                    "stats",
                    "CDS_transcripts_gtf",
                    "CDS_renamed_transcripts_gtf",
                    "heur_filtered_transcripts_gtf",
                    "heur_filtered_renamed_transcripts_gtf",
                    "CDS_expression_gff",
                    "heur_filtered_expression_gff",
                    "intron_expression_csv",
                    "gene_expression_csv",
                    "exon_expression_csv",
                    "heur_filtered_gene_expression_csv",
                    "heur_filtered_exon_expression_csv",
                    "transcript_sources" ]

Element = namedtuple( "Element", ["data_type", "sample_type", \
                                  "sample_id", "strand", "fname", "status"] )
ElementType = namedtuple( "ElementType", \
                              ["data_type", "sample_type", \
                                   "sample_id", "strand"] )


def GSTN( sample_type ):
    """Get Sample Type Name ( from sample type str )
    """
    if sample_type == '*':
        return "merged_input"
    elif sample_type == 'M':
        return "merged"
    else:
        return sample_type

GSIN = GSTN


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
            line = line.strip()
            # skip empty lines
            if len( line ) == 0: continue
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
        for dependency_fname in dependencies:
            c.execute( q, [element.fname, dependency_fname] )
        
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
        try:
            assert data_type in RnaSeqDataTypes
        except:
            print data_type
            raise
        
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
    
    def remove_element( self, fname, commit=True ):
        print "REMOVING", fname
        # first, find the elements that depend on this elements
        # and remove them. 
        deps_q = "SELECT sample_fname \n" \
               + "  FROM dependencies \n" \
               + " WHERE dependency_fname = '{0}';"
        for dependent_fname, in self.conn.execute( deps_q.format( fname ) ):
            self.remove_element( dependent_fname, commit=False )
        
        # next, remove this element from the deps table
        q = "DELETE FROM dependencies WHERE sample_fname = '{0}';".format(fname)
        self.conn.execute( q )

        # remove this element from the samples table
        q = "DELETE FROM samples WHERE fname = '{0}';".format( fname )
        self.conn.execute( q )
        
        if commit:
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
                self.remove_element( fname, commit=False )
            
            if status == 'finished' and not os.path.exists( fname ):
                self.remove_element( fname, commit=False )
        
        # we set commit to false in remove_element to speed up bulk deletes,
        # so we need to commit now.
        self.conn.commit()
        
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
                    self.chr_sizes = ChrmSizes( fp.name )
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
                print self.running_commands
                for entry in self.pending_commands:
                    print "="*80
                    print entry[0]
                    for dep in entry[0].dependencies:
                        print "\t", dep
                    
                    print 
                
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
        self.max_available_resources = available_resources
        self.available_resources = available_resources
        
    def add_process( self, cmd, min_resources, max_resources=None ):        
        if max_resources == None:
            max_resources = min_resources
        
        assert min_resources <= max_resources
        if min_resources > self.max_available_resources:
            error_str = "Can't run '%i' threads on a server with '%i' threads."
            raise ValueError, error_str % ( \
                max_resources, self.max_available_resources )

        max_resources = min( max_resources, self.max_available_resources )
        
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
        # cmd += "--merge-read-ends"
        
        dependencies.extend( ( bam_fname, ) )
        output_fnames = []
        for extension in ('rd1.plus', 'rd1.minus', 'rd2.plus', 'rd2.minus'):
            output_fnames.append( os.path.join( \
                    derived_output_dir, basename + ".%s.bedGraph" % extension ))
        
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

            et1 = ElementType( 
                "rnaseq_cov_pair1_bedgraph", 
                sample.sample_type, sample.sample_id, "+")
            et2 = ElementType( 
                "rnaseq_cov_pair1_bedgraph", 
                sample.sample_type, sample.sample_id, "-")
            et3 = ElementType( 
                "rnaseq_cov_pair2_bedgraph", 
                sample.sample_type, sample.sample_id, "+")
            et4 = ElementType( 
                "rnaseq_cov_pair2_bedgraph", 
                sample.sample_type, sample.sample_id, "-")
            
            cmd_str += " --threads={threads}"

            
            cmd = Cmd( cmd_str, [et1, et2, et3, et4], 
                       output_fnames, dependencies )
            
            process_server.add_process( cmd, Resource(1), Resource(6) )

    # first get all of the read bams
    cov_wig_fnames = defaultdict( list )
    res = elements.get_elements_from_db( "rnaseq_polyaplus_bam" )
    build_build_read_cov_cmds( res, cov_wig_fnames, rev_strand=False )
    res = elements.get_elements_from_db( "rnaseq_polyaplus_rev_bam" )
    build_build_read_cov_cmds( res, cov_wig_fnames, rev_strand=True )
    
    # add the wiggle merger process
    def add_merge_cmd( sample_type, strand, read_num, cov_bedgraph_fnames ):
        if sample_type == "*": 
            base_name = "rnaseq_cov"
        else:
            base_name =  "rnaseq_cov." + sample_type
        
        if read_num == 1:
            element_type_str = "rnaseq_cov_pair1_bedgraph"
            base_name += ".rd1"
        elif read_num == 2:
            element_type_str = "rnaseq_cov_pair2_bedgraph"
            base_name += ".rd2"
        else:
            element_type_str = "rnaseq_cov_bedgraph"
            assert read_num ==  None
        
        assert strand in '+-'
        if strand == '+':
            base_name += ".plus"
        else:
            base_name += ".minus"
        ofname = base_name + ".bedGraph"
            
        #### Get the subprocess command info together
            
        # choose with chromosome because the input wiggles are always have a chr
        # ( to make them easier to visualize in UCSC ), and so we want to choose
        # a chrm sizes with chr as well
        chrm_sizes_fname = elements.chr_sizes.get_tmp_fname( with_chr = True )
        
        # get the output filename
        ofname = os.path.join( derived_output_dir, ofname )
        
        # get the track name
        track_name = base_name
        
        cmd_str_template = "python {0} {1} --out-fname {2} "
        
        ##### get the output data types
        output_fnames = [ ofname, ]
        
        output_element_types = [ 
            ElementType( element_type_str, sample_type, "*", strand ),  \
        ]
        
        dependency_fnames = cov_bedgraph_fnames
        
        cmd_str = cmd_str_template.format( \
            MERGE_BEDGRAPHS_CMD, " ".join(cov_bedgraph_fnames), ofname )
        
        cmd = Cmd( cmd_str, output_element_types, \
                   output_fnames, dependency_fnames )
        
        process_server.add_process( cmd, Resource(1) )
    
    # build the merge sample type wiggles    
    cov_fnames = defaultdict( lambda: defaultdict(list) )
    
    res = elements.get_elements_from_db( "rnaseq_cov_pair1_bedgraph" )
    if len( res ) > 0:
        for entry in res:
            cov_fnames[(1, entry.strand)][entry.sample_type].append( entry.fname )
    
    res = elements.get_elements_from_db( "rnaseq_cov_pair2_bedgraph" )
    if len( res ) > 0:
        rd2_cov_fnames = defaultdict( list )
        for entry in res:
            cov_fnames[(2, entry.strand)][entry.sample_type].append( entry.fname )

    for (rd_num, strand), data in cov_fnames.iteritems():
        for sample_type, fnames in data.iteritems():
            add_merge_cmd( sample_type, strand, rd_num, fnames  )
        merged_fnames = list( itertools.chain( *data.values() ) )
        add_merge_cmd( "*", strand, rd_num, merged_fnames )
    
    """
    # comments this out because we no longer usee the merged RNAseq covergae,
    # but we may wish to in the future.
    
    res = elements.get_elements_from_db( "rnaseq_cov_bedgraph" )
    if len( res ) > 0:
        rd_cov_fnames = defaultdict( list )
        for entry in res:
            rd_cov_fnames[entry.sample_type].append( entry.fname )
        for sample_type, fnames in rd_cov_fnames.iteritems():
            add_merge_cmd( sample_type, fnames, None )
        add_merge_cmd("*", list( itertools.chain(*rd_cov_fnames.values() ) ), None)
    """
    
    return 

def merged_cage_wiggles( elements, process_server, derived_output_dir ):
    # create the output directory
    try:
        os.makedirs( derived_output_dir )
    except OSError:
        # assume that the directory already exists, and continue
        pass
    
    # add the wiggle merger process
    def add_merge_cmd( sample_type, strand, cov_fnames ):
        if sample_type == "*": 
            base_name = "cage_cov"
        else:
            base_name =  "cage_cov." + sample_type
        
        assert strand in '+-'
        if strand == '+':
            base_name += ".plus"
        else:
            base_name += ".minus"
        ofname = base_name + ".bedGraph"
            
        #### Get the subprocess command info together
            
        # choose with chromosome because the input wiggles are always have a chr
        # ( to make them easier to visualize in UCSC ), and so we want to choose
        # a chrm sizes with chr as well
        chrm_sizes_fname = elements.chr_sizes.get_tmp_fname( with_chr = True )
        
        # get the output filename
        ofname = os.path.join( derived_output_dir, ofname )
        
        # get the track name
        track_name = base_name
        
        cmd_str_template = "python {0} {1} --out-fname {2} "
        
        ##### get the output data types
        output_fnames = [ ofname, ]
        
        output_element_types = [ 
            ElementType( "cage_cov_wig", sample_type, "*", strand ),  \
        ]
        
        dependency_fnames = cov_fnames
        
        cmd_str = cmd_str_template.format( \
            MERGE_BEDGRAPHS_CMD, " ".join(cov_fnames), ofname )
        
        cmd = Cmd( cmd_str, output_element_types, \
                   output_fnames, dependency_fnames )
        
        process_server.add_process( cmd, Resource(1) )
    
    # build the merge sample type wiggles    
    cov_fnames = defaultdict( list )
    all_cov_fnames = { '+': [], '-': []}
    res = elements.get_elements_from_db( "cage_cov_wig" )
    if len( res ) > 0:
        for entry in res:
            cov_fnames[(entry.sample_type, entry.strand)].append( entry.fname )
            all_cov_fnames[ entry.strand ].append( entry.fname )
    
    for (sample_type, strand), fnames in cov_fnames.iteritems():
        add_merge_cmd( sample_type, strand, fnames )
    
    add_merge_cmd( "*", "+", all_cov_fnames['+'] )
    add_merge_cmd( "*", "-", all_cov_fnames['-'] )

    """
    # comments this out because we no longer usee the merged RNAseq covergae,
    # but we may wish to in the future.
    
    res = elements.get_elements_from_db( "rnaseq_cov_bedgraph" )
    if len( res ) > 0:
        rd_cov_fnames = defaultdict( list )
        for entry in res:
            rd_cov_fnames[entry.sample_type].append( entry.fname )
        for sample_type, fnames in rd_cov_fnames.iteritems():
            add_merge_cmd( sample_type, fnames, None )
        add_merge_cmd("*", list( itertools.chain(*rd_cov_fnames.values() ) ), None)
    """
    
    return 


def extract_all_junctions( elements, pserver, output_prefix ):
    try:
        os.makedirs( output_prefix )
    except OSError:
        # assume that the directory already exists, and continue
        pass
    
    def build_extract_junctions_cmd( bam_element, stranded, rev_strand ):
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
        stranded_str = '--stranded' if stranded else ""
        
        call_template = "python %s {0} --fasta {1} {2} {3} > {4} " \
            % EXTRACT_JNS_CMD
        
        call = call_template.format( \
            bam_fn, elements.genome_fname, rev_strand_str, stranded_str, output_fname ) 
        
        dependency_fnames = [ bam_fn, ]
        output_fnames = [ output_fname, ]
        output_element_type = ElementType( \
            "jns_gff", bam_element.sample_type, \
                bam_element.sample_id, bam_element.strand )
        output_element_types = [ output_element_type, ]
        
        return Cmd(call, output_element_types, output_fnames, dependency_fnames)
    
    # build all of the raw junctions
    res = elements.get_elements_from_db( "rnaseq_polyaplus_bam" )
    for e in res:
        cmd = build_extract_junctions_cmd(e, stranded=True, rev_strand = False)
        pserver.add_process( cmd, Resource(1) )

    res = elements.get_elements_from_db( "rnaseq_polyaplus_rev_bam" )
    for e in res:
        cmd = build_extract_junctions_cmd(e, stranded=True, rev_strand = True)
        pserver.add_process( cmd, Resource(1) )

    res = elements.get_elements_from_db( "rnaseq_polyaplus_unstranded_bam" )
    for e in res:
        cmd = build_extract_junctions_cmd(e, stranded=False, rev_strand = False)
        pserver.add_process( cmd, Resource(1) )

    res = elements.get_elements_from_db( "rnaseq_total_bam" )
    for e in res:
        cmd = build_extract_junctions_cmd(e, stranded=True, rev_strand=True)
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
    call = "python %s --filter --threads {threads} --fasta %s %s > %s"\
        % ( MERGE_AND_FILTER_JNS_CMD, elements.genome_fname, 
            " ".join( raw_jn_fnames ), merged_output_fname  )

    dependency_fnames = raw_jn_fnames
    output_fnames = [ merged_output_fname, ]    
    output_element_types = [ ElementType("filtered_jns_gff", "*", "*", "."), ]
    
    cmd = Cmd( call, output_element_types, output_fnames, dependency_fnames )
    pserver.add_process( cmd, Resource(pserver.max_available_resources/2), 
                              Resource( pserver.max_available_resources)  )
       
    return

def merge_sample_type_junctions( elements, pserver, output_prefix ):
    def merge_junctions( sample_type ):
        res = elements.get_elements_from_db( "filtered_jns_gff", sample_type )
        raw_jn_fnames = []
        for e in res:
            # skip the merged forms
            if e.sample_id == '*': continue
            raw_jn_fnames.append( e.fname )

        merged_output_fname = os.path.join( \
            output_prefix, "%s.filtered.jns.gff" % sample_type )
        
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
        "filtered_jns_gff", get_merged=False )

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

def get_rnaseqcov_bedgraphs( elements, sample_type, sample_id ):
    # find all of the merged rnaseq coverage files
    res = elements.get_elements_from_db( \
        "rnaseq_cov_pair1_bedgraph", sample_type, sample_id )
    rd1_bedgraph_fnames = [ i.fname for i in res ]
    if len( rd1_bedgraph_fnames ) != 2:
        print elements
        raise ValueError, "Can't find the rnaseq coverage files " \
                          + "necessary to build exons from."

    res = elements.get_elements_from_db( \
        "rnaseq_cov_pair2_bedgraph", sample_type, sample_id )
    rd2_bedgraph_fnames = [ i.fname for i in res ]
    if len( rd2_bedgraph_fnames ) != 2:
        raise ValueError, "Can't find the rnaseq coverage files " \
                          + "necessary to build exons from."

    return rd1_bedgraph_fnames + rd2_bedgraph_fnames

    pass


def build_all_exon_files( elements, pserver, output_prefix ):
    try:
        os.makedirs( output_prefix )
    except OSError:
        # assume that the directory already exists, and continue
        pass
        
    """
    python ~/grit/elements/exons/discover_exons.py 
      read_cov_bedgraphs/rnaseq_cov.plus.bedGraph 
      read_cov_bedgraphs/rnaseq_cov.minus.bedGraph 
      junctions/merged.filtered.jns.gff 
      /media/scratch/genomes/drosophila/dm3.chrom.sizes 
      --cage-wigs ../chr2L_test_region_raw_data/cage/cage_marginal.chr2L.plus.wig 
                  ../chr2L_test_region_raw_data/cage/cage_marginal.chr2L.minus.wig
      --polya-reads-gffs ../chr2L_test_region_raw_data/merged_polyA.chr2L_4605421_5481770.gff
      --out-file-prefix=test
    """
    def build_find_exons_cmds( sample_type, sample_id ):
        call_template = "python {0} {1} {2} {3} --cage-wigs {4} {5} --out-file-prefix {6}"
        
        bedgraph_fnames = get_rnaseqcov_bedgraphs( 
            elements, sample_type, sample_id )
        
        # find the merged, filtered jns file
        if USE_MERGED_JNS_FOR_EXONS:
            res = elements.get_elements_from_db( \
                "filtered_jns_gff", "*" , "*", "." )
        else:
            res = elements.get_elements_from_db( \
                "filtered_jns_gff", sample_type , "*", "." )
        if len( res ) < 1:
            raise ValueError, "Can't find the jns gff file " \
                              + "necessary to build exons."
        assert len( res ) == 1
        jns_fname = res[0].fname
        
        # find the merged CAGE signal files
        if USE_MERGED_CAGE_SIGNAL_FOR_EXONS:
            ress = elements.get_elements_from_db( \
                "cage_cov_wig", "*" , "*" )
        else:
            ress = elements.get_elements_from_db( \
                "cage_cov_wig", sample_type , "*" )
        if len( ress ) != 2:
            ress = elements.get_elements_from_db( \
                "cage_cov_wig", "*" , "*" )
        
        if len( ress ) != 2:
            raise ValueError, "Expected exactly 2 cage coverage files for " \
                + "sample %s. Found %i." % (sample_type, len( ress ))
        assert len( ress ) == 2
        cage_wig_fnames = [ res.fname for res in ress ]
        
        # find the candidate polya files
        ress = elements.get_elements_from_db( \
            "cand_polya_sites", "*" , "*" )
        assert len( ress ) in ( 0, 2 )
        
        candidate_polya_site_fnames = [ res.fname for res in ress ]
        
        chrm_sizes_fname = elements.chr_sizes.get_tmp_fname( with_chr = True )

        output_fn_prefix = os.path.join( output_prefix, \
            "discovered_exons.%s.%s" % ( sample_type, sample_id ) )
        output_fn_prefix = output_fn_prefix.replace( ".*", "" )
        
        call = call_template.format( 
            FIND_EXONS_CMD, jns_fname, chrm_sizes_fname,
            " ".join(bedgraph_fnames),   
            " ".join(cage_wig_fnames), 
            "--polya-candidate-sites "+" ".join( candidate_polya_site_fnames ),
            output_fn_prefix )
        
        call += " --threads {threads}"
        
        extensions = [ ".internal_exons.gff",    
                       ".tss_exons.gff", 
                       ".tes_exons.gff" ,
                       ".single_exon_genes.gff" ]
        
        output_fnames = [output_fn_prefix+extension for extension in extensions]
        output_element_types = [
              ElementType( "exons_gff", sample_type, sample_id, "." ),
              ElementType( "cage_tss_exons_gff", sample_type, sample_id, "." ),
              ElementType( "polya_tes_exons_gff", sample_type, sample_id, "." ),
              ElementType( "single_exon_genes_gff", sample_type, sample_id, ".")
        ]
        
        dependency_fnames = [ jns_fname ]
        dependency_fnames.extend( bedgraph_fnames )
        dependency_fnames.extend( cage_wig_fnames )
        
        cmd = Cmd(call, output_element_types, output_fnames, dependency_fnames)
        pserver.add_process( cmd, Resource(1), Resource(6) )

    # build the merged sample exons
    build_find_exons_cmds( "*", "*" )
    
    if BUILD_SAMPLE_SPECIFIC_EXONS:
        element_types = [ "filtered_jns_gff", "rnaseq_cov_pair1_bedgraph", 
                         "rnaseq_cov_pair2_bedgraph"]
        sample_types_and_ids = elements.get_distinct_element_types_and_ids( \
            element_types, get_merged=True )
        
        # buidl the build eoxn cmds
        for sample_type, sample_id in sample_types_and_ids:
            if USE_MERGED_RNASEQ_FOR_EXONS and sample_id == '*':
                build_find_exons_cmds( sample_type, sample_id )
            elif not USE_MERGED_RNASEQ_FOR_EXONS and sample_id != '*':
                build_find_exons_cmds( sample_type, sample_id )
    
    
    return

def build_transcripts( elements, pserver, output_prefix, use_TF_elements=False ):
    """
    build_transcripts.py [-h] --exons EXONS [EXONS ...] 
                              --tss TSS [TSS ...]
                              --tes TES [TES ...] 
                              --junctions JUNCTIONS [JUNCTIONS ...] 
                              [ --single-exon-genes ]
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
        exon_input_type = 'TF_exons_gff'
        jn_input_type = 'TF_jns_gff'
        transcript_output_type = "TF_filt_element_transcripts_gtf"
    else:
        exon_input_type = 'exons_gff'
        jn_input_type = 'filtered_jns_gff'
        transcript_output_type = "transcripts_gtf"
        
    def build_transcripts_for_sample( sample_type, sample_id ):
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
        if not USE_MERGED_EXONS_FOR_TRANSCRIPTS and BUILD_SAMPLE_SPECIFIC_EXONS:
            ress = elements.get_elements_from_db(\
                "cage_tss_exons_gff", sample_type, sample_id)
        else:
            ress = elements.get_elements_from_db(\
                "cage_tss_exons_gff", "*", "*" )

        assert len( ress ) == 1
        
        tss_exons_fnames = [ ress[0].fname, ]
        ress = elements.get_elements_from_db( "cdna_tss_gff" )
        tss_exons_fnames.extend( e.fname for e in ress )
        
        # get tes exons
        if not USE_MERGED_EXONS_FOR_TRANSCRIPTS and BUILD_SAMPLE_SPECIFIC_EXONS:
            ress = elements.get_elements_from_db( \
                "polya_tes_exons_gff", sample_type, sample_id  )
        else:
            ress = elements.get_elements_from_db( \
                "polya_tes_exons_gff", "*", "*"  )

        assert len( ress ) == 1
        tes_exons_fnames = [ ress[0].fname, ]
        ress = elements.get_elements_from_db( "cdna_tes_gff" )
        tes_exons_fnames.extend( e.fname for e in ress )

        # get single exon genes
        if not USE_MERGED_EXONS_FOR_TRANSCRIPTS and BUILD_SAMPLE_SPECIFIC_EXONS:
            ress = elements.get_elements_from_db( \
                "single_exon_genes_gff", sample_type, sample_id  )
        else:
            ress = elements.get_elements_from_db( \
                "single_exon_genes_gff", "*", "*"  )

        assert len( ress ) == 1
        single_exon_genes_fnames = [ ress[0].fname, ]

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
                     + "--single-exon-genes {4} --out-fname {5} --log-fname {6}"
        call = call % BUILD_TRANSCRIPTS_CMD
        
        call = call.format( exons_fname, 
                            " ".join( tss_exons_fnames ),
                            " ".join( tes_exons_fnames ),
                            jns_fname,
                            " ".join( single_exon_genes_fnames ),
                            op_fname, log_fname  )
        call +=  " --threads={threads}"
        
        op_element_types = [ \
            ElementType( transcript_output_type, sample_type, sample_id, "." ),]
        op_fnames = [ op_fname,]

        dependencies = [ exons_fname, jns_fname ]
        dependencies.extend( tss_exons_fnames )
        dependencies.extend( tes_exons_fnames )
        
        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        
        max_res = min(Resource(8), Resource(pserver.max_available_resources))
        pserver.add_process( cmd, Resource(1), max_res )
        return

    ress = elements.get_elements_from_db( exon_input_type )
    for e in ress[:]:
        if USE_MERGED_RNASEQ_FOR_EXONS and e.sample_id != '*': 
            continue
        elif not USE_MERGED_RNASEQ_FOR_EXONS and e.sample_id == '*':
            continue
        
        build_transcripts_for_sample( e.sample_type, e.sample_id  )
    
    build_transcripts_for_sample( "*", "*"  )
    
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
            assert False
            suffix = ".TFfilt_transcripts.gtf"
        elif op_e_type == 'transcripts_gtf':
            suffix = ".transcripts.gtf"
        else:
            assert False
        
        #if sample_type == "*": 
        #    raise ValueError, "Can't merge merged transcripts."

        sample_type_name = "merged" if sample_type == "*" else sample_type
        op_fname = os.path.join(output_prefix, \
                                    sample_type_name + suffix)

        sources_fname = os.path.join(output_prefix, \
                                    sample_type_name + ".sources")
        
        call = "python %s %s --sources-fname %s --threads={threads} > %s" \
            % ( MERGE_TRANSCRIPTS_CMD, " ".join( input_fnames ), \
                    sources_fname, op_fname )
        
        op_element_types = [ \
            ElementType( op_e_type, sample_type, "M", "." ),
            ElementType( "transcript_sources", sample_type, "M", "." ) ]
        op_fnames = [ op_fname, sources_fname ]
        dependencies = input_fnames
        
        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        print cmd
        
        res = min(pserver.max_available_resources, 16, len(input_fnames) )
        pserver.add_process( cmd, Resource( res ) )
    
    # add the merge all command first, since it takes the 
    # longest to complete
    in_es = elements.get_elements_from_db( input_e_type )
    input_fnames = [ e.fname for e in in_es \
                     if e.sample_id in "*" and e.sample_type != '*']
    if USE_MERGED_INPUT:
        merged_input_fname = elements.get_elements_from_db( \
            input_e_type, "*", "*" )[0].fname
        input_fnames.append( merged_input_fname )
    
    add_merge_transcripts_sample_cmd( "*", input_fnames )

    """
    # get the merged input transcripts    
    if USE_MERGED_INPUT:
        merged_input_fname = elements.get_elements_from_db( \
            input_e_type, "*", "*" )[0].fname
        input_fnames.append( merged_input_fname )


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
    """
    
    return

def estimate_fl_dists( elements, pserver, output_prefix ):
    try:
        os.makedirs( output_prefix )
    except OSError:
        # assume that the directory already exists, and continue
        pass
    
    """
     python ~/Desktop/grit/sparsify/frag_len.py 
     ~/Desktop/grit_test_data/chr4_transcriptome/exons/
           discovered_exons.AdVirginF_Ecl_1day_Heads.sample_1.tes_exons.gff 
     /home/nboley/Desktop/grit_test_data/chr4_test_data/raw_data/bams/
           AdVirginF_Ecl_1day_Heads.497_BS234.chr4.bam  
     AdVirginF_Ecl_1day_Heads.sample_1.fldist 
     --analyze
    """
    def build_fl_dist_cmd( sample_type, sample_id ):
        # get the bam file for this sample_type, id combo
        ress = elements.get_elements_from_db( \
            "rnaseq_polyaplus_bam", sample_type, sample_id )
        if len( ress ) == 0:
            ress = elements.get_elements_from_db( \
                "rnaseq_polyaplus_rev_bam", sample_type, sample_id )
        assert len( ress ) == 1
        bam_fname = ress[0].fname
        
        # get the exons file ( we always use distal exons, because reads that 
        # fall entirely within them are necessarily full length
        ress = elements.get_elements_from_db( \
            "cage_tss_exons_gff", sample_type, "*" )
        assert len( ress ) == 1
        exons_fname = ress[0].fname
        
        # get the output file name
        op_fname = os.path.join( output_prefix, sample_type + "." \
                                     + sample_id + ".fldist.obj" )
        
        cmd_str_template = "python %s --analyze {0} {1} {2} " % BUILD_FL_DIST_CMD
        call = cmd_str_template.format( exons_fname, bam_fname, op_fname )
        
        op_element_types = [ \
            ElementType( "fl_dist", sample_type, sample_id, "." ),]
        op_fnames = [ op_fname,]

        dependencies = [ bam_fname, exons_fname ]
        
        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        pserver.add_process( cmd, Resource(1) )
    

    sample_types_and_ids = set()
    sample_types_and_ids.update( elements.get_distinct_element_types_and_ids( \
                [ "rnaseq_polyaplus_bam", ], False ) )
    sample_types_and_ids.update( elements.get_distinct_element_types_and_ids( \
                [ "rnaseq_polyaplus_rev_bam", ], False ) )

    for sample_type, sample_id in sorted( sample_types_and_ids ):
        build_fl_dist_cmd( sample_type, sample_id )
    
    return


def sparsify_transcripts( elements, pserver, output_prefix ):
    try:
        os.makedirs( output_prefix )
    except OSError:
        # assume that the directory already exists, and continue
        pass
    
    """
    python ~/Desktop/grit/sparsify/sparsify_transcripts.py 
    tmp.transcripts.sparse.gtf                            # ofname
    transcripts/merged_input.merged_input.transcripts.gtf #candidate transcripts
    ../chr4_test_data/raw_data/bams/L3_ImaginalDiscs.56*.bam 
    --fl-dists fl_dists/*.obj 
    --threads=3 
    -m
    """
    
    # get all of the bam files
    ress = elements.get_elements_from_db( "rnaseq_polyaplus_bam" )
    bam_fnames = [ res.fname for res in ress ]
    
    ress = elements.get_elements_from_db( "fl_dist" )
    fldist_fnames = [ res.fname for res in ress ]
    
    ress = elements.get_elements_from_db( "transcripts_gtf", "*", "M" )
    assert( len(ress) == 1 )
    merged_transcript_fname = ress[0].fname
    
    op_fname = os.path.join( output_prefix, "merged.sparse.transcripts.gtf" )

    cmd_str_template  = "python %s %s %s %s --fl-dists %s -m "   \
        % ( SPARSIFY_TRANSCRIPTS_CMD,                            \
            op_fname, merged_transcript_fname,                    \
            " ".join(bam_fnames), " ".join( fldist_fnames)  )
    cmd_str_template += "--threads={threads}"
    call = cmd_str_template
    
    op_element_types = [ \
        ElementType( "sparse_transcripts_gtf", "M", "M", "." ),]
    op_fnames = [ op_fname,]
    
    dependencies = [ merged_transcript_fname, ]
    dependencies.extend( bam_fnames )
    dependencies.extend( fldist_fnames )

    cmd = Cmd( call, op_element_types, op_fnames, dependencies )
    pserver.add_process( cmd,  Resource(pserver.max_available_resources) )
        
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
    
    return


def call_orfs( elements, pserver, output_prefix ):
    """Find orfs for each transript file.

    Since ORF finding can be slow, and doing it sample by sample is very 
    redundant, we call orfs on the merged samples and then split them
    using the merged samples file.
    """
    
    try:
        os.makedirs( output_prefix )
    except OSError:
        # assume that the directory already exists, and continue
        pass
    
    """
    time python ~/grit/proteomics/ORF_finder.py \
        ./transcripts/merged_input.merged_input.transcripts.gtf  \
        /media/scratch/genomes/drosophila/BDGP_5/all_fa/fly.fa \
        --threads=50 \
        --output-filename output.cds.gtf
    
    OUTPUTS:
    output.cds.gtf

    python ~/Desktop/grit/elements/transcripts/split_merged_transcripts.py \
        transcripts/merged.transcripts.gtf 
        transcripts/merged.sources 
        --suffix .cds.gtf
        --output-dir output_prefix
    
    OUTPUTS:
    *.CDS.gtf
    
    
    """
    
    def call_orfs( sample_type, sample_id  ):
        # get the merged transcripts from the database
        res = elements.get_elements_from_db(
            "transcripts_gtf", sample_type, sample_id)
        assert len( res ) == 1
        merged_trans_fname = res[0].fname

        # get the fasta file from the DB
        fasta_fn = elements.genome_fname

        # get the CODING annotation file 
        res = elements.get_elements_from_db( "annotation_gtf" )
        assert len( res ) == 1
        coding_ann_fname = res[0].fname

        # just append CDS.gtf to the transcripts gtf, for the CDS annotated GTF
        basename = os.path.basename(merged_trans_fname)
        if basename.endswith( '.gtf' ): basename = basename[:-4]
        op_fname = os.path.join( output_prefix, basename + ".CDS.gtf")
        
        call = "python %s {0} {1} --only-longest-orf --output-filename {2}" % FIND_ORFS_CMD
        call = call.format( merged_trans_fname, fasta_fn, op_fname )
        call += " --threads {threads}"
        
        op_element_types = [ 
            ElementType( "CDS_transcripts_gtf", sample_type, sample_id, "." ), ]
        op_fnames = [ op_fname, ]
        dependencies = [ merged_trans_fname, fasta_fn ]
        
        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        pserver.add_process( cmd, 
                             Resource(pserver.max_available_resources-2), 
                             Resource(pserver.max_available_resources) )


    # TODO - add split ORFs cmd
    def split_orf_results_by_sample( elements, pserver, output_prefix ):
        try:
            os.makedirs( output_prefix )
        except OSError:
            # assume that the directory already exists, and continue
            pass

        """pass orf_transcripts to split_orf_results_by_sample.py and produce 
        one orf results per merge_input_fn 

        python split_orf_results_by_sample.py 
             tmp.annotated.gtf 
             merged.sources 
             --out-prefix ./proteomics/tmp.annotated.gtf
        """
        orf_gtf_e = elements.get_elements_from_db( 
            "coding_transcripts_gtf", "M", "*")
        assert len( orf_gtf_e ) == 1
        orf_gtf_fn = orf_gtf_e[0].fname

        sources_e = elements.get_elements_from_db("transcript_sources", "M", "*")
        assert len( sources_e ) == 1
        sources_fn = sources_e[0].fname

        # find the transcript filenames, by source. 
        in_es = elements.get_elements_from_db( "transcripts_gtf" )
        transcript_fnames_and_sample_types = [ 
            (e.sample_type, e.fname) for e in in_es 
            if e.sample_id == "*" and e.sample_type not in 'M*' ]

        # build the output
        op_fns = []
        op_element_types = []
        for sample_type, sample_fn in transcript_fnames_and_sample_types:
            op_fns.append( sample_fn + '.annotated.gtf' )
            op_element_types.append( ElementType( "coding_transcripts_gtf", 
                                                  sample_type, "*", "." ) )

        dependencies = [ orf_gtf_fn, sources_fn ]

        call = "python {0} {1} {2} --out-prefix {3}".format( 
            SPLIT_ORFS_CMD, orf_gtf_fn, sources_fn, output_prefix + '/' )

        cmd = Cmd( call, op_element_types, op_fns, dependencies )
        pserver.add_process( cmd, Resource( 1 ) )

        return
    
    # call on merged input transcripts
    call_orfs("*", "*")
    # call on merged transcripts
    call_orfs("*", "M")
    
    return

def produce_final_annotation( elements, pserver, output_prefix ):
    """Produce a final annotation.

    This means running everything through the heuristic transcript 
    filtering scripts, and then running the renaming script. The 
    renaming script also does some filtering, particularly for gene 
    merges, and then re-adds any transcript models from the annotation
    that we ended up missing.
    """
    try:
        os.makedirs( output_prefix )
        os.makedirs( os.path.join( output_prefix, "heuristic_filtered" ) )
    except OSError:
        # assume that the directory already exists, and continue
        pass    
    """
    python ~/Desktop/grit/sparsify/filter_transcript_set.py 
        CDS_transcripts/merged_input.merged_input.transcripts.CDS.gtf 
        ../chr4_test_data/raw_data/chr4_flybase.gtf  
        ../chr4_test_data/raw_data/chr4.polya_cdna.gtf
    
    """
    def add_heur_filter_transcripts_cmd( sample_type, sample_id  ):
        assert False
        # get the merged transcripts from the database
        res = elements.get_elements_from_db(
            "CDS_transcripts_gtf", sample_type, sample_id)
        assert len( res ) == 1
        trans_fname = res[0].fname

        # get the CODING annotation file 
        res = elements.get_elements_from_db( "annotation_gtf" )
        assert len( res ) == 1
        ann_fname = res[0].fname

        # get the full length cdnas file 
        res = elements.get_elements_from_db( "cdna_gtf" )
        assert len( res ) == 1
        cdnas_fname = res[0].fname
        
        op_fname = os.path.join( output_prefix, "heuristic_filtered",
                                 "%s.%s.transcripts.cds.heur_filt.gtf" %
                                 ( GSTN(sample_type), GSTN(sample_id) ) )
        
        call = "python %s {0} {1} {2} > {3}" % HEURISTIC_FILTER_TRANS_CMD
        call = call.format( trans_fname, ann_fname, cdnas_fname, op_fname )
        
        op_element_types = [ 
            ElementType( "heur_filtered_transcripts_gtf", 
                         sample_type, sample_id, "." ), ]
        op_fnames = [ op_fname, ]
        dependencies = [ trans_fname, ann_fname, cdnas_fname ]
        
        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        pserver.add_process( cmd, Resource(1) )

    """
        python ~/Desktop/grit/utilities/group_transcripts_by_reference.py 
        final_ann/heuristic_filtered/MI.MI.transcripts.cds.heur_filt.gtf 
        ../chr4_test_data/raw_data/chr4_flybase.gtf 
        --output-file tmp.renamed.gtf 
        --good-gene-merges ../chr4_test_data/good_gene_merges.txt 
        --tss-exons exons/discovered_exons.tss_exons.gff

    """
    def add_rename_transcripts_cmd( 
            sample_type, sample_id, use_heuristic_filtered  ):
        
        if use_heuristic_filtered:
            ip_element_type = "heur_filtered_transcripts_gtf"
            op_element_type_prefix = "heur_filtered"
        else:
            ip_element_type = "CDS_transcripts_gtf"
            op_element_type_prefix = "CDS"

        # get the merged transcripts from the database
        res = elements.get_elements_from_db( 
            ip_element_type, sample_type, sample_id)
        assert len( res ) == 1
        trans_fname = res[0].fname
        
        # get the CODING annotation file 
        res = elements.get_elements_from_db( "annotation_gtf" )
        assert len( res ) == 1
        ann_fname = res[0].fname
        
        op_fname = os.path.join( output_prefix, "%s.%s.%s.renamed.gtf" %(
            GSTN(sample_type), GSTN(sample_id), op_element_type_prefix ) )
        
        # get the tss exons
        res = elements.get_elements_from_db( "cage_tss_exons_gff", sample_type )
        assert len( res ) == 1
        tss_exons_fname = res[0].fname

        res = elements.get_elements_from_db( "gene_merges_white_list" )
        assert len( res ) in ( 0, 1 )
        good_gene_merges_fname = res[0].fname if len(res) == 1 else None
        
        call  = "python {0} {1} {2}".format( 
            RENAME_TRANS_CMD, trans_fname, ann_fname  )
        call += " --output-file {0}".format( op_fname )
        dependencies = [ trans_fname, ann_fname ]

        if good_gene_merges_fname != None:
            call += " --good-gene-merges {0}".format( good_gene_merges_fname )
            dependencies.append( good_gene_merges_fname )
        
        if tss_exons_fname != None:
            call += " --tss-exons {0}".format( tss_exons_fname )
            dependencies.append( tss_exons_fname )
        
        op_element_types = [ 
            ElementType( op_element_type_prefix + "_renamed_transcripts_gtf", 
                         sample_type, sample_id, "." ), ]

        op_fnames = [ op_fname,  ]
        
        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        pserver.add_process( cmd, Resource(1) )

    if USE_HEURISTIC_FILTERING:
        add_heur_filter_transcripts_cmd( "*", "M" )
        add_rename_transcripts_cmd( "*", "M", True )
        add_heur_filter_transcripts_cmd( "*", "*" )
        add_rename_transcripts_cmd( "*", "*", True )

    
    add_rename_transcripts_cmd( "*", "M", False )
    add_rename_transcripts_cmd( "*", "*", False )
    
    return

def calc_expression_scores( elements, pserver, output_prefix ):
    try:
        os.makedirs( output_prefix )
        os.makedirs( os.path.join( output_prefix, "expression_gffs" ) )
    except OSError:
        # assume that the directory already exists, and continue
        pass

    """
        python ~/Desktop/grit/expression/estimate_exon_expression.py 
        final_ann/merged_input.merged_input.CDS.missing.gtf 
        /home/nboley/Desktop/grit_test_data/genomes/drosophila/dm3.chrom.sizes 
        read_cov_bedgraphs/L3_ImaginalDiscs.566_BS303.*
        
    """
    def add_calc_exon_expression_cmd(sample_type, sample_id, 
                                     t_element_type, t_sample_type, t_sample_id,
                                     op_type, op_e_type):
        bedgraph_fnames = get_rnaseqcov_bedgraphs( 
            elements, sample_type, sample_id )
        
        res = elements.get_elements_from_db( 
            t_element_type, t_sample_type, t_sample_id )
        assert len( res ) == 1
        trans_fname = res[0].fname
        
        chrm_sizes_fname = elements.chr_sizes.get_tmp_fname( with_chr = True )
        
        op_fname = os.path.join( output_prefix, "expression_gffs", 
                                 "%s.%s.%s.expression.gff" % ( 
                GSTN(sample_type), GSTN(sample_id), op_type ) )
                
        cmd_str_template  = "python %s {0} {1} {2} > {3}" % EST_EXON_EXP_CMD
        call = cmd_str_template.format( trans_fname, chrm_sizes_fname, 
                                        " ".join( bedgraph_fnames ), op_fname )
        
        op_fnames = [ op_fname, ]
        op_element_types = [ \
            ElementType( op_e_type, sample_type, sample_id, "." ),]
        
        dependencies = [ trans_fname, ]
        dependencies.extend( bedgraph_fnames )
        
        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        pserver.add_process( cmd,  Resource(1) )
    
    def add_build_csv_cmd( element_type, use_heuristic_filtering=True ):
        """
        python ~/Desktop/grit/utilities/make_csv_from_gtfs.py 
        expression/*.CDS.* 
        junctions/*.* 
        --intron-output-filename 
        intron.tmp.csv
        """
        assert element_type in ['intron', 'exon', 'gene']
        if element_type == 'intron':
            input_element_type = 'filtered_jns_gff'
            op_fname = os.path.join( output_prefix, "intron_expression.csv" )
            output_flag = "--intron-output-filename "
            op_fname_prefix = ""
        else:
            if use_heuristic_filtering:
                input_element_type = 'heur_filtered_expression_gff'
                op_fname_prefix = "heur_filtered_"
            else:
                input_element_type = 'CDS_expression_gff'
                op_fname_prefix = "" 

            if element_type == 'exon':
                op_fname = os.path.join( 
                    output_prefix, op_fname_prefix + "exon_expression.csv" )
                output_flag = "--exon-output-filename "
            elif element_type == 'gene':
                op_fname = os.path.join( 
                    output_prefix, op_fname_prefix + "gene_expression.csv" )
                output_flag = "--gene-output-filename "
            else:
                raise ValueError, "Unrecognized element type for " \
                    + "quantification '%s'." % element_type
        
        input_fnames = [ element.fname for element in 
                         elements.get_elements_from_db(input_element_type) ]
        assert len( input_fnames ) > 0
        cmd = "python %s {0} {1} {2}" % AGGREGATE_EXP_CMD
        call = cmd.format( " ".join(input_fnames), output_flag, op_fname )
        
        op_fnames = [ op_fname, ]
        op_element_types = [ \
            ElementType( op_fname_prefix + element_type + "_expression_csv", 
                         "M", "M", "." ),]
        
        dependencies = []
        dependencies.extend( input_fnames )
        
        cmd = Cmd( call, op_element_types, op_fnames, dependencies )
        pserver.add_process( cmd,  Resource(1) )
        
        return
    
    sample_types_and_ids = elements.get_distinct_element_types_and_ids( 
        ["rnaseq_cov_pair1_bedgraph", "rnaseq_cov_pair2_bedgraph"], 
        get_merged=True )
    
    for sample_type, sample_id in sample_types_and_ids:
        add_calc_exon_expression_cmd(sample_type, sample_id, 
                                     "CDS_renamed_transcripts_gtf", "*", "*", 
                                     "CDS", "CDS_expression_gff")
        if USE_HEURISTIC_FILTERING:
            add_calc_exon_expression_cmd(
                sample_type, sample_id, 
                "heur_filtered_renamed_transcripts_gtf", 
                "*", "*", 
                "heur_filtered",
                "heur_filtered_expression_gff")
    
    # add build csv command
    add_build_csv_cmd( "intron" )
    add_build_csv_cmd( "exon", False )
    add_build_csv_cmd( "gene", False )
    
    if USE_HEURISTIC_FILTERING:
        add_build_csv_cmd( "exon", True )
        add_build_csv_cmd( "gene", True )
        
    return


def main():
    elements_fname = sys.argv[1]
    num_threads = int( sys.argv[2] )
    base_dir = sys.argv[3]
    
    with open( elements_fname ) as fp:
        elements = Elements( fp )
    
    pserver = ProcessServer( elements, Resource( num_threads )  )

    merged_cage_wiggles( elements, pserver, base_dir + "cage_wiggles/" )
    
    build_all_cov_wiggles( elements, pserver, base_dir + "read_cov_bedgraphs" )
    
    extract_all_junctions( elements, pserver, base_dir + "junctions/" )
    build_high_quality_junctions( elements, pserver, base_dir + "junctions/" )
    intersect_all_junctions( elements, pserver, base_dir + "junctions/" )
    merge_sample_type_junctions( elements, pserver, base_dir + "junctions/" )
        
    build_all_exon_files( elements, pserver, base_dir + "exons" )
    estimate_fl_dists( elements, pserver, base_dir + "fl_dists" )
    
    build_transcripts( elements, pserver, base_dir + "transcripts" )

    merge_transcripts( elements, pserver, base_dir + "transcripts" )
    
    call_orfs( elements, pserver, base_dir + "CDS_transcripts" )

    produce_final_annotation( elements, pserver, base_dir + "final_ann" )

    calc_expression_scores( elements, pserver, base_dir + "expression" )
    
    run_all_slide_compares( elements, pserver, base_dir + "stats" )
    return
    #sparsify_transcripts( elements, pserver, base_dir + "transcripts" )
    
    pserver.process_queue()                
    return

main()

