#!/usr/bin/python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep


################################################################
def parse_options( argv, tag, default):
    value = default
    if argv.count( "-"+tag ):
        pos = argv.index( "-"+tag )
        try:
            if isinstance( default, int ):
                value = int( argv[ pos + 1 ] )
            elif isinstance( default, float ):
                value = float( argv[ pos + 1 ] )
            else:
                value = argv[ pos + 1 ]
        except:
            value = 1
    return value

fasta_file = parse_options( argv, "fasta", "1shf.fasta" )
assert( exists( fasta_file ) )
sequence = open( fasta_file  ).readlines()[1][:-1]

MIN_RES = parse_options( argv, "min_res", 1 )
MAX_RES = parse_options( argv, "max_res", len( sequence ) )
ZIGZAG = parse_options( argv, "zigzag", 0 )
N_SAMPLE = parse_options( argv, "n_sample", 18 )
FINAL_NUMBER = parse_options( argv, "final_number", 40 )
SCORE_WEIGHTS = parse_options( argv, "weights", "score12.wts" )
PACK_WEIGHTS = parse_options( argv, "pack_weights", "pack.wts" )
N_MINIMIZE = parse_options( argv, "n_minimize", 100 )
FILTER_RMSD = parse_options( argv, "filter_rmsd", 999.999 )
CLUSTER_RADIUS = parse_options( argv, "cluster_radius", 2.0 )
filter_native_big_bins = parse_options( argv, "filter_native_big_bins", 0 )
native_pdb = parse_options( argv, "native", "1shf.pdb" )
cst_file = parse_options( argv, "cst_file", "" )

assert( exists( SCORE_WEIGHTS ) )
assert( exists( PACK_WEIGHTS ) )
assert( exists( native_pdb ) ) # Get rid of this later...

###############################################################
# Where's the executable?
###############################################################
HOMEDIR = expanduser('~')

EXE = HOMEDIR+'/src/mini/bin/stepwise_test.macosgccrelease'
if not( exists( EXE )):
    EXE = HOMEDIR+'/src/mini/bin/stepwise_test.linuxgccrelease'
assert( exists( EXE ) )

DB = HOMEDIR+'/minirosetta_database'
assert( exists( DB ) )

PYDIR = HOMEDIR+'/python'
assert( exists( PYDIR ) )


PRE_PROCESS_SETUP_SCRIPT = PYDIR+"/stepwise_pre_process_setup_dirs.py"
assert( exists( PRE_PROCESS_SETUP_SCRIPT ) )

POST_PROCESS_FILTER_SCRIPT = PYDIR+"/stepwise_post_process_combine_and_filter_outfiles.py"
assert( exists( POST_PROCESS_FILTER_SCRIPT ) )

POST_PROCESS_CLUSTER_SCRIPT = PYDIR+"/stepwise_post_process_cluster.py"
assert( exists( POST_PROCESS_CLUSTER_SCRIPT ) )


fid_dag = open( "full_build.dag", 'w' )
fid_dag.write("DOT dag.dot\n")

###############################################################
# MAIN LOOP
###############################################################

# Loop over fragment lengths.
# Here make them in chunks of two to simplify this first calculation.

BLOCK_SIZE = 1

system( 'mkdir -p CONDOR' )
def make_condor_submit_file( condor_submit_file, arguments, queue_number ):
    fid = open( condor_submit_file, 'w' )
    fid.write('+TGProject = TG-MCB090153\n')
    fid.write('universe = vanilla\n')
    fid.write('executable = %s\n' % EXE )
    fid.write('arguments = %s\n' % arguments)

    job_tag = basename(condor_submit_file).replace('.condor','')

    subdir = 'CONDOR/'+job_tag
    if not exists( subdir ): system( 'mkdir -p '+subdir )

    fid.write('output = CONDOR/%s/$(Process).out\n' % job_tag )
    fid.write('log = CONDOR/%s.log\n' % job_tag )
    fid.write('error = CONDOR/%s/$(Process).err\n' % job_tag)
    fid.write('notification = never\n')
    fid.write('Queue %d\n' % queue_number )
    fid.close()


def  parse_cst_file( cst_file, i, j, cst_file_for_step):
    assert( exists( cst_file ) )
    lines = open( cst_file ).readlines()
    fid = open( cst_file_for_step, 'w' )
    fid.write( '[ atompairs ]\n' )
    num_cst = 0
    for line in lines[1:]:
        if len( line ) > 10:
            cols = string.split( line )
            atom_name1 = cols[0]
            res_num1 = int( cols[1] )
            atom_name2 = cols[2]
            res_num2 = int( cols[3] )
            if ( res_num1 in range(i,j+1)  ) and \
               ( res_num2 in range(i,j+1)  ):
                res_num1 -= (i-1)
                res_num2 -= (i-1)
                fid.write('%s %d %s %d %s\n' % (atom_name1, res_num1, atom_name2, res_num2, string.join( cols[4:]) ) )
                num_cst += 1
    fid.close()
    return num_cst

def parse_fasta_file( fasta_file, i, j, fasta_for_step):
    lines = open( fasta_file ).readlines()
    fid = open( fasta_for_step, 'w' )
    fid.write( lines[0] )
    fid.write( lines[1][(i-1):j] + '\n' )
    fid.close()

all_job_tags = []

for L in range( 2, len(sequence)/BLOCK_SIZE + 1 ):

    chunk_length = BLOCK_SIZE * L;
    num_chunks = ( len( sequence) - chunk_length) / BLOCK_SIZE + 1

    for k in range( 1, num_chunks+1 ) :
        i = BLOCK_SIZE * ( k - 1 ) + 1
        j = i + chunk_length - 1

        if ( i < MIN_RES or j > MAX_RES ): continue

        #ZIGZAG!!
        if ( ZIGZAG and abs( ( i - MIN_RES ) - ( MAX_RES - j ) ) > 1 ) : continue

        print 'DO_CHUNK',i,j

        # Native PDB.
        prefix = 'region_%d_%d_' % (i,j)

        # This jobs is maybe already done...
        outfile_cluster = prefix+'sample.cluster.out'
        #if exists( outfile_cluster ):
        #    continue

        native_pdb_for_step = prefix + native_pdb
        if not exists( native_pdb_for_step ):
            command = 'pdbslice.py %s %d %d %s' % ( native_pdb, i, j, prefix )
            print( command )
            system( command )

        ###########################################
        # OUTPUT DIRECTORY
        outdir = 'REGION_%d_%d' % (i,j)
        if not( exists( outdir) ):
            system( 'mkdir -p ' + outdir )

        fasta_for_step = outdir + '/' + prefix + fasta_file
        if not exists( fasta_for_step ):
            parse_fasta_file( fasta_file, i, j, fasta_for_step)


        # BASIC COMMAND
        extraflags = '-extrachi_cutoff 0 -ex1 -ex2 -score:weights %s -pack_weights %s' % (SCORE_WEIGHTS, PACK_WEIGHTS )
        args = ' -out:file:silent_struct_type binary -database %s  -rebuild -native %s -fasta %s -n_sample %d -n_minimize %d -minimize  -fullatom %s  -filter_rmsd %8.3f  ' % ( DB, native_pdb_for_step, fasta_for_step, N_SAMPLE, N_MINIMIZE, extraflags, FILTER_RMSD )

        if filter_native_big_bins:  args+= " -filter_native_big_bins "

        if len( cst_file ) > 0:
            cst_file_for_step = outdir + '/' + prefix + cst_file
            num_cst = parse_cst_file( cst_file, i, j, cst_file_for_step)
            if (num_cst > 0 ): args += ' -cst_file %s ' % cst_file_for_step


        overall_job_tag = 'REGION_%d_%d' % (i,j)

        ###########################################
        # DO THE JOBS
        start_regions = []
        if ( i + BLOCK_SIZE < j ):
            start_regions.append( [i, j - BLOCK_SIZE] )
            start_regions.append( [i + BLOCK_SIZE,j] )

        job_tags = []
        combine_files = []


        if len( start_regions ) == 0:
            ##########################################
            # START FROM SCRATCH
            ##########################################
            outfiles = []
            newdir = outdir+'/START_FROM_SCRATCH'
            if not exists( newdir ): system( 'mkdir -p '+newdir )
            outfile = newdir + '/' + prefix + 'sample.out'
            start_tag = ' -start_from_scratch '

            job_tag = 'REGION_%d_%d_START_FROM_SCRATCH' % (i,j)
            condor_submit_file = 'CONDOR/%s.condor' %  job_tag.lower()
            fid_dag.write('\nJOB %s %s\n' % (job_tag,condor_submit_file) )
            args2 = '%s -out:file:silent %s %s' % (args, outfile, start_tag )
            make_condor_submit_file( condor_submit_file, args2, 1 )
            fid_dag.write('SCRIPT POST %s %s %s\n' % (job_tag, POST_PROCESS_FILTER_SCRIPT,newdir) )

            job_tags.append( job_tag )
            combine_files.append( '%s/start_from_scratch_sample_minimize.low4000.out' % outdir )
        else:
            ##########################################
            # APPEND OR PREPEND TO PREVIOUS PDB
            ##########################################
            for start_region in start_regions:
                i_prev = start_region[0]
                j_prev = start_region[1]

                dir_prev = 'REGION_%d_%d' % (i_prev, j_prev )
                prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)

                if ZIGZAG and (prev_job_tag not in all_job_tags): #Note previous job may have been accomplished in a prior run -- not in the current DAG.
                    continue

                pdbfile = '%s/region_%d_%d_sample.cluster.out.$(Process).pdb' % (dir_prev,i_prev,j_prev)
                tag = basename( pdbfile ).replace( '.pdb' ,'' )
                newdir = outdir+'/START_FROM_'+tag.upper()
                outfile = newdir + '/' + prefix + 'sample.out'

                start_tag = ' -s %s ' % pdbfile
                args2 = '%s -out:file:silent %s %s' % (args, outfile, start_tag )

                job_tag = 'REGION_%d_%d_START_FROM_REGION_%d_%d' % (i,j,i_prev,j_prev)
                condor_submit_file = 'CONDOR/%s.condor' %  job_tag.lower()
                fid_dag.write('\nJOB %s %s\n' % (job_tag, condor_submit_file) )

                make_condor_submit_file( condor_submit_file, args2, FINAL_NUMBER )

                if prev_job_tag in all_job_tags: #Note previous job may have been accomplished in a prior run -- not in the current DAG.
                    fid_dag.write('PARENT %s  CHILD %s\n' % (prev_job_tag, job_tag) )

                fid_dag.write('SCRIPT PRE %s   %s %s %s %s\n' % (job_tag, PRE_PROCESS_SETUP_SCRIPT,outdir,dir_prev,condor_submit_file) )
                fid_dag.write('SCRIPT POST %s %s %s/START_FROM_REGION_%d_%d\n' % (job_tag, POST_PROCESS_FILTER_SCRIPT,outdir,i_prev,j_prev ) )

                job_tags.append( job_tag )
                combine_files.append( '%s/start_from_region_%d_%d_sample_minimize.low4000.out' % ( outdir, i_prev,j_prev) )


        ##########################################
        # CLUSTER! And keep a small number of representatives (400)
        ##########################################

        outfile_cluster = prefix+'sample.cluster.out'
        args_cluster = ' -cluster_test -in:file:silent %s  -in:file:silent_struct_type binary  -database %s  -radius %f -out:file:silent %s -nstruct %d ' % (string.join( combine_files ), DB,  CLUSTER_RADIUS, outfile_cluster, FINAL_NUMBER )
        condor_submit_cluster_file = 'CONDOR/region_%d_%d_cluster.condor' % (i,j)
        make_condor_submit_file( condor_submit_cluster_file, args_cluster, 1 )

        fid_dag.write('\nJOB %s %s\n' % (overall_job_tag,condor_submit_cluster_file) )
        fid_dag.write('PARENT %s CHILD %s\n' % (string.join(job_tags),overall_job_tag) )
        fid_dag.write('SCRIPT POST %s %s %s %s\n' % (overall_job_tag, POST_PROCESS_CLUSTER_SCRIPT, outfile_cluster, outdir ) )

        all_job_tags.append(  overall_job_tag )
