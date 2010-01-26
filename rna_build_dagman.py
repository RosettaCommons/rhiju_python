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

# Define sequence
fasta_file = parse_options( argv, "fasta", "1shf.fasta" )
assert( exists( fasta_file ) )
sequence = open( fasta_file  ).readlines()[1][:-1]

# Defines loop residues.
MIN_RES = parse_options( argv, "min_res", 1 )
MAX_RES = parse_options( argv, "max_res", len( sequence ) )

# This is a dummy number -- typical number of pdbs in an ensemble.
# For each job, it will be replaced with a more accurate number, based
#  on the results of parent jobs.
FINAL_NUMBER = 108

# Starting PDB files.
input_pdb = parse_options( argv, "s", "1shf.pdb" )
native_pdb = parse_options( argv, "native", "1shf.pdb" )
score_diff_cut = parse_options( argv, "score_diff_cut", 1000000.0 )
CLUSTER_RADIUS = parse_options( argv, "cluster_radius", 0.25 )
cutpoint_open = parse_options( argv, "cutpoint_open", 0 )

cutpoints_open = []
if (cutpoint_open > 0 ): cutpoints_open.append( cutpoint_open )

assert( exists( native_pdb ) ) # Get rid of this later...
assert( exists( input_pdb ) ) # Get rid of this later...


###############################################################
# Where's the executable?
###############################################################
HOMEDIR = expanduser('~')

EXE = HOMEDIR+'/src/mini_swa_rna/bin/rna_swa_test.macosgccrelease'
if not( exists( EXE )):
    EXE = HOMEDIR+'/src/mini_swa_rna/bin/rna_swa_test.linuxgccrelease'
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

fid_dag = open( "rna_build.dag", 'w' )
fid_dag.write("DOT dag.dot\n")

#############################################
# Should probably put in a check to make sure
# that starting pdb ("input_pdb") has the correct
# sequence!
starting_sequence = sequence[:(MIN_RES-1)] + sequence[MAX_RES:]
input_pdb_sequence = popen( '%s/pdb2fasta.py %s ' % (PYDIR,input_pdb) ).readlines()[1][:-1]
assert( starting_sequence == input_pdb_sequence )

###############################################################
# MAIN LOOP
###############################################################

system( 'mkdir -p CONDOR' )
def make_condor_submit_file( condor_submit_file, arguments, queue_number, universe="vanilla" ):

    fid = open( condor_submit_file, 'w' )
    fid.write('+TGProject = TG-MCB090153\n')
    fid.write('universe = %s\n' % universe)
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

# Keep a list of all regions that we are building.
# Initialize with starting pdb -- already built!
all_job_tags = []
jobs_done = []

input_file_tag = 'REGION_%d_%d' % (MIN_RES-1, MAX_RES+1)
all_job_tags.append( input_file_tag )
jobs_done.append( input_file_tag )

loop_length = MAX_RES - MIN_RES + 1

# Order calculation based on number of residues built so far.
for L in range( 1, loop_length+1 ):

    # k is the number of residues built from the 5' side. Goes from 0 to L.
    for k in range( L+1 ) :

        i = MIN_RES - 1 + k
        j = MAX_RES + 1 - (L - k)

        # Later, we can copy/paste in a specific path.
        # Not implemented currently
        # if follow_path and ( [i,j] not in pathway_regions ): continue

        # Native PDB.
        prefix = 'region_%d_%d_' % (i,j)
        print 'Loop boundaries --> %d-%d and %d-%d' % (1,i,j,len(sequence))

        # This job is maybe already done...
        outfile_cluster = prefix+'sample.cluster.out'
        overall_job_tag = 'REGION_%d_%d' % (i,j)



        if exists( outfile_cluster ):
            all_job_tags.append(  overall_job_tag )
            jobs_done.append( overall_job_tag   )
            continue

        ###########################################
        # OUTPUT DIRECTORY
        outdir = 'REGION_%d_%d' % (i,j)
        if not( exists( outdir) ):
            system( 'mkdir -p ' + outdir )

        # BASIC COMMAND
        args = ' -algorithm rna_resample_test  -database %s -fasta %s -native %s -output_virtual ' %  ( DB, fasta_file, native_pdb )

        if len( cutpoints_open ) > 0:
            args += ' -cutpoint_open '
            for cutpos in cutpoints_open: args += '%d ' % cutpos

        ###########################################
        # DO THE JOBS
        start_regions = []

        i_prev = i - 1
        j_prev = j
        prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
        if (prev_job_tag in all_job_tags) and (i_prev not in cutpoints_open):   start_regions.append( [i_prev, j_prev ] )

        i_prev = i
        j_prev = j+1
        prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
        if (prev_job_tag in all_job_tags) and (j not in cutpoints_open):   start_regions.append( [i_prev, j_prev ] )

        job_tags = []
        combine_files = []

        for start_region in start_regions:
            i_prev = start_region[0]
            j_prev = start_region[1]

            dir_prev = 'REGION_%d_%d' % (i_prev, j_prev )

            # prev_job_tag and dir_prev are currently the same, but may change directory structure in the future.
            prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)

            num_jobs_to_queue = FINAL_NUMBER

            if prev_job_tag == input_file_tag:
                args2 = '%s -s %s ' % (args, input_pdb)
                newdir = outdir+'/START_FROM_REGION_%d_%d' % (i_prev, j_prev )
                if not exists( newdir ): system( 'mkdir -p '+newdir )
                outfile = newdir + '/' + prefix + 'sample.out'
                num_jobs_to_queue = 1
            else:
                infile = 'region_%d_%d_sample.cluster.out' % (i_prev,j_prev)
                tag = 'S_$(Process)'
                args2 = '%s -in:file:silent_struct_type binary_rna -in:file:silent %s -tags %s ' % (args, infile, tag )
                newdir = outdir+'/START_FROM_REGION_%d_%d_%s' % (i_prev, j_prev, tag.upper() )
                outfile = newdir + '/' + prefix + 'sample.out'

            args2 += '-out:file:silent %s ' % outfile


            # What is already built? What will move?
            args2 += ' -input_res '
            for m in range(1,i_prev+1): args2 += ' %d' % m
            for m in range(j_prev, len(sequence)+1): args2 += ' %d' % m

            if ( i == i_prev ):
                moving_res = j
            else:
                moving_res = i
            args2 += ' -sample_res %d ' % moving_res

            # Special --> close cutpoint.
            if ( L == loop_length ):
                args2 += ' -cutpoint_closed %d ' % i

            job_tag = 'REGION_%d_%d_START_FROM_REGION_%d_%d' % (i,j,i_prev,j_prev)
            condor_submit_file = 'CONDOR/%s.condor' %  job_tag
            fid_dag.write('\nJOB %s %s\n' % (job_tag, condor_submit_file) )

            if not exists( condor_submit_file ):
                make_condor_submit_file( condor_submit_file, args2, num_jobs_to_queue )

            if (prev_job_tag in all_job_tags)  and   (prev_job_tag not in jobs_done): #Note previous job may have been accomplished in a prior run -- not in the current DAG.
                fid_dag.write('PARENT %s  CHILD %s\n' % (prev_job_tag, job_tag) )

            # The pre process script finds out how many jobs there actually are...
            if ( prev_job_tag != input_file_tag ):
                fid_dag.write('SCRIPT PRE %s   %s %s %s %s\n' % (job_tag, PRE_PROCESS_SETUP_SCRIPT,outdir,dir_prev,condor_submit_file) )
            fid_dag.write('SCRIPT POST %s %s %s/START_FROM_REGION_%d_%d\n' % (job_tag, POST_PROCESS_FILTER_SCRIPT,outdir,i_prev,j_prev ) )

            job_tags.append( job_tag )
            combine_files.append( '%s/start_from_region_%d_%d_sample.low4000.out' % ( outdir, i_prev,j_prev) )


        ##########################################
        # CLUSTER! And keep a small number of representatives (400)
        ##########################################

        if len( combine_files ) == 0: continue

        outfile_cluster = prefix+'sample.cluster.out'
        args_cluster = ' -algorithm cluster_old -in:file:silent %s  -in:file:silent_struct_type binary_rna  -database %s  -radius %f -out:file:silent %s  -score_diff_cut %8.3f' % (string.join( combine_files ), DB,  CLUSTER_RADIUS, outfile_cluster, score_diff_cut )

        condor_submit_cluster_file = 'CONDOR/REGION_%d_%d_cluster.condor' % (i,j)
        make_condor_submit_file( condor_submit_cluster_file, args_cluster, 1, "scheduler" )

        fid_dag.write('\nJOB %s %s\n' % (overall_job_tag,condor_submit_cluster_file) )
        fid_dag.write('PARENT %s CHILD %s\n' % (string.join(job_tags),overall_job_tag) )
        fid_dag.write('SCRIPT POST %s %s %s %s\n' % (overall_job_tag, POST_PROCESS_CLUSTER_SCRIPT, outfile_cluster, outdir ) )

        all_job_tags.append(  overall_job_tag )
