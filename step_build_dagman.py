#!/usr/bin/python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
from parse_options import parse_options

fasta_file = parse_options( argv, "fasta", "1shf.fasta" )
assert( exists( fasta_file ) )
sequence = open( fasta_file  ).readlines()[1][:-1]

MIN_RES = parse_options( argv, "min_res", 1 )
MAX_RES = parse_options( argv, "max_res", len( sequence ) )
N_SAMPLE = parse_options( argv, "n_sample", 18 )
FINAL_NUMBER = parse_options( argv, "final_number", 40 )
SCORE_WEIGHTS = parse_options( argv, "weights", "score12.wts" )
PACK_WEIGHTS = parse_options( argv, "pack_weights", "pack.wts" )
NSTRUCT = parse_options( argv, "nstruct", 100 )
FILTER_RMSD = parse_options( argv, "filter_rmsd", 999.999 )
CLUSTER_RADIUS = parse_options( argv, "cluster_radius", 0.5 )
AUTO_TUNE = parse_options( argv, "auto_tune", 0 )
cst_file = parse_options( argv, "cst_file", "" )
cluster_by_all_atom_rmsd = parse_options( argv, "cluster_by_all_atom_rmsd", 0 )
score_diff_cut = parse_options( argv, "score_diff_cut", 1000000.0 )
add_peptide_plane = parse_options( argv, "add_peptide_plane", 0 )
make_native = parse_options( argv, "make_native", 0 )
align_pdb = parse_options( argv, "align_pdb", "" )
steps = parse_options( argv, "steps", [ "" ] )

native_pdb = parse_options( argv, "native", "1shf.pdb" )
assert( exists( native_pdb ) ) # Get rid of this later...

###############################################################
# Where's the executable?
###############################################################
HOMEDIR = expanduser('~rhiju')

EXE = HOMEDIR+'/src/mini/bin/stepwise_protein_test.macosgccrelease'
if not( exists( EXE )):
    EXE = HOMEDIR+'/src/mini/bin/stepwise_protein_test.linuxgccrelease'
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


assert( exists( SCORE_WEIGHTS ) or exists( DB + "/scoring/weights/"+SCORE_WEIGHTS) )
assert( exists( PACK_WEIGHTS ) or exists( DB + "/scoring/weights/"+PACK_WEIGHTS) )

fid_dag = open( "protein_build.dag", 'w' )
fid_dag.write("DOT dag.dot\n")

###############################################################
# MAIN LOOP
###############################################################

# Loop over fragment lengths.
# Here make them in chunks of two to simplify this first calculation.

BLOCK_SIZE = 1

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


if AUTO_TUNE:
    cluster_tag = ' -auto_tune '
else:
    cluster_tag = ' -cluster:radius %s ' % CLUSTER_RADIUS

all_job_tags = []
jobs_done = []

def what_is_the_step( step ):
    # later A can mean antiparallel; different tokens for parallel...
    if ( step.find('A') > 0 ):
        action = 'beta'
        cols = step.split('A')
        assert( len( cols ) == 2 )
        sample_residue = int( cols[0] )
        partner_residue = int( cols[1] )
    else:
        action = 'build'
        sample_residue = int( step )
        partner_residue = 0
    return ( sample_residue, action, partner_residue )

# Figure out jump_res from steps...
jump_res =[]
for s in steps:
    ( sample_res, action, partner_res ) = what_is_the_step( s )
    if action == "beta":
        jump_res.append( sample_res )
        jump_res.append( partner_res )

# BASIC COMMAND
extraflags = '-extrachi_cutoff 0 -ex1 -ex2 -score:weights %s -pack_weights %s' % (SCORE_WEIGHTS, PACK_WEIGHTS )
args = ' -out:file:silent_struct_type binary -database %s  -rebuild -native %s -fasta %s -n_sample %d -nstruct %d -minimize  -fullatom %s  -filter_rmsd %8.3f -cluster:radius 0.25 ' % ( DB, native_pdb, fasta_file, N_SAMPLE, NSTRUCT, extraflags, FILTER_RMSD )

if add_peptide_plane: args += ' -add_peptide_plane '
if len( cst_file ) > 0:
    assert( exists( cst_file ) )
    args += ' -cst_file %s ' % cst_file
if len( align_pdb ) > 0:
    assert( exists( align_pdb ) )
    args += ' -align_pdb %s ' % align_pdb
if len( jump_res ) > 0:
    args += ' -jump_res'
    for m in jump_res: args += ' %d' % m

step_native_pdbs = []

res_build = []
for s in range( 0, len(steps) ):

    prefix = 'step_%d_' % (s)
    overall_job_tag = 'STEP_%d' % s

    ###########################################
    # OUTPUT DIRECTORY
    (sample_res, action, partner_res ) = what_is_the_step( steps[ s ] )

    res_prev_build = []
    for m in res_build: res_prev_build.append( m )

    assert( sample_res not in res_prev_build )
    res_build.append( sample_res )
    res_build.sort()

    partner_res_tag = ''
    if partner_res > 0: partner_res_tag = '%d' % partner_res
    print '%10s ==> %3d  [%s %s]' % ( overall_job_tag, sample_res, action, partner_res_tag ),

    if make_native and len( native_pdb ) > 0:
        step_native_pdb = 'step_%d_%s' % (s,native_pdb)
        step_native_pdbs.append( step_native_pdb )
        if not exists( 'NATIVES/'+step_native_pdb ):
            if not exists( 'NATIVES' ): system( 'mkdir -p NATIVES')
            command = PYDIR+'/pdbslice.py '+native_pdb
            command += ' -subset'
            for m in res_build: command += ' %d' % m
            command += ' step_%d_' % s
            print command
            system( command )
            system( 'mv %s NATIVES/%s' % (step_native_pdb, step_native_pdb) )

    # This job is maybe already done...
    outfile_cluster = prefix+'sample.cluster.out'
    if exists( outfile_cluster ):
        all_job_tags.append(  overall_job_tag )
        jobs_done.append( overall_job_tag   )
        print 'DONE'
        continue
    else:
        print

    outdir = overall_job_tag
    if not( exists( outdir) ):
        system( 'mkdir -p ' + outdir )

    #overall_job_tag = 'REGION_%d_%d' % (i,j)

    ###########################################
    # DO THE JOBS
    job_tags = []
    combine_files = []

    if len( res_build ) <= 1:
        ##########################################
        # START FROM SCRATCH
        ##########################################
        outfiles = []
        newdir = outdir+'/START_FROM_SCRATCH'
        if not exists( newdir ): system( 'mkdir -p '+newdir )
        outfile = newdir + '/' + prefix + 'sample.out'

        job_tag = '%s_START_FROM_SCRATCH' % overall_job_tag
        condor_submit_file = 'CONDOR/%s.condor' %  job_tag
        fid_dag.write('\nJOB %s %s\n' % (job_tag,condor_submit_file) )

        args2 = '%s -out:file:silent %s ' % (args, outfile )
        args2 += ' -sample_res '
        args2 += ' %d' % sample_res

        args2 += ' -superimpose_res %d ' % sample_res

        make_condor_submit_file( condor_submit_file, args2, 1 )
        fid_dag.write('SCRIPT POST %s %s %s\n' % (job_tag, POST_PROCESS_FILTER_SCRIPT,newdir) )

        job_tags.append( job_tag )
        combine_files.append( '%s/start_from_scratch_sample.low4000.out' % outdir )

    else:

        s_prev = s - 1

        prev_job_tag = 'STEP_%d' % (s_prev)
        dir_prev = prev_job_tag

        infile = 'step_%d_sample.cluster.out' % (s_prev)

        tag = 'S_$(Process)'

        newdir = outdir+'/START_FROM_STEP_%d_%s' % (s_prev, tag.upper() )
        outfile = newdir + '/' + prefix + 'sample.out'

        args2 = '%s -out:file:silent %s -silent1 %s -tags1 %s' % (args, outfile, infile, tag )

        args2 += ' -input_res1 '
        for m in res_prev_build: args2 += ' %d' % m

        args2 += ' -superimpose_res %d ' % sample_res
        for m in res_prev_build: args2 += ' %d' % m

        args2 += ' -sample_res '

        all_sample_res = [ sample_res ]
        pos = res_build.index( sample_res )

        if ( pos > 0 and res_build[ pos - 1] == sample_res-1 ):
            all_sample_res.append( sample_res - 1)
        elif ( pos < len(res_build)-1 and res_build[ pos + 1 ] == sample_res+1 ):
            all_sample_res.append( sample_res + 1)
        else:
            assert( partner_res in res_prev_build )

        all_sample_res.sort()
        for m in all_sample_res: args2 += ' %d' % m

        if action == 'beta':
            assert( len(all_sample_res) == 1 )
            args2 += ' -sample_beta'

        job_tag = '%s_START_FROM_%s' % (overall_job_tag, prev_job_tag )
        condor_submit_file = 'CONDOR/%s.condor' %  job_tag
        fid_dag.write('\nJOB %s %s\n' % (job_tag, condor_submit_file) )

        if not exists( condor_submit_file ):  make_condor_submit_file( condor_submit_file, args2, FINAL_NUMBER )

        if (prev_job_tag in all_job_tags)  and   (prev_job_tag not in jobs_done): #Note previous job may have been accomplished in a prior run -- not in the current DAG.
            fid_dag.write('PARENT %s  CHILD %s\n' % (prev_job_tag, job_tag) )

        # The pre process script finds out how many jobs there actually are...
        fid_dag.write('SCRIPT PRE %s   %s %s %s %s\n' % (job_tag, PRE_PROCESS_SETUP_SCRIPT,outdir,dir_prev,condor_submit_file) )
        fid_dag.write('SCRIPT POST %s %s %s/START_FROM_%s\n' % (job_tag, POST_PROCESS_FILTER_SCRIPT,outdir,prev_job_tag ) )

        job_tags.append( job_tag )
        combine_files.append( '%s/start_from_%s_sample.low4000.out' % ( outdir,prev_job_tag.lower() ) )


    ##########################################
    # CLUSTER! And keep a small number of representatives (400)
    ##########################################

    cluster_by_all_atom_rmsd_tag = ''
    if cluster_by_all_atom_rmsd: cluster_by_all_atom_rmsd_tag = ' -cluster_by_all_atom_rmsd '

    outfile_cluster = prefix+'sample.cluster.out'
    args_cluster = ' -cluster_test -in:file:silent %s  -in:file:silent_struct_type binary  -database %s  %s -out:file:silent %s -nstruct %d %s -score_diff_cut %8.3f' % (string.join( combine_files ), DB,  cluster_tag, outfile_cluster, FINAL_NUMBER, cluster_by_all_atom_rmsd_tag, score_diff_cut )
    condor_submit_cluster_file = 'CONDOR/STEP_%d_cluster.condor' % (s)
    make_condor_submit_file( condor_submit_cluster_file, args_cluster, 1, "scheduler" )

    fid_dag.write('\nJOB %s %s\n' % (overall_job_tag,condor_submit_cluster_file) )
    fid_dag.write('PARENT %s CHILD %s\n' % (string.join(job_tags),overall_job_tag) )
    fid_dag.write('SCRIPT POST %s %s %s %s\n' % (overall_job_tag, POST_PROCESS_CLUSTER_SCRIPT, outfile_cluster, outdir ) )

    all_job_tags.append(  overall_job_tag )


command = PYDIR+'/catpdb.py '
for file in step_native_pdbs: command += ' NATIVES/%s' % file
command += ' > all_natives.pdb'
print( command )
system( command )
