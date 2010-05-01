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
FINAL_NUMBER = parse_options( argv, "final_number", 40 )
SCORE_WEIGHTS = parse_options( argv, "weights", "score12.wts" )
PACK_WEIGHTS = parse_options( argv, "pack_weights", "pack.wts" )
CLUSTER_RADIUS = parse_options( argv, "cluster_radius", 0.5 )
AUTO_TUNE = parse_options( argv, "auto_tune", 0 )
template_pdb = parse_options( argv, "s", "1shf.pdb" )
native_pdb = parse_options( argv, "native", "1shf.pdb" )
cst_file = parse_options( argv, "cst_file", "" )
align_pdb = parse_options( argv, "align_pdb", "" )
pathway_file = parse_options( argv, "pathway_file", "" )
cluster_by_backbone_rmsd = parse_options( argv, "cluster_by_backbone_rmsd", 0 )
score_diff_cut = parse_options( argv, "score_diff_cut", 1000000.0 )
MAX_FRAG = parse_options( argv, "max_frag", len( sequence ) )


assert( exists( template_pdb ) )
assert( exists( native_pdb ) ) # Get rid of this later...

###############################################################
# Where's the executable?
###############################################################
HOMEDIR = expanduser('~')

EXE = HOMEDIR+'/src/mini/bin/stepwise_template_test.macosgccrelease'
if not( exists( EXE )):
    EXE = HOMEDIR+'/src/mini/bin/stepwise_template_test.linuxgccrelease'
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

fid_dag = open( "template_build.dag", 'w' )
fid_dag.write("DOT dag.dot\n")

###############################################################
# MAIN LOOP
###############################################################

# Loop over fragment lengths.
# Here make them in chunks of two to simplify this first calculation.

BLOCK_SIZE = 1

system( 'mkdir -p CONDOR' )
def make_condor_submit_file( condor_submit_file, arguments, queue_number, universe="vanilla" ):

    #print 'making condor_submit_file', condor_submit_file
    if not exists( dirname( condor_submit_file ) ):  system( 'mkdir -p '+ dirname( condor_submit_file ) )

    fid = open( condor_submit_file, 'w' )
    fid.write('+TGProject = TG-MCB090153\n')
    fid.write('universe = %s\n' % universe)
    fid.write('executable = %s\n' % EXE )
    fid.write('arguments = %s\n' % arguments)

    job_tag = basename(condor_submit_file).replace('.condor','')

    subdir = dirname( condor_submit_file ) + '/' + job_tag
    if not exists( subdir ):  system( 'mkdir -p '+subdir )

    fid.write('output = %s/$(Process).out\n' % subdir )
    fid.write('log = %s.log\n' % subdir )
    fid.write('error = %s/$(Process).err\n' % subdir)
    fid.write('notification = never\n')
    fid.write('Queue %d\n' % queue_number )
    fid.close()


def get_start_end( line ):
    in_seq = 0
    start_res = 1
    end_res = 0
    for k in range( len(line) ):
        if line[ k ] == ' ' or line[ k ] == '\n':
            if (in_seq):
                end_res = k
                break
        else:
            if not in_seq:
                start_res = k+1
                in_seq = 1

    if end_res == 0: end_res  = len( line )
    return ( start_res, end_res )

#####################################################
follow_path = 0
if len(pathway_file) > 0:
    follow_path = 1
    lines = open( pathway_file ).readlines()
    pathway_regions = []

    for line in lines:
        # don't yet know how to handle "merges" (i.e., inter-domain docking)
        if len( line ) > 5  and line[:5] == 'MERGE': break
        assert( line[:4] == 'PATH' )

        cols = map( lambda x:int(x), string.split( line )[1:] )
        i = cols[0]
        j = cols[0]
        cols = cols[1:]
        for m in cols:
            if ( m == i-1 ):
                i = m
            else:
                assert( m == j+1 )
                j = m
            pathway_regions.append( [i, j] )

    # old pathway file format -- way too verbose
    #for i in range( len( lines ) - 1 ) :
    #    line = lines[ i+1 ]
    #    ( start_res, end_res ) = get_start_end( line )
    #    region = [start_res, end_res]
    #    pathway_regions.append( region )

#####################################################
if AUTO_TUNE:
    cluster_tag = ' -auto_tune '
else:
    cluster_tag = ' -cluster:radius %s ' % CLUSTER_RADIUS

#####################################################
all_job_tags = []
real_compute_job_tags = []
jobs_done = []

# BASIC COMMAND
args0 = "-database %s  -s1 %s -native %s  -pack_weights %s -score:weights %s -ex1 -ex2 -extrachi_cutoff 0 -fasta %s   -constant_seed " % (DB, template_pdb, native_pdb, PACK_WEIGHTS, SCORE_WEIGHTS, fasta_file )

if len( cst_file ) > 0:
    assert( exists( cst_file ) )
    args0 += ' -cst_file %s ' % cst_file

if len( align_pdb ) > 0:
    assert( exists( align_pdb ) )
    args0 += ' -align_pdb %s ' % align_pdb

for L in range( 2, len(sequence)/BLOCK_SIZE + 1 ):

    chunk_length = BLOCK_SIZE * L;
    num_chunks = ( len( sequence) - chunk_length) / BLOCK_SIZE + 1

    for k in range( 1, num_chunks+1 ) :
        i = BLOCK_SIZE * ( k - 1 ) + 1
        j = i + chunk_length - 1

        if ( i < MIN_RES or j > MAX_RES ): continue

        if follow_path and ( [i,j] not in pathway_regions ): continue

        # Native PDB.
        prefix = 'region_%d_%d_' % (i,j)

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

        args = args0
        args += " -input_res1 "
        for k in range( i, j+1 ): args += "%d " % k
        args = args + " -slice_res1 "
        for k in range( i, j+1 ): args += "%d " % k


        ###########################################
        # DO THE JOBS
        start_regions = []

        ################################
        # WILL NEED TO UPDATE THIS!!
        ################################
        for k in range( 1, MAX_FRAG):
            i_prev = i
            j_prev = j - k
            prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
            if prev_job_tag in all_job_tags:   start_regions.append( [i_prev, j_prev ] )

        for k in range( 1, MAX_FRAG):
            i_prev = i + k
            j_prev = j
            prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
            if prev_job_tag in all_job_tags:   start_regions.append( [i_prev, j_prev ] )


        job_tags = []
        combine_files = []

        ##########################################
        # START FROM TEMPLATE
        ##########################################
        outfiles = []
        newdir = outdir+'/START_FROM_TEMPLATE'
        if not exists( newdir ): system( 'mkdir -p '+newdir )
        outfile = newdir + '/' + prefix + 'sample.out'
        start_tag = ' -start_from_template '

        job_tag = 'REGION_%d_%d_START_FROM_TEMPLATE' % (i,j)
        condor_submit_file = 'CONDOR/%s/%s.condor' %  ( overall_job_tag, job_tag )
        fid_dag.write('\nJOB %s %s\n' % (job_tag,condor_submit_file) )
        args2 = args + ' -out:file:silent %s ' % ( outfile )
        make_condor_submit_file( condor_submit_file, args2, 1 )
        fid_dag.write('SCRIPT POST %s %s %s\n' % (job_tag, POST_PROCESS_FILTER_SCRIPT,newdir) )

        job_tags.append( job_tag )
        real_compute_job_tags.append( job_tag )
        combine_files.append( '%s/start_from_template_sample.low4000.out' % outdir )


        ##########################################
        # APPEND OR PREPEND TO PREVIOUS PDB
        ##########################################
        for start_region in start_regions:
            i_prev = start_region[0]
            j_prev = start_region[1]

            dir_prev = 'REGION_%d_%d' % (i_prev, j_prev )
            prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)

            infile = 'region_%d_%d_sample.cluster.out' % (i_prev,j_prev)

            tag = 'S_$(Process)'

            #newdir = outdir+'/START_FROM_REGION_%d_%d' % (i_prev, j_prev )
            newdir = outdir+'/START_FROM_REGION_%d_%d_%s' % (i_prev, j_prev, tag )
            #if not exists( newdir ): system( 'mkdir -p '+newdir ) # In other dag-scripts, a pre-processing script "PRE" makes these directories.
            outfile = newdir + '/' + prefix + 'sample.out'

            args2 = "%s  -out:file:silent %s   -silent2 %s  -tags2 %s " % (args, outfile, infile, tag )
            args2 += " -input_res2 "
            for k in range( i_prev, j_prev+1 ): args2 += " %d " % k

            job_tag = 'REGION_%d_%d_START_FROM_REGION_%d_%d' % (i,j,i_prev,j_prev)
            condor_submit_file = 'CONDOR/%s/%s.condor' %  (overall_job_tag, job_tag)
            fid_dag.write('\nJOB %s %s\n' % (job_tag, condor_submit_file) )

            if not exists( condor_submit_file ):  make_condor_submit_file( condor_submit_file, args2, FINAL_NUMBER )

            if (prev_job_tag in all_job_tags)  and   (prev_job_tag not in jobs_done): #Note previous job may have been accomplished in a prior run -- not in the current DAG.
                fid_dag.write('PARENT %s  CHILD %s\n' % (prev_job_tag, job_tag) )

            # The pre process script finds out how many jobs there actually are...
            fid_dag.write('SCRIPT PRE %s   %s %s %s %s\n' % (job_tag, PRE_PROCESS_SETUP_SCRIPT,outdir,dir_prev,condor_submit_file) )
            fid_dag.write('SCRIPT POST %s %s %s/START_FROM_REGION_%d_%d\n' % (job_tag, POST_PROCESS_FILTER_SCRIPT,outdir,i_prev,j_prev ) )

            job_tags.append( job_tag )
            real_compute_job_tags.append( job_tag )
            combine_files.append( '%s/start_from_region_%d_%d_sample.low4000.out' % ( outdir, i_prev,j_prev) )



        if len( combine_files) > 0 : print 'DO_CHUNK',i,j

        ##########################################
        # CLUSTER! And keep a small number of representatives (400)
        ##########################################

        cluster_by_backbone_rmsd_tag = ''
        if cluster_by_backbone_rmsd: cluster_by_backbone_rmsd_tag = ' -cluster_by_backbone_rmsd '

        outfile_cluster = prefix+'sample.cluster.out'
        args_cluster = ' -cluster_test -in:file:silent %s  -in:file:silent_struct_type binary  -database %s  %s -out:file:silent %s -nstruct %d %s -score_diff_cut %8.3f' % (string.join( combine_files ), DB, cluster_tag, outfile_cluster, FINAL_NUMBER, cluster_by_backbone_rmsd_tag, score_diff_cut )
        condor_submit_cluster_file = 'CONDOR/%s/REGION_%d_%d_cluster.condor' % (overall_job_tag, i,j)

        make_condor_submit_file( condor_submit_cluster_file, args_cluster, 1, "scheduler" )

        fid_dag.write('\nJOB %s %s\n' % (overall_job_tag,condor_submit_cluster_file) )
        fid_dag.write('PARENT %s CHILD %s\n' % (string.join(job_tags),overall_job_tag) )
        fid_dag.write('SCRIPT POST %s %s %s %s\n' % (overall_job_tag, POST_PROCESS_CLUSTER_SCRIPT, outfile_cluster, outdir ) )

        all_job_tags.append(  overall_job_tag )


print "Number of rebuild jobs: ", len( real_compute_job_tags )
