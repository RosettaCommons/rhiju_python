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
ZIGZAG = parse_options( argv, "zigzag", 0 )
N_SAMPLE = parse_options( argv, "n_sample", 18 )
FINAL_NUMBER = parse_options( argv, "final_number", 100 )
SCORE_WEIGHTS = parse_options( argv, "weights", "score12_no_hb_env_dep.wts" )
PACK_WEIGHTS = parse_options( argv, "pack_weights", "pack_no_hb_env_dep.wts" )
NSTRUCT = parse_options( argv, "nstruct", 100 )
FILTER_RMSD = parse_options( argv, "filter_rmsd", -1.0 )
CLUSTER_RADIUS = parse_options( argv, "cluster_radius", 0.25 )
CLUSTER_RADIUS_SAMPLE = parse_options( argv, "cluster_radius_sample", 0.1 )
AUTO_TUNE = parse_options( argv, "auto_tune", 0 )
filter_native_big_bins = parse_options( argv, "filter_native_big_bins", 0 )
score_diff_cut = parse_options( argv, "score_diff_cut", 10.0 )
max_res_to_add_denovo = parse_options( argv, "denovo", 0 )

native_pdb = parse_options( argv, "native", "" )
template_pdb = parse_options( argv, "s", "" )
cst_file = parse_options( argv, "cst_file", "" )
pathway_file = parse_options( argv, "pathway_file", "" )
cluster_by_all_atom_rmsd = parse_options( argv, "cluster_by_all_atom_rmsd", 0 )
add_peptide_plane = parse_options( argv, "add_peptide_plane", 0 )
BUILD_BOTH_TERMINI = parse_options( argv, "build_both_termini", 0 )
MAX_FRAG = parse_options( argv, "max_frag", 16 )
max_length = parse_options( argv, "max_length", 0 )
superimpose_res = parse_options( argv, "superimpose_res", [ -1 ] )
align_pdb = parse_options( argv, "align_pdb", "" )


DENOVO = ( max_res_to_add_denovo > 0 )
TEMPLATE = len( template_pdb ) > 0

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

if not exists( 'CONDOR/' ):
    system( 'mkdir -p CONDOR' )


#########################################################
# Some useful functions (move somewhere else?)
#########################################################
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

def setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tag, args2, decoy_tag,\
                                         fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files):
    newdir = overall_job_tag+'/'+sub_job_tag
    if len( decoy_tag ) > 0:
        newdir += '_' + decoy_tag
    outfile = newdir + '/' + overall_job_tag.lower() + '_sample.out'

    args2 += ' -out:file:silent %s ' % outfile

    job_tag = '%s_%s' % ( overall_job_tag,  sub_job_tag )
    condor_submit_file = 'CONDOR/%s.condor' %  job_tag

    fid_dag.write('\nJOB %s %s\n' % (job_tag, condor_submit_file) )

    if not exists( condor_submit_file ):
        make_condor_submit_file( condor_submit_file, args2, FINAL_NUMBER )

    if (prev_job_tag in all_job_tags)  and  (prev_job_tag not in jobs_done): #Note previous job may have been accomplished in a prior run -- not in the current DAG.
        fid_dag.write('PARENT %s  CHILD %s\n' % (prev_job_tag, job_tag) )

    # The pre process script finds out how many jobs there actually are...
    if len( decoy_tag ) > 0:
        fid_dag.write('SCRIPT PRE %s   %s %s %s %s %s\n' % (job_tag, PRE_PROCESS_SETUP_SCRIPT,overall_job_tag,prev_job_tag,condor_submit_file,sub_job_tag) )
    else:
        if not exists( newdir ): system( 'mkdir -p '+newdir )

    fid_dag.write('SCRIPT POST %s %s %s/%s\n' % (job_tag, POST_PROCESS_FILTER_SCRIPT,overall_job_tag,sub_job_tag ) )

    # In python, lists are passed by reference... these should get updated for the outside world.
    job_tags.append( job_tag )
    real_compute_job_tags.append( job_tag )
    combine_files.append( '%s/%s_sample.low4000.out' % ( overall_job_tag, sub_job_tag.lower() ) )


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

##########################
# BASIC COMMAND
##########################

args = ' -database %s  -rebuild -out:file:silent_struct_type binary  -fasta %s -n_sample %d -nstruct %d -cluster:radius %8.3f' % ( DB, fasta_file, N_SAMPLE, NSTRUCT, CLUSTER_RADIUS_SAMPLE )

args += ' -extrachi_cutoff 0 -ex1 -ex2' # These may be redundant actually.

args += ' -score:weights %s -pack_weights %s' % (SCORE_WEIGHTS, PACK_WEIGHTS )

if ( FILTER_RMSD > 0.0 ):
    args += ' -filter_rmsd %8.3f' % FILTER_RMSD

if add_peptide_plane: args += ' -add_peptide_plane'
if filter_native_big_bins:  args+= ' -filter_native_big_bins' # this is defunct now, I think
if len( cst_file ) > 0:
    assert( exists( cst_file ) )
    args += ' -cst_file %s' % cst_file
if len( align_pdb ) > 0:
    assert( exists( align_pdb ) )
    args += ' -align_pdb %s' % align_pdb
if len( native_pdb ) > 0:
    assert( exists( native_pdb ) )
    args += ' -native %s' % native_pdb
if len( superimpose_res ) > 0:
    args += ' -superimpose_res '
    for k in superimpose_res: args += '%d ' % k

args_START = args

if AUTO_TUNE:
    cluster_tag = ' -auto_tune '
else:
    cluster_tag = ' -cluster:radius %s ' % CLUSTER_RADIUS


################################
# MAIN LOOP
################################
# Loop over fragment lengths.

all_job_tags = []
real_compute_job_tags = []
jobs_done = []

for L in range( 2, len(sequence) + 1 ):

    chunk_length = L;
    num_chunks = ( len( sequence) - chunk_length) + 1

    if ( max_length > 0 and L > max_length ): continue

    for k in range( 1, num_chunks+1 ) :
        i = k
        j = i + chunk_length - 1

        if ( i < MIN_RES or j > MAX_RES ): continue

        #ZIGZAG!! special case for beta hairpins.
        if ( ZIGZAG and abs( ( i - MIN_RES ) - ( MAX_RES - j ) ) > 1 ) : continue

        if follow_path and ( [i,j] not in pathway_regions ): continue

        overall_job_tag = 'REGION_%d_%d' % (i,j)

        print 'Do region ==> ',i,j,

        # This job is maybe already done...
        outfile_cluster = overall_job_tag.lower()+'_sample.cluster.out'
        if exists( outfile_cluster ):
            all_job_tags.append(  overall_job_tag )
            jobs_done.append( overall_job_tag   )
            print 'DONE'
            continue
        else:
            print

        # OUTPUT DIRECTORY
        if not( exists( overall_job_tag) ):
            system( 'mkdir -p ' + overall_job_tag )

        termini_tag = ""
        if ( i == 1 ): termini_tag += " -n_terminus"
        if ( j == len(sequence)  ): termini_tag += " -c_terminus"
        args = args_START + termini_tag


        ###########################################
        # DO THE JOBS
        ###########################################
        start_regions = []

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

        if BUILD_BOTH_TERMINI:
            i_prev = i + 1
            j_prev = j - 1
            prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
            if prev_job_tag in all_job_tags:   start_regions.append( [i_prev, j_prev ] )

        job_tags = []
        combine_files = []

        ########################################################
        # Several "modes" -- different stuff for the grinder!
        #
        #  currently have -- DENOVO, TEMPLATE
        ########################################################

        if TEMPLATE:
            # Good to just build the whole thing off template.
            sub_job_tag = "START_FROM_TEMPLATE"

            args2 = args
            args2 += ' -prepack -s1 ' + template_pdb
            args2 += ' -input_res1 '
            for k in range(i,j+1): args2 += ' %d' % k
            args2 += ' -slice_res1 '
            for k in range(i,j+1): args2 += ' %d' % k

            setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tag, args2, '', \
                                                 fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)

        if len( start_regions ) == 0 and DENOVO: # This happens for two-residue fragments.

            sub_job_tag = "START_FROM_SCRATCH"

            args2 = args
            args2 += ' -sample_res %d %d ' % (i,j)

            setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tag, args2, '', \
                                                 fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)


        ##########################################
        # APPEND OR PREPEND TO PREVIOUS PDB
        ##########################################
        for start_region in start_regions:

            i_prev = start_region[0]
            j_prev = start_region[1]

            prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)

            infile = 'region_%d_%d_sample.cluster.out' % (i_prev,j_prev)

            if ( DENOVO and \
                 abs(i - i_prev ) <= max_res_to_add_denovo and \
                 abs(j - j_prev ) <= max_res_to_add_denovo) :

                sub_job_tag = 'START_FROM_REGION_%d_%d_DENOVO' % ( i_prev, j_prev )

                decoy_tag = 'S_$(Process)'
                args2 = '%s  -silent1 %s -tags1 %s' % (args, infile, decoy_tag )

                args2 += ' -input_res1 '
                for m in range(i_prev,j_prev+1): args2 += ' %d' % m

                args2 += ' -sample_res '
                if ( i < i_prev and j == j_prev):
                    for m in range(i,i_prev+1): args2 += ' %d' % m
                elif ( i == i_prev and j > j_prev ):
                    for m in range(j_prev,j+1): args2 += ' %d' % m
                else:
                    for m in [i,j]: args2 += ' %d' % m

                setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tag, args2, decoy_tag, \
                                                     fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)



            if TEMPLATE:
                sub_job_tag = 'START_FROM_REGION_%d_%d_TEMPLATE' % ( i_prev, j_prev )

                args2 = "%s  -silent1 %s " % (args, infile )
                args2 += " -input_res1 "
                input_res1 =  range( i_prev, j_prev+1 )
                for m in input_res1: args2 += ' %d' % m

                args2 += " -s2 %s" % template_pdb
                input_string = ''
                for m in range( i, j+1 ):
                    if m not in input_res1:
                        input_string += ' %d' % m
                args2 += " -input_res2 %s -slice_res2 %s" % (input_string, input_string)

                args2 += ' -sample_res '
                for m in range( i, j+1 ): args2 += ' %d' % m

                setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tag, args2, '', \
                                                     fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)



        if len( combine_files ) == 0:
            print "No jobs for ", overall_job_tag
            print "Something does not look right! Exiting."
            exit( 0 )

        ################################################################
        # CLUSTER! And keep a small number of representatives (400)
        ################################################################

        cluster_by_all_atom_rmsd_tag = ''
        if cluster_by_all_atom_rmsd: cluster_by_all_atom_rmsd_tag = ' -cluster_by_all_atom_rmsd '

        outfile_cluster = overall_job_tag.lower()+'_sample.cluster.out'
        args_cluster = ' -cluster_test -in:file:silent %s  -in:file:silent_struct_type binary  -database %s  %s -out:file:silent %s -nstruct %d %s -score_diff_cut %8.3f' % (string.join( combine_files ), DB,  cluster_tag, outfile_cluster, FINAL_NUMBER, cluster_by_all_atom_rmsd_tag, score_diff_cut )
        condor_submit_cluster_file = 'CONDOR/REGION_%d_%d_cluster.condor' % (i,j)
        make_condor_submit_file( condor_submit_cluster_file, args_cluster, 1, "scheduler" )

        fid_dag.write('\nJOB %s %s\n' % (overall_job_tag,condor_submit_cluster_file) )
        fid_dag.write('PARENT %s CHILD %s\n' % (string.join(job_tags),overall_job_tag) )
        fid_dag.write('SCRIPT POST %s %s %s %s\n' % (overall_job_tag, POST_PROCESS_CLUSTER_SCRIPT, outfile_cluster, overall_job_tag ) )

        all_job_tags.append(  overall_job_tag )

print
print "Total number of jobs to run (not counting clustering):", len( real_compute_job_tags )
print "Total number of final outfiles (and clustering jobs):", len( all_job_tags )
