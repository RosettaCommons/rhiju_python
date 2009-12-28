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
FINAL_NUMBER = parse_options( argv, "final_number", 400 )
NSTRUCT = parse_options( argv, "nstruct", FINAL_NUMBER )
CLUSTER_RADIUS = parse_options( argv, "cluster_radius", 2.0 )
filter_native_big_bins = parse_options( argv, "filter_native_big_bins", 0 )
native_pdb = parse_options( argv, "native", "1shf.pdb" )
cst_file = parse_options( argv, "cst_file", "" )
pathway_file = parse_options( argv, "pathway_file", "" )
internal_score_diff_cut = parse_options( argv, "internal_score_diff_cut", 5.0 )
score_diff_cut = parse_options( argv, "score_diff_cut", 10.0 )
max_overlap = parse_options( argv, "max_overlap", 10 )
centroid_weights = parse_options( argv, "centroid_weights", "score3.wts" )
skip_internal = parse_options( argv, "skip_internal", 0 )

if not centroid_weights == "score3.wts": assert( exists( centroid_weights) )

frag_files = {}
frag3 = parse_options( argv, "frag3", "" )
frag9 = parse_options( argv, "frag9", "" )

if len( frag3 ) > 0:
    frag_files[3] = frag3
    assert( exists( frag3 )  )
if len( frag9 ) > 0:
    frag_files[9] = frag9
    assert( exists( frag9 )  )

frag_lengths = frag_files.keys()


assert( exists( native_pdb ) ) # Can get rid of this later.

###############################################################
# Where's the executable?
###############################################################
HOMEDIR = expanduser('~')

EXE = HOMEDIR+'/src/mini/bin/stepwise_centroid_test.macosgccrelease'
if not( exists( EXE )):
    EXE = HOMEDIR+'/src/mini/bin/stepwise_centroid_test.linuxgccrelease'
assert( exists( EXE ) )

DB = HOMEDIR+'/minirosetta_database'
assert( exists( DB ) )

#PYDIR = HOMEDIR+'/python'
#assert( exists( PYDIR ) )

fid_dag = open( "centroid_build.dag", 'w' )
fid_dag.write("DOT dag.dot\n")

###############################################################
# MAIN LOOP
###############################################################

# Loop over fragment lengths.
# Here make them in chunks of two to simplify this first calculation.

BLOCK_SIZE = 1

system( 'mkdir -p CONDOR' )
def make_condor_submit_file( condor_submit_file, arguments, queue_number ):

    subdir =  dirname( condor_submit_file )+'/'
    if not exists( subdir ): system( 'mkdir -p '+subdir )

    fid = open( condor_submit_file, 'w' )
    fid.write('+TGProject = TG-MCB090153\n')
    fid.write('universe = vanilla\n')
    fid.write('executable = %s\n' % EXE )
    fid.write('arguments = %s\n' % arguments)

    job_tag = basename(condor_submit_file).replace('.condor','')


    fid.write('output = %s/$(Process).out\n' % subdir )
    fid.write('log = %s/%s.log\n' % (subdir,job_tag) )
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


follow_path = 0
if len(pathway_file) > 0:
    follow_path = 1
    lines = open( pathway_file ).readlines()
    pathway_regions = []
    parent_region = {}
    for i in range( len( lines ) - 1 ) :
        line = lines[ i+1 ]
        ( start_res, end_res ) = get_start_end( line )
        region = [start_res, end_res]
        pathway_regions.append( region )

        line_prev = lines[ i ]
        ( start_res_prev, end_res_prev ) = get_start_end( line_prev )
        if start_res_prev < end_res_prev:
            region_prev = [start_res_prev, end_res_prev]
            region_tag = 'REGION_%d_%d' % (region[0],region[1])
            parent_region[ region_tag ] = region_prev

all_job_tags = []
jobs_done = []

# Assume we have 3 mer and 9 mer frag files. So minimum pose size is 3.
frag_lengths_sort = frag_lengths
frag_lengths_sort.sort()
MIN_LENGTH = ( frag_lengths_sort )[0]
#print "MIN_LENGTH", MIN_LENGTH

for L in range( MIN_LENGTH, len(sequence) + 1 ):

    chunk_length = L;
    num_chunks = ( len( sequence) - chunk_length) + 1

    for k in range( 1, num_chunks+1 ) :
        i = k;
        j = i + chunk_length - 1

        the_very_last_region = 0
        if ( i == 1 and j == len( sequence ) ): the_very_last_region = 1

        if ( i < MIN_RES or j > MAX_RES ): continue

        if follow_path and ( [i,j] not in pathway_regions ): continue

        # Native PDB.
        prefix = 'region_%d_%d_' % (i,j)
        print 'DO_CHUNK',i,j

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
        args = ' -database %s -native %s -fasta %s -nstruct %d  -radius 1.0 -min_res %d -max_res %d -score_diff_cut %8.3f -centroid_weights %s' %\
               ( DB, native_pdb, fasta_file,  NSTRUCT, i, j, internal_score_diff_cut, centroid_weights )

        if len( cst_file ) > 0:
            assert( exists( cst_file ) )
            args += ' -cst_file %s ' % cst_file


        #overall_job_tag = 'REGION_%d_%d' % (i,j)

        ###########################################
        # DO THE JOBS
        start_regions = []

        # A job is defined by:
        #      - starting region
        #      - which frag file to use.
        #      [ later --> starting regions and insertion point for frag file ]

        job_parameters = []
        if follow_path:
            region_tag = overall_job_tag
            if region_tag in parent_region.keys():
                region_prev = parent_region[ region_tag ]
                start_regions.append( [ region_prev[0], region_prev[1] ] )
        else:

            # Later do a loop over 3mers, 9mers... 1mers?
            for frag_length in frag_lengths:

                if ( i + frag_length - 1 == j ):
                    # start from scratch
                    job_parameters.append( [ [] , 1, frag_length] )
                else:
                    # prepend.
                    j_prev = j
                    insert_res = i
                    last_frag_res = insert_res + frag_length - 1
                    for i_prev in range( i + 1, i + frag_length ):
                        prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
                        if ( last_frag_res - i_prev ) > max_overlap and ( not the_very_last_region ): continue
                        if ( prev_job_tag in all_job_tags )  and  ( insert_res >= i ) and ( last_frag_res <= j ) :
                            job_parameters.append( [ [[i_prev, j_prev]], insert_res, frag_length ] )

                    # append
                    i_prev = i
                    insert_res = j - frag_length + 1
                    last_frag_res = insert_res + frag_length - 1
                    for j_prev in range( j - frag_length + 1, j):
                        prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
                        if ( j_prev - insert_res ) > max_overlap and ( not the_very_last_region ): continue
                        if ( prev_job_tag in all_job_tags ) and ( insert_res >= i ) and ( last_frag_res <= j ):
                            job_parameters.append( [ [[i_prev, j_prev]], insert_res, frag_length ] )

                    if (skip_internal): continue

                    # internal fragment insertion
                    # assume previous start pdb's overlap with first 4 and last 4 residues of 9mer fragment...
                    # this could be generalized with different overlaps, different fragment lengths, though
                    # at the expense of more computation
                    if ( frag_length == 9):
                        for insert_res in range( i + 1, j - frag_length + 1):
                            i_prev1 = i;
                            j_prev1 = insert_res + 3
                            prev_job_tag1 = 'REGION_%d_%d' % (i_prev1,j_prev1)

                            i_prev2 = insert_res + 5
                            j_prev2 = j
                            prev_job_tag2 = 'REGION_%d_%d' % (i_prev2,j_prev2)

                            if prev_job_tag1 not in all_job_tags: continue
                            if prev_job_tag2 not in all_job_tags: continue

                            job_parameters.append( [ [[i_prev1, j_prev1],[i_prev2,j_prev2]], insert_res, frag_length ] )


        job_tags = []
        combine_files = []

        for job_parameter in job_parameters:
            start_regions = job_parameter[ 0 ]
            insert_res = job_parameter[ 1 ]
            frag_length = job_parameter[ 2 ]

            if len( start_regions ) == 0:
                ##########################################
                # START FROM SCRATCH
                ##########################################
                outfiles = []
                outfile = outdir + '/start_from_scratch.out'

                job_tag = 'REGION_%d_%d_START_FROM_SCRATCH' % (i,j)
                condor_submit_file = 'CONDOR/%s/%s/job.condor' %  (outdir,job_tag)
                fid_dag.write('\nJOB %s %s\n' % (job_tag,condor_submit_file) )

                args2 = '%s -out:file:silent %s -start_from_scratch -in:file:frag_files %s -insert_res %d' % (args, outfile, frag_files[ frag_length ],i )
                make_condor_submit_file( condor_submit_file, args2, 1 )

                job_tags.append( job_tag )
                combine_files.append( outfile )

            else:
                ##################################################
                # APPEND OR PREPEND OR INSERT INTO PREVIOUS PDB
                ##################################################
                infiles = []
                prev_job_tags = []
                for start_region in start_regions:
                    i_prev = start_region[0]
                    j_prev = start_region[1]

                    prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)

                    infile = 'region_%d_%d_sample.cluster.out' % (i_prev,j_prev)
                    infiles.append( infile )

                    prev_job_tags.append( prev_job_tag )

                job_tag = 'REGION_%d_%d_START_FROM_%s_FRAG%d' % (i,j,string.join( prev_job_tags,'_'), frag_length)

                outfile = '%s/start_from_%s_frag%d.out' % (outdir, string.join(prev_job_tags,'_'), frag_length )

                args2 = '%s -out:file:silent %s -in:file:silent %s -insert_res %d -in:file:frag_files %s' % (args, outfile, string.join( infiles ), insert_res, frag_files[ frag_length] )

                if len( start_regions ) > 1: args2 += " -max_input 20 "

                condor_submit_file = 'CONDOR/%s/%s/job.condor' %  (outdir,job_tag)
                fid_dag.write('\nJOB %s %s\n' % (job_tag, condor_submit_file) )

                if not exists( condor_submit_file ):  make_condor_submit_file( condor_submit_file, args2, 1 )

                for prev_job_tag in prev_job_tags:
                    if (prev_job_tag in all_job_tags)  and  (prev_job_tag not in jobs_done): #Note previous job may have been accomplished in a prior run -- not in the current DAG.
                        fid_dag.write('PARENT %s  CHILD %s\n' % (prev_job_tag, job_tag) )

                job_tags.append( job_tag )
                combine_files.append( outfile )

        ##########################################
        # CLUSTER! And keep a small number of representatives (400)
        ##########################################

        if len( combine_files ) == 0: continue

        outfile_cluster = prefix+'sample.cluster.out'
        args_cluster = ' -cluster_test -in:file:silent %s  -database %s -radius %f -out:file:silent %s -nstruct %d -score_diff_cut %8.3f' % (string.join( combine_files ), DB,  CLUSTER_RADIUS, outfile_cluster, FINAL_NUMBER, score_diff_cut )

        condor_submit_cluster_file = 'CONDOR/%s/CLUSTER/job.condor' %  (outdir)
        make_condor_submit_file( condor_submit_cluster_file, args_cluster, 1 )

        fid_dag.write('\nJOB %s %s\n' % (overall_job_tag,condor_submit_cluster_file) )
        fid_dag.write('PARENT %s CHILD %s\n' % (string.join(job_tags),overall_job_tag) )

        all_job_tags.append(  overall_job_tag )
