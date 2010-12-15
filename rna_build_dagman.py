#!/usr/bin/python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
from parse_options import parse_options

# Define sequence
fasta_file = parse_options( argv, "fasta", "1shf.fasta" )
assert( exists( fasta_file ) )
sequence = open( fasta_file  ).readlines()[1][:-1]

# Defines loop residues.
MIN_RES = parse_options( argv, "min_res", 1 )
MAX_RES = parse_options( argv, "max_res", len( sequence ) )

# Starting PDB files.
native_pdb = parse_options( argv, "native", "" )
input_pdbs = parse_options( argv, "s", ["1shf.pdb"] )
input_res_full = parse_options( argv, "input_res", [ 0 ] )
score_diff_cut = parse_options( argv, "score_diff_cut", 1000000.0 )
CLUSTER_RADIUS = parse_options( argv, "cluster_radius", 0.25 )
CLUSTER_RADIUS_SAMPLE = parse_options( argv, "cluster_radius_sample", 0.1 )
AUTO_TUNE = parse_options( argv, "auto_tune", 0 )
NSTRUCT = parse_options( argv, "nstruct", 400 )
cutpoints_open = parse_options( argv, "cutpoint_open", [ -1 ] )
fixed_res = parse_options( argv, "fixed_res", [ -1 ] )
superimpose_res = parse_options( argv, "superimpose_res", [ -1 ] )
calc_rms_res = parse_options( argv, "calc_rms_res", [ -1 ] )
terminal_res = parse_options( argv, "terminal_res", [ -1 ] )
bulge_res = parse_options( argv, "bulge_res", [ -1 ] )
internal_loop = parse_options( argv, "internal_loop", 0 )
native_rmsd_screen =  parse_options( argv, "native_rmsd_screen", 0 )
max_slip =  parse_options( argv, "max_slip", -1 )
cst_file = parse_options( argv, "cst_file", "" )
score_weights = parse_options( argv, "score_weights", "" )
no_rm_files= parse_options( argv, "no_rm_files", 0 )

if ( len( argv ) > 1 ): # Should remain with just the first element, the name of this script.
    print " Unrecognized flags?"
    print "   ",string.join(argv[1:] )
    exit( 0 )

assert( len( input_pdbs ) > 0 ) # Later can create a mode that builds from scratch
# Assert uniqueness -- not overlapping input pdbs!
for i in range( len(input_res_full) ): assert( input_res_full[i] not in input_res_full[:i] )



###############################################################
# Where's the executable?
###############################################################
HOMEDIR = expanduser('~rhiju')

EXE = HOMEDIR+'/src/mini/bin/rna_swa_test.macosgccrelease'
if not( exists( EXE )):
    EXE = HOMEDIR+'/src/mini/bin/rna_swa_test.linuxgccrelease'
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

if len( score_weights ) > 0 : assert( exists( score_weights ) or exists( DB + "/scoring/weights/"+score_weights) )

if no_rm_files:
    POST_PROCESS_CLUSTER_SCRIPT += ' -no_rm_files 1'
    POST_PROCESS_FILTER_SCRIPT += ' -no_rm_files 1'

fid_dag = open( "rna_build.dag", 'w' )
fid_dag.write("DOT dag.dot\n")

#############################################
# Should probably put in a check to make sure
# that starting pdb ("input_pdb") has the correct
# sequence!
#############################################
count = 0
#print input_pdbs
input_res  = []
for input_pdb in input_pdbs:

    assert( exists( input_pdb ) )

    input_pdb_sequence = popen( '%s/pdb2fasta.py %s ' % (PYDIR,input_pdb) ).readlines()[1][:-1]
    input_seq_length = len( input_pdb_sequence )

    starting_sequence = ""
    input_res_particular = []
    for i in range( count, count+input_seq_length ):
        starting_sequence += sequence[ input_res_full[ i ] - 1  ]
        input_res_particular.append( input_res_full[ i ] )
    input_res.append( input_res_particular )

    print starting_sequence, input_pdb_sequence
    assert( starting_sequence == input_pdb_sequence )

    count += input_seq_length


###############################################
# Now define the elements that need to be rebuilt.
#  -- most are going to be single residues
#  -- some are going to be input_pdbs though.
#
# "coloring" --> -1 if free residue, otherwise 0,1,2,... if assigned to an input pdb.
#   I'm not entirely sure if we need to go through this coloring step...
coloring = []
for i in range( len( sequence ) ):
    color = -1
    for j in range( len( input_res ) ):
        if (i+1) in input_res[j]:
            color = j
    coloring.append( color )

#print coloring
# for helix cap case:
#
# rna_build_dagman.py  -native helix_full.out.1.pdb -fasta helix_full.fasta  -s helix_half.out.1.pdb au.pdb -input_res 1 2 7 8 4 5 -cutpoint_open 4
#
# coloring   is        0, 0, -1, 1, 1, -1, 0, 0
# assigned_element is  0, 0,  1, 2, 2,  3, 0, 0
#
#  i.e., there are actually four moving elements, numbered 0, 1, 2, and 3
#   for assigned_element, its convenient to start numbering at zero because of
#   varius modulo operations coming up...
#
colors_so_far = []
num_elements = -1
color_to_element = {}
element_to_color = {}
element_definition = {}
assigned_element = []
for i in range( len(sequence ) ):

    if ( coloring[ i ] == -1 ):
        num_elements += 1
        assign_element = num_elements
    elif  (coloring[ i ]  not in colors_so_far ):
        num_elements += 1
        colors_so_far.append( coloring[ i ] )
        color_to_element[ coloring[ i ] ] = num_elements
        element_to_color[ num_elements ] = coloring[ i ]
        assign_element = num_elements
    else:
        assign_element = color_to_element[ coloring[ i ] ]

    if assign_element not in element_definition:
        element_definition[ assign_element ] = []
    element_definition[ assign_element ].append( i+1  )

    assigned_element.append( assign_element )

num_elements += 1

#print color_to_element
#print element_to_color

print "#################################"
print "       DEFINITION OF ELEMENTS    "
print "#################################"
print element_definition
print " or: "
print assigned_element
print

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
real_compute_job_tags = []
jobs_done = []
last_jobs = []
last_outfiles = []

# Tags, etc. are now indexed by "element", 0, 1, 2, 3, etc...
input_file_tags = []
for input_num in range( len( input_pdbs ) ):
    element_assignment_for_input_pdb = color_to_element[ input_num ]
    input_file_tag = 'REGION_%d_%d' % (element_assignment_for_input_pdb, element_assignment_for_input_pdb)
    all_job_tags.append( input_file_tag )
    jobs_done.append( input_file_tag )
    input_file_tags.append( input_file_tag )

def get_modeled_elements( i, j ):
    modeled_elements = []
    count = i
    modeled_elements.append( count )
    while not( count == j ):
        count += 1
        count = count % num_elements
        modeled_elements.append( count )
    return modeled_elements

def get_boundary_res( i, assigned_element ):
    sample_res = -1
    seq_length = len( assigned_element )
    num_elements = max( assigned_element ) + 1
    target_element = i % num_elements
    target_element_next = (i+1) % num_elements
    for k in range( seq_length  ):
        if (assigned_element[k] == target_element) and \
               (assigned_element[ (k+1) % seq_length ] == target_element_next ):
            sample_res = k
            break
    return (sample_res+1)


# BASIC COMMAND
args = ' -algorithm rna_resample_test -database %s -fasta %s -output_virtual  ' %  ( DB, fasta_file  )
args +=  ' -cluster:radius %8.3f ' % CLUSTER_RADIUS_SAMPLE
args +=  ' -num_pose_kept %d ' % NSTRUCT

if len( score_weights ) > 0: args += ' -score:weights %s ' % score_weights

if len( native_pdb ) > 0: args += '-native %s ' % native_pdb


if len( cst_file ) > 0:
    assert( exists( cst_file ) )
    args += ' -cst_file %s' % cst_file

if len( cutpoints_open ) > 0:
    args += ' -cutpoint_open '
    for cutpos in cutpoints_open: args += '%d ' % cutpos

if len( fixed_res ) > 0:
    args += ' -fixed_res '
    for k in fixed_res: args += '%d ' % k

if len( superimpose_res ) == 0: superimpose_res = fixed_res

if len( superimpose_res ) > 0:
    args += ' -superimpose_res '
    for k in superimpose_res: args += '%d ' % k

if len( terminal_res ) > 0:
    args += ' -terminal_res '
    for k in terminal_res: args += '%d ' % k

if len( bulge_res ) > 0:
    args += ' -bulge_res '
    for k in bulge_res: args += '%d ' % k

if native_rmsd_screen: args += ' -sampler_native_rmsd_screen '

if AUTO_TUNE:
    cluster_tag = ' -auto_tune '
else:
    cluster_tag = ' -cluster:radius %s ' % CLUSTER_RADIUS

if len( calc_rms_res ) > 0:
    cluster_tag += ' -calc_rms_res'
    for k in calc_rms_res: cluster_tag += ' %d' % k

# Order calculation based on number of elements modeled -- smaller fragments first.
for L in range( 2, num_elements+1 ):

    for k in range( num_elements ) :

        i = k
        j = ( k + L - 1 ) % num_elements

        slip = abs( i + j ) % num_elements
        if abs(slip-num_elements) < abs(slip): slip = slip - num_elements
        if max_slip > 0 and ( abs(slip) > max_slip ) : continue

        # Native PDB.
        prefix = 'region_%d_%d_' % (i,j)

        modeled_elements = get_modeled_elements( i, j )

        # SPECIAL --> For now, I only want to begin with an input pdb.
        #if ( 0 not in modeled_elements ): continue

        # This job is maybe already done...
        outfile_cluster = prefix+'sample.cluster.out'
        overall_job_tag = 'REGION_%d_%d' % (i,j)

        if exists( outfile_cluster ):
            all_job_tags.append(  overall_job_tag )
            jobs_done.append( overall_job_tag   )

            if ( L == num_elements ):
                last_jobs.append( overall_job_tag )
                last_outfiles.append( outfile_cluster )

            continue

        ###########################################
        # OUTPUT DIRECTORY
        outdir = 'REGION_%d_%d' % (i,j)

        ###########################################
        # DO THE JOBS
        start_regions = []

        i_prev = (i + 1) % num_elements
        j_prev = j
        prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
        if (prev_job_tag in all_job_tags) and \
           ( max( element_definition[i]) not in cutpoints_open) and \
           ( (not internal_loop) or (i not in element_to_color.keys()) or L == num_elements ) :
            start_regions.append( [i_prev, j_prev ] )

        i_prev = i
        j_prev = (j - 1) % num_elements
        prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
        if (prev_job_tag in all_job_tags) and \
           ( max( element_definition[j_prev])  not in cutpoints_open) and \
           ( (not internal_loop) or (j not in element_to_color.keys()) or L == num_elements ) :
            start_regions.append( [i_prev, j_prev ] )

        job_tags = []
        combine_files = []

        if len( start_regions ) > 0 and not exists( outdir ): system( 'mkdir -p '+outdir )

        for start_region in start_regions:
            i_prev = start_region[0]
            j_prev = start_region[1]

            dir_prev = 'REGION_%d_%d' % (i_prev, j_prev )

            # prev_job_tag and dir_prev are currently the same, but may change directory structure in the future.
            prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)

            num_jobs_to_queue = NSTRUCT

            ######################################################
            # Start with "fixed", pre-constructed pose
            ######################################################
            # Job input depends on whether we are starting from pdb, silent file, etc.
            if prev_job_tag in input_file_tags:
                input_num = input_file_tags.index( prev_job_tag )
                newdir = outdir+'/START_FROM_REGION_%d_%d' % (i_prev, j_prev )
                if not exists( newdir ): system( 'mkdir -p '+newdir )
                outfile = newdir + '/' + prefix + 'sample.out'
                args2 = args

                args2 += " -s1 %s " % input_pdbs[ input_num ]
                args2 += ' -input_res1 '
                for k in input_res[ input_num ]: args2 += '%d ' % k

                num_jobs_to_queue = 1
            else:

                modeled_elements_prev = get_modeled_elements( i_prev, j_prev )
                #if the fixed region is a single element, for now it needs to be read in from disk -- an input file, covered above.
                if len( modeled_elements_prev ) == 1:
                    continue

                infile = 'region_%d_%d_sample.cluster.out' % (i_prev,j_prev)
                tag = 'S_$(Process)'
                args2 = '%s -in:file:silent_struct_type binary_rna -silent1 %s -tags1 %s ' % (args, infile, tag )
                newdir = outdir+'/START_FROM_REGION_%d_%d_%s' % (i_prev, j_prev, tag.upper() )
                outfile = newdir + '/' + prefix + 'sample.out'

                # What is already built? What will move? I think the code guarantees that any silent file
                # will have residues modeled in appropriate sequential order. I hope.
                modeled_res_prev = []
                for m in modeled_elements_prev:
                    for k in element_definition[ m ]: modeled_res_prev.append( k )
                modeled_res_prev.sort()
                # This defines input res for only the first pdb/silent_struct.
                args2 += ' -input_res1 '
                for k in modeled_res_prev: args2 += '%d ' % k


            args2 += '-out:file:silent %s ' % outfile


            ######################################################
            # Then define the moving element.
            ######################################################
            if ( i == i_prev ):
                moving_element = j
            else:
                moving_element = i

            if len( element_definition[ moving_element ] ) == 1: # Just moving a single residue, rebuild it from scratch.
                args2 += ' -sample_res %d ' % element_definition[ moving_element ][0]
            else:
                # Its an input file.
                input_num = element_to_color[ moving_element ]
                args2 += " -s2 %s " % input_pdbs[ input_num ]
                args2 += ' -input_res2 '
                for k in input_res[ input_num ]: args2 += '%d ' % k

                sample_res = -1
                if ( i == i_prev ):
                    sample_res = get_boundary_res( j-1, assigned_element )
                else:
                    sample_res = get_boundary_res( i, assigned_element )
                assert( sample_res > 0 )
                args2 += ' -sample_res %d ' % sample_res


            # Special --> close cutpoint.
            if ( L == num_elements ):
                boundary_res = get_boundary_res( j, assigned_element )
                if ( boundary_res not in cutpoints_open ):
                    args2 += ' -cutpoint_closed %d ' % boundary_res

            job_tag = 'REGION_%d_%d_START_FROM_REGION_%d_%d' % (i,j,i_prev,j_prev)
            condor_submit_file = 'CONDOR/%s.condor' %  job_tag
            fid_dag.write('\nJOB %s %s\n' % (job_tag, condor_submit_file) )

            if not exists( condor_submit_file ):
                make_condor_submit_file( condor_submit_file, args2, num_jobs_to_queue )

            if (prev_job_tag in all_job_tags)  and   (prev_job_tag not in jobs_done): #Note previous job may have been accomplished in a prior run -- not in the current DAG.
                fid_dag.write('PARENT %s  CHILD %s\n' % (prev_job_tag, job_tag) )

            # The pre process script finds out how many jobs there actually are...
            if ( prev_job_tag not in input_file_tags ):
                fid_dag.write('SCRIPT PRE %s   %s %s %s %s\n' % (job_tag, PRE_PROCESS_SETUP_SCRIPT,outdir,dir_prev,condor_submit_file) )
            fid_dag.write('SCRIPT POST %s %s %s/START_FROM_REGION_%d_%d\n' % (job_tag, POST_PROCESS_FILTER_SCRIPT,outdir,i_prev,j_prev ) )

            job_tags.append( job_tag )
            real_compute_job_tags.append( job_tag )
            combine_files.append( '%s/start_from_region_%d_%d_sample.low4000.out' % ( outdir, i_prev,j_prev) )


        ##########################################
        # CLUSTER! And keep a small number of representatives (400)
        ##########################################

        if len( combine_files ) == 0: continue

        print 'Modeling elements: ', modeled_elements

        outfile_cluster = prefix+'sample.cluster.out'
        args_cluster = ' -algorithm cluster_old -in:file:silent %s  -in:file:silent_struct_type binary_rna  -database %s  %s -out:file:silent %s  -score_diff_cut %8.3f -silent_read_through_errors  -nstruct %d ' % (string.join( combine_files ), DB,  cluster_tag , outfile_cluster, score_diff_cut , NSTRUCT )


        condor_submit_cluster_file = 'CONDOR/REGION_%d_%d_cluster.condor' % (i,j)
        make_condor_submit_file( condor_submit_cluster_file, args_cluster, 1 )

        fid_dag.write('\nJOB %s %s\n' % (overall_job_tag,condor_submit_cluster_file) )
        fid_dag.write('PARENT %s CHILD %s\n' % (string.join(job_tags),overall_job_tag) )
        fid_dag.write('SCRIPT POST %s %s %s %s\n' % (overall_job_tag, POST_PROCESS_CLUSTER_SCRIPT, outfile_cluster, overall_job_tag ) )

        all_job_tags.append(  overall_job_tag )

        if ( L == num_elements ):
            last_jobs.append( overall_job_tag )
            last_outfiles.append( outfile_cluster )


#####################################################################################
# New loop -- "circularize" two segments:
#
#     Segment 1. i+1 to j
#     Segment 2. j+1 to i
#
# the combination closes to model the whole molecule.
#  In the future, we could actually include i to k   and   k to j   (N^3)
#  In that case, we should unify with the main loop above.

# Order calculation based on number of elements modeled -- smaller fragments first.
CLOSE_CIRCLE = 0
if CLOSE_CIRCLE:

    for i in range( num_elements ): # sample res
        for j in range( num_elements ) : # cutpoint_closed

            if ( i == j ): continue

            # Native PDB.
            prefix = 'region_%d_close_%d_' % (i,j)

            modeled_elements1 = get_modeled_elements( (i+1) % num_elements, j )
            modeled_elements2 = get_modeled_elements( (j+1) % num_elements, i )

            # these cases are already covered in main loop.
            if ( len( modeled_elements1 ) == 1 ): continue
            if ( len( modeled_elements2 ) == 1 ): continue

            print 'Modeling: circle close ',i,j

            # This job is maybe already done...
            outfile_cluster = prefix+'sample.cluster.out'
            overall_job_tag = 'REGION_%d_CLOSE_%d' % (i,j)

            if exists( outfile_cluster ):
                all_job_tags.append(  overall_job_tag )
                jobs_done.append( overall_job_tag   )

                last_jobs.append( overall_job_tag )
                last_outfiles.append( outfile_cluster )
                continue

            ###########################################
            # OUTPUT DIRECTORY
            outdir = overall_job_tag
            combine_files = []
            job_tags = []
            ###########################################
            # DO THE JOB
            if not exists( outdir ): system( 'mkdir -p '+outdir )
            prev_job_tag1 = 'REGION_%d_%d' % ( (i+1) % num_elements ,j)
            prev_job_tag2 = 'REGION_%d_%d' % ( (j+1) % num_elements ,i)

            num_jobs_to_queue = NSTRUCT

            infile1 = 'region_%d_%d_sample.cluster.out' % ( (i+1) % num_elements, j)
            infile2 = 'region_%d_%d_sample.cluster.out' % ( (j+1) % num_elements, i)

            tag = '$(Process)'
            args2 = '%s -in:file:silent_struct_type binary_rna -in:file:silent %s %s -combo %s ' % (args, infile1, infile2, tag )
            newdir = outdir+'/REGION_%d_CLOSE_%d_COMBO_%s' % (i, j, tag )

            #tag = 'S_$(Process)'
            #args2 = '%s -in:file:silent_struct_type binary_rna -in:file:silent %s %s -tags %s S_0 ' % (args, infile1, infile2, tag )
            #newdir = outdir+'/START_FROM_REGION_%d_%d_%s' % ( (i+1) % num_elements, j, tag )

            outfile = newdir + '/' + prefix + 'sample.out'
            for m in range( NSTRUCT ):
                actual_dir = newdir.replace( tag, '%d' % m )
                if not exists( actual_dir ): system( 'mkdir -p '+actual_dir )

            # What is already built? What will move? I think the code guarantees that any silent file
            # will have residues modeled in appropriate sequential order. I hope.
            modeled_res1 = []
            args2 += ' -input_res1 '
            for m in modeled_elements1:
                for k in element_definition[ m ]: modeled_res1.append( k )
                modeled_res1.sort()
            for k in modeled_res1: args2 += '%d ' % k

            modeled_res2 = []
            args2 += ' -input_res2 '
            for m in modeled_elements2:
                for k in element_definition[ m ]: modeled_res2.append( k )
                modeled_res2.sort()
            for k in modeled_res2: args2 += '%d ' % k

            args2 += '-out:file:silent %s ' % outfile

            ######################################################
            # Then define the moving element.
            ######################################################
            sample_res = get_boundary_res( i, assigned_element )
            args2 += ' -sample_res %d ' % sample_res

            boundary_res = get_boundary_res( j, assigned_element )
            if ( boundary_res not in cutpoints_open ):
                args2 += ' -cutpoint_closed %d ' % boundary_res

            job_tag = 'START_FROM_REGION_%d_CLOSE_%d' % (i,j)
            condor_submit_file = 'CONDOR/%s.condor' %  job_tag
            fid_dag.write('\nJOB %s %s\n' % (job_tag, condor_submit_file) )

            if not exists( condor_submit_file ):
                make_condor_submit_file( condor_submit_file, args2, num_jobs_to_queue )

            if (prev_job_tag1 in all_job_tags)  and   (prev_job_tag1 not in jobs_done): #Note previous job may have been accomplished in a prior run -- not in the current DAG.
                fid_dag.write('PARENT %s  CHILD %s\n' % (prev_job_tag1, job_tag) )
            if (prev_job_tag2 in all_job_tags)  and   (prev_job_tag2 not in jobs_done): #Note previous job may have been accomplished in a prior run -- not in the current DAG.
                fid_dag.write('PARENT %s  CHILD %s\n' % (prev_job_tag2, job_tag) )

            #if ( prev_job_tag not in input_file_tags ):
            #fid_dag.write('SCRIPT PRE %s   %s %s %s %s\n' % (job_tag, PRE_PROCESS_SETUP_SCRIPT,outdir,prev_job_tag1,condor_submit_file) )
            fid_dag.write('SCRIPT POST %s %s %s/REGION_%d_CLOSE_%d\n' % (job_tag, POST_PROCESS_FILTER_SCRIPT,outdir,i,j ) )

            job_tags.append( job_tag )
            real_compute_job_tags.append( job_tag )
            combine_files.append( '%s/region_%d_close_%d_sample.low4000.out' % ( outdir, i,j) )

            ##########################################
            # CLUSTER! And keep a small number of representatives (400)
            ##########################################

            if len( combine_files ) == 0: continue

            outfile_cluster = prefix+'sample.cluster.out'
            args_cluster = ' -algorithm cluster_old -in:file:silent %s  -in:file:silent_struct_type binary_rna  -database %s  %s -out:file:silent %s  -score_diff_cut %8.3f -silent_read_through_errors  -nstruct %d ' % (string.join( combine_files ), DB,  cluster_tag , outfile_cluster, score_diff_cut , NSTRUCT )


            condor_submit_cluster_file = 'CONDOR/REGION_%d_CLOSE_%d_cluster.condor' % (i,j)
            make_condor_submit_file( condor_submit_cluster_file, args_cluster, 1 )

            fid_dag.write('\nJOB %s %s\n' % (overall_job_tag,condor_submit_cluster_file) )
            fid_dag.write('PARENT %s CHILD %s\n' % ( string.join( job_tags ),overall_job_tag) )
            all_job_tags.append(  overall_job_tag )

            last_jobs.append( overall_job_tag )
            last_outfiles.append( outfile_cluster )


#####################################################################################
final_outfile = "region_FINAL.out"
if not exists( final_outfile ):
    args_cluster = ' -algorithm cluster_old -in:file:silent %s  -in:file:silent_struct_type binary_rna  -database %s  %s -out:file:silent %s  -score_diff_cut %8.3f -silent_read_through_errors  -nstruct %d ' % (string.join( last_outfiles ), DB,  cluster_tag, final_outfile, 2 * score_diff_cut, 10000 )

    condor_submit_cluster_file = 'CONDOR/REGION_FINAL_cluster.condor'
    make_condor_submit_file( condor_submit_cluster_file, args_cluster, 1 )

    final_job_tag = "REGION_FINAL"
    fid_dag.write('\nJOB %s %s\n' % ( final_job_tag,condor_submit_cluster_file) )
    for prev_job_tag in last_jobs:
        if ( prev_job_tag not in jobs_done ):
            fid_dag.write('PARENT %s  CHILD %s\n' % (prev_job_tag, final_job_tag) )
    #fid_dag.write('SCRIPT POST %s %s %s %s\n' % (overall_job_tag, POST_PROCESS_CLUSTER_SCRIPT, outfile_cluster, outdir ) )


print
print "Total number of jobs to run (not counting clustering):", len( real_compute_job_tags )
print "Total number of final outfiles (and clustering jobs):", len( all_job_tags )


