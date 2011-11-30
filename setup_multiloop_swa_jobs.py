#!/usr/bin/python

from sys import argv
from os.path import exists
from os import system, chdir, getcwd
from get_sequence import get_sequence
from make_tag import make_tag, make_tag_with_dashes
from parse_options import parse_options
import string

args = argv

loop_force_Nsquared = parse_options( args, "loop_force_Nsquared", 0 )
nstruct = parse_options( args, "nstruct", 1000 )

fasta_file = args[1]
loops_file = args[2]
pdb_file   = args[3]

sequence = open( fasta_file ).readlines()[-1][:-1] # no endline character
print len( sequence )

start_sequence = get_sequence( pdb_file )


# book-keeping -- where are these loops?
in_a_loop = []
for m in range( len( sequence ) ): in_a_loop.append( -1 )

lines = open( loops_file ).readlines()
count = 0
loops = []
loop_jump_edges = []
for line in lines:
    cols = string.split( line )
    loop_start = int( cols[0] )
    loop_stop = int( cols[1] )
    loops.append( [loop_start, loop_stop ] )

    loop_jump_edges.append( loop_start - 1 )
    loop_jump_edges.append( loop_stop   +1 )

    for m in range( loop_start-1, loop_stop ): in_a_loop[ m ] = count
    count += 1

print loops
print in_a_loop

start_sequence_CHECK = ''
in_res = []
for m in range( len(sequence)  ):
    if ( in_a_loop[m] < 0 ) :
        in_res.append( m+1 )
        start_sequence_CHECK += sequence[m]
print start_sequence
print start_sequence_CHECK
in_res_tag = make_tag_with_dashes( in_res )
CWD = getcwd()

# need to prepack structure -- try to be uniform across all cases.

prepack_silent_file = 'global_prepack.out'
if not exists( prepack_silent_file ):
    EXE = '/home/rhiju/src/mini/bin/stepwise_protein_test.linuxgccrelease'
    assert( exists( EXE ) )
    DATABASE =  '/home/rhiju/minirosetta_database'
    assert( exists( DATABASE ) )

    command = ' %s -database %s  -rebuild -out:file:silent_struct_type binary  -fasta %s -n_sample 18 -nstruct 100 -cluster:radius    0.100 -extrachi_cutoff 0 -ex1 -ex2 -score:weights score12.wts -pack_weights pack_no_hb_env_dep.wts -add_peptide_plane -superimpose_res  %s  -fixed_res %s  -calc_rms_res %s  -jump_res  %s  -mute all -s1 %s  -input_res1  %s -use_packer_instead_of_rotamer_trials -out:file:silent %s' % (  EXE, DATABASE, fasta_file, in_res_tag, in_res_tag, in_res_tag, make_tag(loop_jump_edges), pdb_file, in_res_tag, prepack_silent_file )
    print command
    system( command )
    assert( exists( prepack_silent_file ) )

# prepare directories
letters = 'ABCDEFG'
for n in range( len( loops ) ):

    looptag = 'loop'+letters[n]
    if not exists( looptag ): system( 'mkdir -p '+looptag )
    chdir( looptag )

    loop_start = loops[n][0]
    loop_stop  = loops[n][1]

    # need to define global sequence, loop numbering in context of a pdb with other loops cutout!
    sequence_puzzle = ''
    loop_def = []
    for m in range( len(sequence)  ):
        if ( in_a_loop[m] < 0 ) or (in_a_loop[m] == n):  sequence_puzzle += sequence[m]
        if ( in_a_loop[m] == n ): loop_def.append(  len( sequence_puzzle ) )
    loop_start_in_puzzle = loop_def[ 0 ]
    loop_stop_in_puzzle = loop_def[-1 ]

    puzzle_fasta_file = looptag+'.fasta'
    if not exists(  puzzle_fasta_file ):
        fid = open( puzzle_fasta_file, 'w' )
        fid.write( '> Puzzle '+looptag+'\n' )
        fid.write( sequence_puzzle )
        fid.write( '\n' )
        fid.close()

    loop_file = looptag+'.loop' # This actually is not necessary for SWA jobs. Anyway...
    if not exists( looptag + '/'+loop_file ):
        fid = open( loop_file, 'w' )
        fid.write( '%d %d %d 0 1\n' % (loop_start_in_puzzle, loop_stop_in_puzzle, loop_stop_in_puzzle ) )
        fid.close()

    start_pdb_file = pdb_file
    if not exists( start_pdb_file ):
        command = 'rsync ../%s .' % pdb_file
        system( command )
        assert( exists( start_pdb_file ) )

    puzzle_prepack_silent_file = 'region_%d_%d_sample.cluster.out' % (loop_stop+1, loop_start-1 )
    if not exists( puzzle_prepack_silent_file ):
        system( 'rsync ../'+prepack_silent_file+' '+puzzle_prepack_silent_file )
        assert( exists( puzzle_prepack_silent_file ) )

    readme_setup_file = 'README_SETUP'
    if True or not exists( readme_setup_file ):
        fid = open( readme_setup_file, 'w' )
        fid.write( 'rm -rf STEP* *~ CONDOR core.* SLAVE*  \n' )
        command = 'grinder_dagman.py  -loop_start_pdb %s -fasta %s -cluster_radius 0.25 -final_number %d   -denovo 1   -loop_res `seq %d %d` -weights score12.wts ' % (start_pdb_file, puzzle_fasta_file, nstruct, loop_start_in_puzzle, loop_stop_in_puzzle)
        #if ( near_native ): command += ' -rmsd_screen %8.3f -cst_file %s ' % ( 2.0, native_cst_file )
        #if ( no_hb_env_dep ): command = command.replace( 'score12.wts', 'score12_no_hb_env_dep.wts' )
        if (loop_force_Nsquared): command += ' -loop_force_Nsquared'
        command += '\n'

        print command

        fid.write( command )
        fid.close()

    readme_sub_file = 'README_SUB'
    fid = open( readme_sub_file, 'w' )
    #fid.write( 'echo "HELLO WORLD" >> region_FINAL.out \n' )
    fid.write( 'rm -rf blah.* \n' )
    fid.write( 'bsub -W 96:0 -o blah.out -e blah.err SWA_pseudo_dagman_continuous.py -j 400 protein_build.dag \n' )
    fid.close()

    chdir( CWD )
