#!/usr/bin/python

from sys import argv
from os.path import exists
from os import system, chdir, getcwd
from parse_options import parse_options

args = argv
no_hb_env_dep = parse_options( args, "no_hb_env_dep", 0 )
near_native = parse_options( args, "near_native", 0 )
loop_force_Nsquared = parse_options( args, "loop_force_Nsquared", 0 )
nstruct = parse_options( args, "nstruct", 1000 )

job_list = args[1]
lines = open( job_list ).readlines()

PUZZLE_DIR = '/home/rhiju/projects/loops/mandell/'
CWD = getcwd()

def make_tag( int_vector ):
    tag = ''
    for m in int_vector: tag += ' %d' % m
    return tag

for line in lines:

    pdb = line[:-1]
    if not exists( pdb ): system( 'mkdir -p '+pdb )

    if not exists( pdb + '/'+pdb+'.loop' ):
        puzzle_dir_test = PUZZLE_DIR + 'plop_set/'
        loop_file = puzzle_dir_test + 'loops/%s.loop' % pdb


        if not exists( loop_file  ):
            puzzle_dir_test = PUZZLE_DIR + 'rosetta_set/'
            loop_file = puzzle_dir_test + 'loops/%s.loop' % pdb

        if not exists( loop_file  ):
            puzzle_dir_test = PUZZLE_DIR + '../functional/'
            loop_file = puzzle_dir_test + 'loops/%s.loop' % pdb

        print loop_file
        assert( exists( loop_file ) )

        pdb_file = puzzle_dir_test + 'start/%s_min.pdb' % pdb
        if not exists( pdb_file ): pdb_file = puzzle_dir_test + 'start/%s_H_stripsidechain.pdb' % pdb
        assert( exists( pdb_file ) )

        system( 'rsync '+loop_file+' '+pdb )
        system( 'rsync '+pdb_file+' '+pdb )

        print( 'Copied '+loop_file )
        print( 'Copied '+pdb_file )

    chdir( pdb )

    loop_file = pdb+'.loop'
    pdb_file  = pdb+'_min.pdb'
    if not exists( pdb_file ): pdb_file = '%s_H_stripsidechain.pdb' % pdb

    noloop_start_pdb = 'noloop_'+pdb_file

    loop = open( loop_file  ).readlines()[0]
    loop_start = int( loop.split()[0] )
    loop_stop = int( loop.split()[1] )

    if not exists( noloop_start_pdb ):
        system( 'renumber_pdb_in_place.py '+pdb_file )

        command = 'pdbslice.py %s -excise %s noloop_' %  ( pdb_file, make_tag( range(loop_start, loop_stop+1) ) )
        system( command )
        assert( exists( noloop_start_pdb ) )

        print( 'Made '+noloop_start_pdb )

    fasta_file = pdb+'.fasta'
    if not exists( fasta_file ):
        system( 'pdb2fasta.py %s > %s ' % ( pdb_file, fasta_file ) )
        assert( exists( fasta_file ) )

    native_cst_file = pdb+'_coordinate2.0.cst'
    if near_native and not exists( native_cst_file ):
        system( 'generate_CA_constraints.py %s_min.pdb  -cst_res %s -coord_cst -anchor_res 1 -fade > %s ' % (pdb,make_tag( range( loop_start, loop_stop+1) ), native_cst_file ) )
        assert( exists( native_cst_file ) )

    prepack_start_file = 'region_%d_%d_sample.cluster.out' % (loop_stop+1, loop_start-1 )
    if not exists( prepack_start_file ) and exists( '/home/vanlang/projects/loops/swa/%s/' % pdb):
        checkpath = '/home/vanlang/projects/loops/swa/%s/%s.gz' % (pdb,prepack_start_file)
        print checkpath
        assert(  exists( checkpath  ) )
        system( 'rsync '+checkpath+' .' )
        system( 'gunzip '+prepack_start_file+'.gz' )
        assert( exists( prepack_start_file ) )

    readme_setup_file = 'README_SETUP'
    if not exists( readme_setup_file ):
        fid = open( readme_setup_file, 'w' )
        fid.write( 'rm -rf STEP* *~ CONDOR core.* SLAVE*  \n' )
        command = 'grinder_dagman.py  -loop_start_pdb %s  -native %s -fasta %s -cluster_radius 0.25 -final_number %d  -denovo 1   -loop_res `seq %d %d` -weights score12.wts -disable_sampling_of_loop_takeoff  ' % (noloop_start_pdb, pdb_file, fasta_file, nstruct, loop_start, loop_stop)
        if ( near_native ): command += ' -rmsd_screen %8.3f -cst_file %s ' % ( 2.0, native_cst_file )
        if ( no_hb_env_dep ): command = command.replace( 'score12.wts', 'score12_no_hb_env_dep.wts' )
        if (loop_force_Nsquared): command += ' -loop_force_Nsquared'

        command += '\n'
        fid.write( command )
        fid.close()

    readme_sub_file = 'README_SUB'
    fid = open( readme_sub_file, 'w' )
    #fid.write( 'echo "HELLO WORLD" >> region_FINAL.out \n' )
    fid.write( 'rm -rf blah.* \n' )
    fid.write( 'bsub -W 96:0 -o blah.out -e blah.err SWA_pseudo_dagman_continuous.py -j 400 protein_build.dag \n' )
    fid.close()

    chdir( CWD )
