#!/usr/bin/python

from sys import argv
from os.path import exists,basename,expanduser
from os import system, chdir, getcwd
from parse_options import parse_options

args = argv
no_hb_env_dep = parse_options( args, "no_hb_env_dep", 0 )
near_native = parse_options( args, "near_native", 0 )
loop_force_Nsquared = parse_options( args, "loop_force_Nsquared", 0 )
nstruct = parse_options( args, "nstruct", 400 )

njobs = 412
if near_native: njobs = 100

job_list = args[1]
lines = open( job_list ).readlines()

PUZZLE_DIR = expanduser('~rhiju')+'/projects/loops/mandell/'
puzzle_subdirs = ['plop_set','../difficult','../symm','rosetta_set','../functional','../blind_test2']
CWD = getcwd()

def make_tag( int_vector ):
    tag = ''
    for m in int_vector: tag += ' %d' % m
    return tag

for line in lines:

    tag = line[:-1]
    pdb = tag[:4]

    workdir = tag
    if not exists( workdir ): system( 'mkdir -p '+workdir )

    for subdir in puzzle_subdirs:
        puzzle_dir_test = PUZZLE_DIR + '/' + subdir + '/'
        loop_file = puzzle_dir_test + 'loops/%s.loop' % tag
        if exists( loop_file  ): break

    print loop_file

    assert( exists( loop_file ) )

    pdb_file = puzzle_dir_test + 'start/%s_min.pdb' % pdb
    if not exists( pdb_file ): pdb_file = puzzle_dir_test + 'start/%s_H_stripsidechain.pdb' % pdb
    if not exists( pdb_file ): pdb_file = puzzle_dir_test + 'start/%s.pdb' % pdb
    assert( exists( pdb_file ) )

    if not exists( workdir+'/'+ basename(loop_file) ):
        system( 'rsync '+loop_file+' '+workdir )
        print( 'Copied '+loop_file )

    if not exists( workdir+'/'+ basename(pdb_file) ):
        system( 'rsync '+pdb_file+' '+workdir )
        print( 'Copied '+pdb_file )

    pdb_original_file = puzzle_dir_test + 'start/%s.pdb' % pdb
    if exists( pdb_original_file ) and not exists( workdir+'/'+basename(pdb_original_file ) ):
        system( 'rsync '+pdb_original_file+' '+workdir )
        print( 'Copied '+pdb_original_file )

    disulfide_file = puzzle_dir_test + 'start/%s.disulfides' % pdb
    if exists( disulfide_file ) and not exists( workdir+'/'+basename(disulfide_file ) ):
        system( 'rsync '+disulfide_file+' '+workdir )

    chdir( workdir )

    loop_file = tag+'.loop'
    pdb_file  = pdb+'_min.pdb'
    if not exists( pdb_file ): pdb_file = '%s_H_stripsidechain.pdb' % pdb
    if not exists( pdb_file ): pdb_file = '%s.pdb' % pdb

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
        system( 'generate_CA_constraints.py %s  -cst_res %s -coord_cst -anchor_res 1 -fade > %s ' % (pdb_file,make_tag( range( loop_start, loop_stop+1) ), native_cst_file ) )
        assert( exists( native_cst_file ) )

    prepack_start_file = 'region_%d_%d_sample.cluster.out' % (loop_stop+1, loop_start-1 )
    search_directories = ['//home/vanlang/projects/loops/swa/','/scratch/users/vanlang/projects/loops/swa/','/scratch/users/rhiju/projects/loops/swa_difficult/','/scratch/users/rhiju/projects/loops/swa_difficult_Nsquared','/scratch/users/rhiju/projects/loops/swa/']
    if not exists( prepack_start_file ):
        found_search_dir = 0
        for search_dir in search_directories:
            if exists( search_dir+'/'+tag):
                found_search_dir = 1
                break
        if found_search_dir:
            checkpath = '%s/%s/%s.gz' % (search_dir,workdir,prepack_start_file)
            print checkpath
            if exists( checkpath ):
                print checkpath
                system( 'rsync '+checkpath+' .' )
                system( 'gunzip '+prepack_start_file+'.gz' )
            else:
                checkpath = '%s/%s/%s' % (search_dir,workdir,prepack_start_file)
                if exists( checkpath ):
                    print checkpath
                    system( 'rsync '+checkpath+' .' )

    if near_native and not( exists( prepack_start_file ) ):
        print prepack_start_file, ' missing'
        exit()

    disulfide_file = pdb+'.disulfides'
    readme_setup_file = 'README_SETUP'

    fid = open( readme_setup_file, 'w' )
    fid.write( 'rm -rf STEP* *~ CONDOR core.* SLAVE*  \n' )
    command = 'grinder_dagman.py  -loop_start_pdb %s  -native %s -fasta %s -cluster_radius 0.25 -final_number %d  -denovo 1   -loop_res `seq %d %d` -weights score12.wts -disable_sampling_of_loop_takeoff  ' % (noloop_start_pdb, pdb_file, fasta_file, nstruct, loop_start, loop_stop)
    if ( near_native ): command += ' -rmsd_screen %8.3f -cst_file %s ' % ( 2.0, native_cst_file )
    if ( no_hb_env_dep ): command = command.replace( 'score12.wts', 'score12_no_hb_env_dep.wts' )
    if (loop_force_Nsquared): command += ' -loop_force_Nsquared'
    if (exists( disulfide_file ) ):
        command += ' -disulfide_file '+disulfide_file
        command += ' -exe ~/src/rosetta_TRUNK/rosetta_source/bin/stepwise_protein_test.linuxgccrelease -database ~/src/rosetta_TRUNK/rosetta_database'

    command += '\n'
    fid.write( command )
    fid.close()

    readme_sub_file = 'README_SUB'
    fid = open( readme_sub_file, 'w' )
    fid.write( 'rm -rf blah.* \n' )
    fid.write( 'bsub -W 96:0  -R hname!=node-4-1 -o blah.out -e blah.err SWA_pseudo_dagman_continuous.py -j 400 protein_build.dag \n' )
    fid.close()


    readme_qsub_file = 'README_QSUB'
    fid = open( readme_qsub_file, 'w' )
    fid.write('#!/bin/bash\n')
    fid.write('#PBS -o pseudo_dagman_PBS.out\n')
    fid.write('#PBS -e pseudo_dagman_PBS.err\n')
    fid.write('#PBS -l walltime=48:00:00\n')
    fid.write('\n')
    fid.write('cd $PBS_O_WORKDIR\n')
    fid.write('/home/rhiju/SWA_dagman_python2/SWA_pseudo_dagman_continuous_2.py -j %d protein_build.dag  > blah.out 2> blah.err\n' % njobs)
    fid.close()

    chdir( CWD )
