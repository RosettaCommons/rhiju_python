#!/usr/bin/python

from sys import argv,exit
from os.path import exists
from os import system, chdir, getcwd
from parse_options import parse_options

job_list = argv[ 1 ]
lines = open( job_list ).readlines()

SWA_DIRS = [ '/scratch/users/rhiju/projects/loops/swa_difficult_Nsquared/', '../swa_difficult/', '../swa/']
CWD = getcwd()

master_file = 'bsubMINI'
fid_master = open( master_file, 'w' )

master_file_condor = 'condorMINI'
fid_master_condor = open( master_file_condor, 'w' )

master_file_qsub = 'qsubMINI'
fid_master_qsub = open( master_file_qsub, 'w' )

nhours = parse_options( argv, 'nhours', 24 );
njobs = parse_options( argv, 'njobs', 20 );
nstruct = parse_options( argv, 'nstruct', 1000 );


EXE_DIR = '/home/rhiju/src/rosetta/rosetta_source/bin/'
if not exists( EXE_DIR ):
    EXE_DIR = '/scratch/scratch95/d/dasr/src/rosetta/rosetta_source/bin/'
    DB = '/scratch/scratch95/d/dasr/src/rosetta/rosetta_database'
if not exists( EXE_DIR ):
    EXE_DIR =  '/home/rhiju/src/rosetta_TRUNK/rosetta_source/bin/'


DB = EXE_DIR + '../../rosetta_database'
assert( exists( EXE_DIR ) )
assert( exists( DB ) )

total_jobs = 0

for line in lines:

    tag = line[:-1]
    pdb = tag[:4]
    if not exists( tag ): system( 'mkdir -p '+tag )

    for swa_dir in SWA_DIRS:
        swa_dir_test = swa_dir+'/'+tag+'/'
        loop_file = swa_dir_test + tag + '.loop'
        if exists( loop_file ): break
    if not exists( loop_file ):
        print 'Could not find: ', loop_file
        exit( 0 )

    pdb_file  = pdb+'_min.pdb'
    if not exists( swa_dir_test+'/'+pdb_file ):
        pdb_file = pdb+'.pdb'

    files = []
    files.append( '%s.loop' % tag)
    files.append( pdb_file )
    files.append( '%s.fasta' % pdb )
    files.append( 'region_FINAL.out.1.pdb' )

    use_disulfides = False
    if exists( swa_dir_test+'/'+pdb+'.disulfides' ):
        disulfide_file = '%s.disulfides' % pdb
        files.append( disulfide_file )
        use_disulfides = True

    for file in files:
        if not exists( swa_dir_test+'/'+file ):
            print 'Could not find: ',  swa_dir_test+'/'+file
            exit()
        if not exists( tag+'/'+file ):
            system( 'rsync '+swa_dir_test+'/'+file+' '+tag )
            print( 'Copied '+file )

    chdir( tag )

    loop_file = tag+'.loop'

    loop = open( loop_file  ).readlines()[0]
    loop_start = int( loop.split()[0] )
    loop_stop = int( loop.split()[1] )

    readme_file = 'README'
    fid = open( readme_file, 'w' )

    fid.write( '%s/loopmodel.linuxgccrelease -database %s  -loops:remodel perturb_kic -loops:refine refine_kic -loops:input_pdb region_FINAL.out.1.pdb -in:file:native %s -loops:loop_file %s -loops:max_kic_build_attempts 10000 -in:file:fullatom -out:file:fullatom -out:prefix 1bhe -out:pdb -ex1 -ex2 -ex1aro -extrachi_cutoff 0 -out:nstruct 1000 -out:file:silent_struct_type binary  -out:file:silent %s_kic.out -fix_ca_bond_angles  -kic_use_linear_chainbreak  -allow_omega_move  -sample_omega_at_pre_prolines' % (EXE_DIR,DB,pdb_file,loop_file,pdb) )

    if use_disulfides: fid.write( ' -fix_disulf '+disulfide_file )

    fid.write( '\n' )

    fid.close()

    system( 'rosetta_submit.py README KIC %d   %d' % (njobs,nhours) )

    total_jobs += njobs

    chdir( CWD )

    fid_master.write( 'cd %s; source bsubMINI; cd %s\n' % (tag,CWD ) )

    fid_master_condor.write( 'cd %s; condor_submit condorMINI; cd %s\n' % (tag,CWD ) )
    fid_master_qsub.write( 'cd %s; source qsubMINI; cd %s\n' % (tag,CWD ) )

fid_master.close()
fid_master_condor.close()
fid_master_qsub.close()

print
print 'Total jobs ready to queue: %d' % total_jobs
print 'To run jobs, type whichever of the following is appropriate:'
print ' source %s ' % master_file
print ' source %s ' % master_file_condor
print ' source %s ' % master_file_qsub
