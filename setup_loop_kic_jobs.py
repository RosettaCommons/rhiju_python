#!/usr/bin/python

from sys import argv
from os.path import exists
from os import system, chdir, getcwd
from parse_options import parse_options

job_list = argv[ 1 ]
lines = open( job_list ).readlines()

SWA_DIR = '../swa/'
CWD = getcwd()

master_file = 'bsubMINI'
fid_master = open( master_file, 'w' )

master_file_condor = 'condorMINI'
fid_master_condor = open( master_file_condor, 'w' )

nhours = parse_options( argv, 'nhours', 24 );


EXE_DIR = '/home/rhiju/src/rosetta/rosetta_source/bin/'
if not exists( EXE_DIR ):
    EXE_DIR = '/scratch/scratch95/d/dasr/src/rosetta/rosetta_source/bin/'
assert( exists( EXE_DIR ) )
    
DB = '/home/rhiju/src/rosetta/rosetta_database'
if not exists( DB ):
    DB = '/scratch/scratch95/d/dasr/src/rosetta/rosetta_database'
assert( exists( DB ) )

for line in lines:

    pdb = line[:-1]
    if not exists( pdb ): system( 'mkdir -p '+pdb )

    pdb_file  = pdb+'_min.pdb'

    swa_dir_test = SWA_DIR+'/%s/' % pdb
    if not exists( swa_dir_test+'/'+pdb_file ):
        pdb_file = pdb+'.pdb'
        

    files = []
    files.append( '%s.loop' % pdb)        
    files.append( pdb_file )
    files.append( '%s.fasta' % pdb )
    files.append( 'region_FINAL.out.1.pdb' )

    use_disulfides = False
    if exists( swa_dir_test+'/'+pdb+'.disulfides' ):
        disulfide_file = '%s.disulfides' % pdb 
        files.append( disulfide_file )
        use_disulfides = True
        
    for file in files:
        assert( exists( swa_dir_test+'/'+file ) )
        if not exists( pdb+'/'+file ):
            system( 'rsync '+swa_dir_test+'/'+file+' '+pdb )
            print( 'Copied '+file )


    chdir( pdb )

    loop_file = pdb+'.loop'

    loop = open( loop_file  ).readlines()[0]
    loop_start = int( loop.split()[0] )
    loop_stop = int( loop.split()[1] )

    readme_file = 'README'
    fid = open( readme_file, 'w' )

    
    fid.write( '%s/loopmodel.linuxgccrelease -database %s  -loops:remodel perturb_kic -loops:refine refine_kic -loops:input_pdb region_FINAL.out.1.pdb -in:file:native %s -loops:loop_file %s.loop -loops:max_kic_build_attempts 10000 -in:file:fullatom -out:file:fullatom -out:prefix 1bhe -out:pdb -ex1 -ex2 -ex1aro -extrachi_cutoff 0 -out:nstruct 1000 -out:file:silent_struct_type binary  -out:file:silent %s_kic.out -fix_ca_bond_angles  -kic_use_linear_chainbreak  -allow_omega_move  -sample_omega_at_pre_prolines' % (EXE_DIR,DB,pdb_file,pdb,pdb) )

    if use_disulfides: fid.write( ' -fix_disulf '+disulfide_file )

    fid.write( '\n' )

    fid.close()

    system( 'rosetta_submit.py README KIC 20   %d' % nhours )

    chdir( CWD )

    fid_master.write( 'cd %s; source bsubMINI; cd %s\n' % (pdb,CWD ) )
    fid_master_condor.write( 'cd %s; condor_submit condorMINI; cd %s\n' % (pdb,CWD ) )

fid_master.close()
print
print
print 'To run jobs, type:'
print ' source %s ' % master_file
print '   ... or ...'
print ' source %s ' % master_file_condor
