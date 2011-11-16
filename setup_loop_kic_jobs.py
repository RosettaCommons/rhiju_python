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

nhours = parse_options( argv, 'nhours', 24 );


for line in lines:

    pdb = line[:-1]
    if not exists( pdb ): system( 'mkdir -p '+pdb )

    if not exists( pdb + '/'+pdb+'.loop' ):
        swa_dir_test = SWA_DIR+'/%s/' % pdb

        files = []
        files.append( '%s.loop' % pdb)
        files.append( '%s_min.pdb' % pdb )
        files.append( '%s.fasta' % pdb )
        files.append( 'region_FINAL.out.1.pdb' )

        for file in files:
            assert( exists( swa_dir_test+'/'+file ) )
            if not exists( pdb+'/'+file ):
                system( 'rsync '+swa_dir_test+'/'+file+' '+pdb )
                print( 'Copied '+file )

    chdir( pdb )

    loop_file = pdb+'.loop'
    pdb_file  = pdb+'_min.pdb'

    loop = open( loop_file  ).readlines()[0]
    loop_start = int( loop.split()[0] )
    loop_stop = int( loop.split()[1] )

    readme_file = 'README'
    fid = open( readme_file, 'w' )
    fid.write( '/home/rhiju/src/mini/bin/loopmodel.linuxgccrelease -database /home/rhiju/minirosetta_database  -loops:remodel perturb_kic -loops:refine refine_kic -loops:input_pdb region_FINAL.out.1.pdb -in:file:native %s_min.pdb -loops:loop_file %s.loop -loops:max_kic_build_attempts 10000 -in:file:fullatom -out:file:fullatom -out:prefix 1bhe -out:pdb -ex1 -ex2 -ex1aro -extrachi_cutoff 0 -out:nstruct 1000 -out:file:silent_struct_type binary  -out:file:silent %s_kic.out -fix_ca_bond_angles  -kic_use_linear_chainbreak  -allow_omega_move  -sample_omega_at_pre_prolines\n' % (pdb,pdb,pdb) )
    fid.close()

    system( 'rosetta_submit.py README KIC 20   %d' % nhours )

    chdir( CWD )

    fid_master.write( 'cd %s; source bsubMINI; cd %s\n' % (pdb,CWD ) )

fid_master.close()
print
print
print 'To run jobs, type:'
print ' source %s ' % master_file
