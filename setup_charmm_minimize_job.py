#!/usr/bin/python

from sys import argv
from os import system, getcwd, chdir
from os.path import exists,abspath
from glob import glob

outfiles = argv[1:]

try:
    NUM_JOBS_PER_NODE = int( outfiles[-1] )
    del( outfiles[ -1 ] )
except:
    NUM_JOBS_PER_NODE = 20


EXE = './run_charmm_minimize.py'

bsub_file = open( 'bsubCHARMM','w' )

CWD = getcwd()

total_count = 0

for outfile in outfiles:

    print outfile

    # Create a directory for extraction
    outdir = outfile.replace( '.out','_OUT')
    if not exists( outdir ):
        command = 'mkdir -p '+outdir
        print( command )
        system( command )


    ####################################################
    # Extract PDB's
    chdir( outdir )
    if ( len(glob( 'S*pdb' )) == 0 and len(glob( 'S*OUT' )) == 0 ) :
        command = 'ln -fs ../'+outfile+' . '
        print( command )
        system( command )

        MINI_EXE = '/work/rhiju/src/mini/bin/rna_extract.linuxgccrelease'
        if not exists( MINI_EXE ):
            MINI_EXE = '~rhiju/src/mini/bin/rna_extract.macosgccrelease'
            if not exists( MINI_EXE ):
                MINI_EXE = '~rhiju/src/mini/bin/rna_extract.linuxgccrelease'

        command = '%s -database ~rhiju/minirosetta_database/ -in::file::silent %s -in::file::silent_struct_type binary_rna' % \
                  ( MINI_EXE, outfile )
        print( command )
        system( command )

        command = 'rm '+outfile # just a soft link anyway
        print( command )
        system( command )

    ####################################################
    # Make subdirectories for each job, and copy in the PDBs, and add line to condor script
    if ( len(glob( 'S*OUT' )) == 0 ) :
        globfiles = glob( 'S_*pdb' )
        for file in globfiles:
            workdir = file.replace( '.pdb', '_OUT' )
            command = 'mkdir -p '+workdir
            print( command )
            system( command )

            command = 'mv '+file+' '+workdir
            print( command )
            system( command )

    count = 0
    start = 0
    globfiles = glob( 'S_*OUT/S*.pdb' )
    globfiles.sort()
    print len( globfiles ),
    if outfile.count( '_native' ) and len( globfiles ) > 1000   : globfiles = globfiles[:1000]
    if outfile.count( '_ideal' ) and len( globfiles ) > 1000    : globfiles = globfiles[:1000]
    if outfile.count( '_nonative' ) and len( globfiles ) > 2000 :   globfiles = globfiles[:2000]
    print len( globfiles )

    for file in globfiles:
        min_pdb_file = file+'.min_pdb'
        if exists( min_pdb_file ):
            continue

        if ( (count % NUM_JOBS_PER_NODE) == 0):
            if ( start == 1 ):
                bsub_file.write( '\n' )
            else:
                start = 1
            bsub_file.write( '\nbsub -W 16:0 %s ' % EXE )
        count += 1
        bsub_file.write( '  '+outdir+'/'+file )

    if (not count % NUM_JOBS_PER_NODE == 0):
        bsub_file.write( '\n')

    chdir( CWD )

    total_count += count

print 'Total number of PDBs to minimize: ', total_count


####################################################
# Create a master script as an "executable" for condor that
# will serially process say, 10 PDBs.
chdir( CWD )
system( 'cp -rf ~rhiju/python/charmm_minimize.py .' )
system( 'cp -rf ~rhiju/python/run_charmm_minimize.py .' )

