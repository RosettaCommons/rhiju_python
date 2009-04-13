#!/usr/bin/python

from sys import argv
from os import getcwd, chdir, system
from os.path import basename, dirname, exists,abspath

CWD = getcwd()
pdbfiles = argv[1:]

for file in pdbfiles:
    chdir( dirname( file ) )

    min_file = basename(file)+'.min_pdb'
    if not exists( min_file ):
        command = '~rhiju/python/fix_chains.py -convert '+basename(file)
        print( command )
        system( command )

        command = " ~rhiju/python/charmm_minimize.py  "+ basename( file )
        print( command )
        system( command )

    if not exists( basename(file)+'.rms.txt'  ):
        if ( abspath( file ).count('chunk') ):
            pos = abspath( file ).index( 'chunk')
            native_rna = abspath( file )[ pos: (pos+13)]
            native_rna =  CWD+'/../bench_final/'+native_rna+'_RNA.pdb'
            if exists( native_rna ):
                command = " ~rhiju/python/charmm_superimpose.py "+native_rna+" "+basename(file)
                print command
                system( command )


    if not exists( basename(min_file)+'.rms.txt'  ):
        if ( abspath( file ).count('chunk') ):
            pos = abspath( file ).index( 'chunk')
            native_rna = abspath( file )[ pos: (pos+13)]
            native_rna =  CWD+'/../bench_final/'+native_rna+'_RNA.pdb'
            if exists( native_rna ):
                command = " ~rhiju/python/charmm_superimpose.py "+native_rna+" "+basename(min_file)
                print command
                system( command )

    chdir( CWD )
