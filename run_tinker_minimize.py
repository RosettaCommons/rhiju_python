#!/usr/bin/python

from sys import argv
from os import getcwd, chdir, system
from os.path import basename, dirname, exists,abspath

CWD = getcwd()
pdbfiles = argv[1:]

for file in pdbfiles:
    chdir( dirname( file ) )

    if not exists( "minimize_"+basename(file) ):
        command = " ~rhiju/python/tinker_minimize.py  "+ basename( file )
        print( command )
        system( command )

    if not exists( "minimize_"+basename(file).replace('.pdb','.rms.txt' ) ):
        if ( abspath( file ).count('chunk') ):
            pos = abspath( file ).index( 'chunk')
            native_rna = abspath( file )[ pos: (pos+13)]
            native_rna =  CWD+'/../bench_final/'+native_rna+'_RNA.pdb'
            if exists( native_rna ):
                command = " ~rhiju/python/tinker_superimpose.py -noclean "+native_rna+" minimize_"+basename(file)
                print command
                system( command )


    if not exists( basename(file).replace('.pdb','.rms.txt' ) ):
        if ( abspath( file ).count('chunk') ):
            pos = abspath( file ).index( 'chunk')
            native_rna = abspath( file )[ pos: (pos+13)]
            native_rna =  CWD+'/../bench_final/'+native_rna+'_RNA.pdb'
            if exists( native_rna ):
                command = " ~rhiju/python/tinker_superimpose.py -noclean "+native_rna+" "+basename(file)
                print command
                system( command )

    chdir( CWD )
