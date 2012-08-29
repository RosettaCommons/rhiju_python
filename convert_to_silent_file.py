#!/usr/bin/python

from sys import argv
from os import system
from os.path import exists,expanduser

EXE = expanduser('~rhiju') + '/src/rosetta_TRUNK/rosetta_source/bin/convert_pdb_to_silent_file.macosgccrelease'

infiles = argv[1:]
for infile in infiles:
    outfile = infile.replace( '.pdb', '.out' )
    if exists( outfile ): system( 'rm -rf '+outfile )
    command = EXE + '  -s %s  -database ~/rosetta_database -output_silent_file %s ' % (infile, outfile )
    print command
    system( command )
