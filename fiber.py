#!/usr/bin/python

from sys import argv
from os import system

sequence = argv[1]
if len( argv ) > 2 : tag = argv[2]
else: tag = sequence
command = 'fiber -r -seq=%s %s.pdb' % (sequence,tag)
print command
system( command )

command = 'make_rna_rosetta_ready.py %s.pdb' % tag
print command
system( command )

command = 'rm -f %s.out' % tag
print command
system( command )

command = 'convert_pdb_to_silent_file.macosgccrelease -database ~/rosetta_database -s %s_RNA.pdb -output_silent_file %s.out' % (tag.lower(), tag )
print command
system( command )

