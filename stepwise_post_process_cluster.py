#!/usr/bin/python

from sys import argv
from glob import glob
import string
from os import system,chdir
from os.path import basename, dirname, exists, expanduser

outfile_cluster = argv[1]
outdir = argv[2]

# Need to be careful to do this WITHIN the directory of interest!
PYDIR = expanduser('~')+'/python'
assert( exists( PYDIR ) )

command = 'mv '+outfile_cluster+' '+outdir
print( command )
system( command )

chdir( outdir )

command = PYDIR+'/extract_lowscore_decoys.py -start_at_zero '+outfile_cluster+' 400'
print( command )
system( command )

command = 'mv '+outfile_cluster+' ..'
print( command )
system( command )


