#!/usr/bin/python

from sys import argv
from glob import glob
import string
from os import system,chdir,getcwd
from os.path import basename, dirname, exists, expanduser
from time import sleep

outfile_cluster = argv[1]
outdir = argv[2]

def wait_for_file( file ):
    for k in range( 10 ):
        if exists( file ):
            break
        else:
            print "waiting for file to show up: ", file
            sleep( 5 )


# Need to be careful to do this WITHIN the directory of interest!
PYDIR = expanduser('~')+'/python'
assert( exists( PYDIR ) )

wait_for_file( outfile_cluster )

command = 'mv '+outfile_cluster+' '+outdir
print( command )
system( command )

chdir( outdir )

wait_for_file( outfile_cluster )

command = PYDIR+'/extract_lowscore_decoys.py -start_at_zero '+outfile_cluster+' 400'
print( command )
system( command )

command = 'mv '+outfile_cluster+' ..'
print( command )
system( command )


