#!/usr/bin/python

from sys import argv,exit
from glob import glob
import string
from os import system
from os.path import basename, dirname, exists, expanduser
from time import sleep

indir_prefix = argv[1]

globstring = indir_prefix+'*/*sample.out'
print globstring

globfiles = glob( globstring )
if len( globfiles ) == 0:
    sleep( 5 )
    globfiles = glob( globstring )
globfiles.sort()

#print globfiles

cat_outfile = dirname( indir_prefix) + '/' + basename(indir_prefix).lower() + '_sample.out'

PYDIR = expanduser('~')+'/python'
assert( exists( PYDIR ) )

command = PYDIR+'/cat_outfiles.py '+string.join( globfiles ) + ' -o ' + cat_outfile
print( command )
system( command )


globstring = indir_prefix+'*'
globfiles = glob( globstring )
globfiles.sort()
for globfile in globfiles:
    command = 'rm -rf '+globfile
    print( command )
    system( command )

##########################################
##########################################
# We don't need to whole file -- just
#  some small fraction (here 4000).
filter_outfile = cat_outfile.replace('.out','.low4000.out')

#wait_for_file( cat_outfile )

command = PYDIR+'/extract_lowscore_decoys_outfile.py '+cat_outfile+' 4000 > '+filter_outfile
print( command )
system( command )


##########################################
# DISK SPACE!
command = 'rm '+cat_outfile
print( command )
system( command )
