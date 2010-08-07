#!/usr/bin/python

from sys import argv,exit
from glob import glob
import string
from os import system
from os.path import basename, dirname, exists, expanduser
from time import sleep
from parse_options import parse_options

# By default remove files.
RM_FILES = not(  parse_options( argv, "no_rm_files", 0 ) )

indir_prefix = argv[1]

globstring = indir_prefix+'/*sample.out'
#print globstring
globfiles = glob( globstring )

globstring = indir_prefix+'_S_*/*sample.out'
#print globstring
globfiles2 = glob( globstring )
for file in globfiles2: globfiles.append( file )


if len( globfiles ) == 0:
    sleep( 5 )
    globfiles = glob( globstring )
globfiles.sort()

print globfiles

cat_outfile = dirname( indir_prefix) + '/' + basename(indir_prefix).lower() + '_sample.out'

PYDIR = expanduser('~')+'/python'
assert( exists( PYDIR ) )

command = PYDIR+'/cat_outfiles.py '+string.join( globfiles ) + ' -o ' + cat_outfile
print( command )
system( command )
system( 'chmod 777 '+cat_outfile )

globstring = indir_prefix
globfiles = glob( globstring )
globfiles.sort()

globstring = indir_prefix+'_S_*'
globfiles2 = glob( globstring )
for file in globfiles2: globfiles.append( file )

if RM_FILES:
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
system( 'chmod 777 '+filter_outfile )

##########################################
# DISK SPACE!
if RM_FILES:
    command = 'rm '+cat_outfile
    print( command )
    system( command )
