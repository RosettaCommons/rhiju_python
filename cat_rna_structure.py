#!/usr/bin/python

from sys import argv,exit
from os import popen, system
from os.path import basename,dirname
from glob import glob

indir = argv[ 1 ]

globfiles = glob( '%s/*/bpp.txt' % indir )
globfiles.sort()
for bpp_file in globfiles:
    new_file = basename( dirname( bpp_file) ).replace( '_PARTITION', '_partition_bpp.txt' )
    command = 'cp %s %s' % (bpp_file, new_file )
    print command
    system( command )

globfiles = glob( '%s/*.seq.ct' % indir )
globfiles.sort()

#exit( 0 )

for ct_file in globfiles:
    cat_file = basename( ct_file ).replace('.seq.ct','_structures.txt' )
    fid = open( cat_file, 'w' )
    lines = popen( 'ct_to_structure.py %s' % ct_file ).readlines()
    fid.write( lines[0] )

    bootfiles = glob( '%s.boot*' %  ct_file )
    bootfiles.sort()
    for bootfile in bootfiles:
        lines = popen( 'ct_to_structure.py %s' % bootfile ).readlines()
        fid.write( lines[0] )

    print 'Outputting data (including %d bootstraps) to: %s' % (len( bootfiles), cat_file )
    fid.close()
