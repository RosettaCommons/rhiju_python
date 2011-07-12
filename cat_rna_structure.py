#!/usr/bin/python

from sys import argv
from os import popen
from os.path import basename
from glob import glob

indir = argv[ 1 ]

globfiles = glob( '%s/*.seq.ct' % indir )
globfiles.sort()

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
