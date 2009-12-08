#!/usr/bin/python

from sys import argv,stdout
from glob import glob
from string import split
from os.path import basename

globfiles = argv[1:]
#globfiles.sort()

N = {}
recovery = {}

for file in globfiles:

    lines = open( file ).readlines()
    tag = basename( file ).split('.')[0]

    N[ tag ] = {}
    recovery[ tag ] = {}

    kinds_of_residues = []
    for line in lines[-4:]:
        cols = split( line )
        N[ tag ][ cols[0] ] = int( cols[1] )
        recovery[ tag ][ cols[0] ] = float( cols[2] )
        kinds_of_residues.append( cols[0] )

tags = N.keys()
tags.sort()

total_residues = {}
total_recovery = {}
for kind in kinds_of_residues:
    total_residues[ kind ] = 0
    total_recovery[ kind ] = 0

stdout.write( '%10s' % 'RNA')
for kind in kinds_of_residues:
    stdout.write( '  %15s' % kind )
stdout.write( '\n')

for tag in tags:
    stdout.write( '%10s ' % tag)

    for kind in kinds_of_residues:
        total_residues[ kind ] += N[ tag ][ kind]
        total_recovery[ kind ] += N[ tag ][ kind] * recovery[ tag ][ kind]

        stdout.write( '    %3d %8.3f ' % (N[tag][ kind ], recovery[tag][kind]) )
    stdout.write( '\n' )


stdout.write( '%10s ' % 'TOTAL')
for kind in kinds_of_residues:
    total_recovery[ kind ] /= total_residues[kind]
    stdout.write( '    %3d %8.3f ' % \
                  (total_residues[kind], total_recovery[kind]) )

stdout.write( '\n' )
