#!/usr/bin/python


from sys import argv
from os import system
from math import sqrt
import string
from parse_options import parse_options

R = parse_options( argv, "R", 0.0 )

infiles = argv[1:]

if R == 0.0:
    cutoff_tag = ''
else:
    cutoff_tag = ' -R %8.3f' % R

command = 'superimpose.py %s %s  > q' % ( string.join( infiles ), cutoff_tag )
#print( command )
system( command )


lines = open( 'q' ).readlines()

CA = []
all_CA = []
all_names = []
count = -1
for i in range( len( lines ) ):
    line = lines[ i ]

    if len( line ) > 6 and line[:5] == 'MODEL':
        if len( CA ) > 0:
            count = count + 1
            all_CA.append( CA )
            all_names.append( infiles[ count ] )
        n = 0
        CA = []
        continue

    if len( line ) > 20 and line[12:16] == ' CA ':
        x = float( line[30:38] )
        y = float( line[38:46] )
        z = float( line[46:54] )
        CA.append( [x,y,z] )


if ( len( CA ) > 0 ):
    count += 1
    all_CA.append( CA )
    all_names.append( infiles[ count ] )


def dist( q, r ):
    return sqrt(  ( q[0] - r[0] ) * ( q[0] - r[0] )  +
                  ( q[1] - r[1] ) * ( q[1] - r[1] ) +
                  ( q[2] - r[2] ) * ( q[2] - r[2] ) )

nres =  len( all_CA[ 0 ] )


all_dists = []
all_names_save = []
for i in range( 1, len( all_CA ) ):
    dists = []
    for n in range( nres ):
        if len( all_CA[i] ) > n:
            dists.append(  dist( all_CA[i][n], all_CA[0][n] ) )

    dists.sort()
    if ( len( dists ) == nres ):
        all_dists.append( dists )
        all_names_save.append( all_names[ i ] )

for n in range( nres ):
    print 100.0 * (n+1.0) / nres ,
    for i in range( len( all_dists ) ):
        print all_dists[ i ][n],
    print

fid = open( 'hubbard_keys.txt','w' )
for i in range( len( all_dists ) ):

    totdev = 0.0
    for n in range( nres ):
        totdev += all_dists[ i ][ n ]
    totdev /= nres
    fid.write( '%s %8.4f\n' %  ( all_names_save[i], totdev ) )
fid.close()
