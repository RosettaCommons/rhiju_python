#!/usr/bin/python

from sys import argv
import string
from math import sqrt

pdbfile = argv[1]


lines = open( pdbfile ).readlines()

oldresnum = ''


N_position = {}
O_position = {}

count = 0
for line in lines:

    if (len(line)>54 and  line[0:4] == 'ATOM' ):

        resnum = line[23:26]
        if not resnum == oldresnum:
            count = count + 1

        atom_name = line[11:15]
        position = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

        if atom_name == '  N ':
            N_position[ count ] = position
        elif atom_name == '  O ':
            O_position[ count ] = position

        oldresnum = resnum

def get_dist( pos1, pos2 ):
    dist2 = 0.0
    for k in range(3):
        dsub = ( pos1[k] - pos2[k] )
        dist2 += dsub*dsub

    return sqrt( dist2 )

cst_file = pdbfile+'.cst'
fid = open( cst_file, 'w' )
print "Generating ...  " , cst_file

fid.write( "[ atompairs ]\n" )

DIST_CUT = 3.1
for i in range( 1,count+1 ):

    assert( i in N_position.keys() )

    for j in range( 1,count+1 ):

        assert( j in O_position.keys() )

        if ( abs( i - j ) <= 1 ): continue
        dist = get_dist( N_position[ i ], O_position[ j ] )
        #print i,j,dist
        if ( dist < DIST_CUT ):
            fid.write( " N %d   O %d  HARMONIC %8.3f 0.5 \n" % (i,j,dist) )

fid.close()
