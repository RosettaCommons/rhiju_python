#!/usr/bin/python

from sys import argv,stdout
import string
from math import sqrt
from parse_options import parse_options

pdbfile = argv[1]

lines = open( pdbfile ).readlines()

oldresnum = ''

fade = parse_options( argv, "fade", 0 )
DIST_CUT = parse_options( argv, "dist_cut", 3.2 )

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
fid = stdout
#fid = open( cst_file, 'w' )
#print "Generating ...  " , cst_file

fid.write( "[ atompairs ]\n" )

STDEV = 0.5
SEQ_SEP_CUTOFF = 2
for i in range( 1,count+1 ):

    assert( i in N_position.keys() )

    for j in range( 1,count+1 ):

        assert( j in O_position.keys() )

        if ( abs( i - j ) <= SEQ_SEP_CUTOFF ): continue
        dist = get_dist( N_position[ i ], O_position[ j ] )
        #print i,j,dist
        if ( dist < DIST_CUT ):
            if fade:
                fid.write( " N %d   O %d  FADE %8.3f %8.3f %8.3f -10.0 10.0 \n" %
                           (i,j, dist - 2*STDEV,
                            dist + 2*STDEV, STDEV ) )
            else:
                fid.write( " N %d   O %d  HARMONIC %8.3f %8.3f \n" % (i,j,dist,STDEV) )
fid.close()
