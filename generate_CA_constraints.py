#!/usr/bin/python

from sys import argv,stdout
import string
from math import sqrt
from parse_options import parse_options

args = argv

pdbfile = args[1]

fade = parse_options( argv, "fade", 0 )
DIST_CUT = parse_options( argv, "dist_cut", 9.0 )
STDEV = parse_options( argv, "stdev", 2.0 )
fixed_res = parse_options( argv, "fixed_res", [-1] )


lines = open( pdbfile ).readlines()

oldresnum = ''

CA_position = {}

count = 0
for line in lines:

    if (len(line)>54 and  line[0:4] == 'ATOM' ):

        resnum = line[23:26]
        if not resnum == oldresnum:
            count = count + 1

        atom_name = line[11:15]
        position = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

        if atom_name == '  CA':
            CA_position[ count ] = position

        oldresnum = resnum

def get_dist( pos1, pos2 ):
    dist2 = 0.0
    for k in range(3):
        dsub = ( pos1[k] - pos2[k] )
        dist2 += dsub*dsub

    return sqrt( dist2 )

cst_file = pdbfile+'.cst'
fid = stdout
#print "Generating ...  " , cst_file

fid.write( "[ atompairs ]\n" )

for i in range( 1,count+1 ):

    for j in range( i,count+1 ):

        if ( abs( i - j ) <= 3 ): continue

        if len( fixed_res ) > 0 and (i in fixed_res) and (j in fixed_res): continue

        dist = get_dist( CA_position[ i ], CA_position[ j ] )
        #print i,j,dist
        if ( dist < DIST_CUT ):
            if fade:
                fid.write( " CA %d   CA %d  FADE %8.3f %8.3f %8.3f %8.3f %8.3f\n" % \
                               (i,j,\
                                    (dist-2*STDEV), (dist+2*STDEV), STDEV, -10.0, 10.0  ) )
            else:
                fid.write( " CA %d   CA %d  HARMONIC %8.3f %8.3f \n" % (i,j,dist,STDEV) )

fid.close()
