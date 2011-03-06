#!/usr/bin/python

from sys import argv
from math import sqrt, acos, pi

file = argv[ 1 ]

coord = {}
resnums = []

backbone_atoms = [ ' N  ',' CA ',' C  '];
for atom in backbone_atoms: coord[ atom ] = {}

lines = open( file ).readlines()

for i in range( len( lines ) ):
    line = lines[ i ]

    if len( line ) > 20 and line[12:16] in backbone_atoms:
        atom = line[12:16]
        x = float( line[30:38] )
        y = float( line[38:46] )
        z = float( line[46:54] )
        resnum = int( line[22:26] )
        if resnum not in resnums: resnums.append( resnum )

        coord[ atom ][ resnum ] = (x,y,z)

def get_dist( coord, i1,j1,i2,j2 ):
    dist = 0.0
    if j1 in coord.keys() and  i1 in coord[ j1 ].keys() \
            and j2 in coord.keys() and  i2 in coord[ j2 ].keys():
        r1 = coord[ j1 ][ i1 ]
        r2 = coord[ j2 ][ i2 ]
        dist2 = ( r1[0] - r2[0] ) * ( r1[0] - r2[0] ) + \
            ( r1[1] - r2[1] ) * ( r1[1] - r2[1] ) + \
            ( r1[2] - r2[2] ) * ( r1[2] - r2[2] )
        dist = sqrt( dist2 )
    return dist

def get_angle( coord, i1,j1,i2,j2,i3,j3 ):
    angle = 0.0

    if  j1 in coord.keys() and  i1 in coord[ j1 ].keys() \
            and j2 in coord.keys() and  i2 in coord[ j2 ].keys() \
            and j3 in coord.keys() and  i3 in coord[ j3 ].keys():
        r1 = coord[ j1 ][ i1 ]
        r2 = coord[ j2 ][ i2 ]
        r3 = coord[ j3 ][ i3 ]

        dotprod = ( r1[0] - r2[0] ) * ( r3[0] - r2[0] ) + \
            ( r1[1] - r2[1] ) * ( r3[1] - r2[1] ) + \
            ( r1[2] - r2[2] ) * ( r3[2] - r2[2] )

        dist2_A = ( r1[0] - r2[0] ) * ( r1[0] - r2[0] ) + \
            ( r1[1] - r2[1] ) * ( r1[1] - r2[1] ) + \
            ( r1[2] - r2[2] ) * ( r1[2] - r2[2] )
        dist2_B = ( r3[0] - r2[0] ) * ( r3[0] - r2[0] ) + \
            ( r3[1] - r2[1] ) * ( r3[1] - r2[1] ) + \
            ( r3[2] - r2[2] ) * ( r3[2] - r2[2] )

        angle = (180/pi) * acos( dotprod / sqrt( dist2_A * dist2_B ) )

    return angle

for i in range( len(resnums) ):
    resnum = resnums[ i ]
    dist_N_CA = get_dist( coord, resnum  , ' N  ', resnum  ,' CA ' );
    dist_CA_C = get_dist( coord, resnum  , ' CA ', resnum  ,' C  ' );
    dist_C_N  = get_dist( coord, resnum  , ' C  ', resnum+1,' N  ' );
    angle_N_CA_C = get_angle( coord, resnum, ' N  ', resnum, ' CA ', resnum, ' C  ');
    angle_CA_C_N = get_angle( coord, resnum, ' CA ', resnum, ' C  ', resnum+1, ' N  ');
    angle_C_N_CA = get_angle( coord, resnum, ' C  ', resnum+1, ' N  ', resnum+1, ' CA ');

    print  '%d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f' % ( resnum, dist_N_CA, dist_CA_C, dist_C_N, angle_N_CA_C, angle_CA_C_N, angle_C_N_CA )



