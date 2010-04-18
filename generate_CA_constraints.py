#!/usr/bin/python

from sys import argv
import string
from math import sqrt

args = argv

pdbfile = args[1]

fade = 0
if args.count( '-fade' ):
    pos = args.index( '-fade' )
    del( args[ pos ] )
    fade = 1


fixed_res = []
if args.count( '-fixed_res' ):
    pos = args.index( '-fixed_res' )
    del( args[ pos ] )
    goodint = 1
    while goodint:
        try:
            fixed_residue = int(args[pos])
            fixed_res.append( fixed_residue )
            del( args[ pos ] )
        except:
            goodint = 0

STDEV = 0.5
if len( args ) > 2:
    STDEV = float( args[2] )

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
fid = open( cst_file, 'w' )
print "Generating ...  " , cst_file

fid.write( "[ atompairs ]\n" )

DIST_CUT = 12.0

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
