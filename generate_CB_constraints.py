#!/usr/bin/python

from sys import argv,stdout
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

STDEV = 4.0
if len( args ) > 2:
    STDEV = float( args[2] )

lines = open( pdbfile ).readlines()

oldresnum = ''

CB_position = {}
resnums = {}
count = 0
for line in lines:

    if (len(line)>54 and  line[0:4] == 'ATOM' ):

        resnum = line[23:26]

        atom_name = line[11:15]
        position = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

        if atom_name == '  CB':
            CB_position[ count ] = position
            resnums[ count ] = int(resnum)
            count = count + 1

        oldresnum = resnum

def get_dist( pos1, pos2 ):
    dist2 = 0.0
    for k in range(3):
        dsub = ( pos1[k] - pos2[k] )
        dist2 += dsub*dsub

    return sqrt( dist2 )

cst_file = pdbfile+'.cst'
#fid = open( cst_file, 'w' )
#print "Generating ...  " , cst_file
fid = stdout
fid.write( "[ atompairs ]\n" )

DIST_CUT = 10.0

for i in range( count ):

    for j in range( count ):

        if ( abs( i - j ) <= 3 ): continue

        if len( fixed_res ) > 0 and (i in fixed_res) and (j in fixed_res): continue

        dist = get_dist( CB_position[ i ], CB_position[ j ] )
        #print i,j,dist
        if ( dist < DIST_CUT ):
            if fade:
                fid.write( " CB %d   CB %d  FADE %8.3f %8.3f %8.3f %8.3f %8.3f\n" % \
                               ( resnums[i],resnums[j],\
                                    (dist-2*STDEV), (dist+2*STDEV), STDEV, -10.0, 10.0  ) )
            else:
                fid.write( " CB %d   CB %d  HARMONIC %8.3f %8.3f \n" % (resnums[i],resnums[j],dist,STDEV) )

fid.close()
