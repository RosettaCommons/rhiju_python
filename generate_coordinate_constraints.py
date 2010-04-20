#!/usr/bin/python

from sys import argv, stdout
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

backbone_atoms = [ ' CA ',' C  ',' N  ',' O  ', ' CB ' ]

backbone_atom_position = []
backbone_atom_name = []
backbone_resnum = []

count = 0
for line in lines:

    if (len(line)>54 and  line[0:4] == 'ATOM' ):

        resnum = line[23:26]
        if not resnum == oldresnum:
            count = count + 1

        atom_name = line[12:16]
        position = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

        if atom_name in backbone_atoms:
            backbone_atom_position.append( position )
            backbone_atom_name.append( atom_name )
            backbone_resnum.append( count )

        oldresnum = resnum

#cst_file = pdbfile+'.cst'
#fid = open( cst_file, 'w' )
fid = stdout
#print "Generating ...  " , cst_file

fid.write( "[ coordinates ]\n" )

anchor_atom_name = " CA "
anchor_resnum    = 1 # Anchor atom.

dist = 0.0
for i in range( len( backbone_atom_name ) ):

    if len( fixed_res ) > 0 and (backbone_resnum[ i ] in fixed_res) : continue
    if fade:
        fid.write( " %s %d   %s %d    %8.3f %8.3f %8.3f    FADE %8.3f %8.3f %8.3f %8.3f %8.3f\n" % \
                       ( backbone_atom_name[ i ],
                         backbone_resnum[ i ],
                         anchor_atom_name,
                         anchor_resnum,
                         backbone_atom_position[ i ][0],
                         backbone_atom_position[ i ][1],
                         backbone_atom_position[ i ][2],
                         (dist-2*STDEV), (dist+2*STDEV), STDEV, -10.0, 10.0  ) )
    else:
        fid.write( " %s %d   %s %d    %8.3f %8.3f %8.3f    HARMONIC %8.3f %8.3f \n" % \
                       ( backbone_atom_name[ i ],
                         backbone_resnum[ i ],
                         anchor_atom_name,
                         anchor_resnum,
                         backbone_atom_position[ i ][0],
                         backbone_atom_position[ i ][1],
                         backbone_atom_position[ i ][2],
                         dist,STDEV) )

fid.close()
