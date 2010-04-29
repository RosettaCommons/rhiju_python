#!/usr/bin/python

from sys import argv, stdout
import string
from math import sqrt
from parse_options import parse_options
args = argv

pdbfile = args[1]

anchor_resnum = parse_options( args, 'anchor_res', 1 )
fixed_res = parse_options( args, 'fixed_res', [-1] )
fade = parse_options( args, 'fade', 0 )
STDEV = parse_options( args, 'stdev', 0.5 )

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
