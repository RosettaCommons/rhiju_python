#!/usr/bin/python

from sys import argv
from os import popen,system
import string

SHOW_CH_BOND = 0
SHOW_H_BOND = 1

infiles = argv[1:]

#superimpose first!
#prefix = 'TEST'
prefix = infiles[0].replace('.pdb','')

if len(infiles) > 1:
    command = "python ~rhiju/python/superimpose.py "
    for pdbfile in infiles: command += " "+pdbfile

    command += " -R 2.0 > "+ prefix+"_superposition.pdb"
    print( command )
    system(command)

#Extract models
    command = "python ~rhiju/python/parse_NMR_models.py "+prefix+"_superposition.pdb"
    system(command)
else:
    command = "cp "+infiles[0]+" "+prefix+"_superposition_001.pdb"
    system(command)

pml_file = prefix + ".pml"
fid = open( pml_file, 'w' )

fid.write('reinitialize\n')

for i in range( len( infiles) ):
    infile = prefix+"_superposition_%03d.pdb" % (i+1)
    fid.write('load %s,decoy%d\n' % (infile,i+1) )

fid.write('\n')
fid.write('hide everything,all\n')
fid.write('\n')
fid.write('select a, resn rA\n')
fid.write('select c, resn rC\n')
fid.write('select g, resn rG\n')
fid.write('select u, resn rU\n')
fid.write('\n')
fid.write('select bases, name c2+c4+c5+c6+c8+n1+n2+n3+n4+n6+n7+n9+o2+o4+o6+n1p\n')
fid.write('select backbone, name o1p+o2p+o3p+p+c1*+c2*+c3*+c5*+c4*+o2*+o3*+o4*+o5*\n')
fid.write('select sugar, name c1*+c2*+c3*+c4*+o2*+o4*\n')
fid.write('select o2star, name o2*')
fid.write('\n')
fid.write('\n')
fid.write('show spheres, o2star\n');
fid.write('alter o2star, vdw=0.25\n')
fid.write('set line_width=2.0\n')
fid.write('show lines, bases\n')
fid.write('\n')
fid.write('select highlight, resi 1-999\n')
fid.write('select asel, resn rA and highlight and bases\n')
fid.write('select csel, resn rC and highlight and bases\n')
fid.write('select gsel, resn rG and highlight and bases\n')
fid.write('select usel, resn rU and highlight and bases\n')
fid.write('color blue,gsel\n')
fid.write('color green,csel\n')
fid.write('color orange,asel\n')
fid.write('color red,usel\n')
fid.write('show sticks, asel+csel+gsel+usel\n')
fid.write('\n')
fid.write('set stick_radius=0.25\n')
fid.write('show sticks, backbone\n')
fid.write('bg_color white\n')
fid.write('\n')
#fid.write('select decoy_backbone, decoy and backbone\n')
#fid.write('cmd.spectrum(selection = "decoy_backbone")\n')
fid.write('\n')
fid.write('show lines, hydro\n')
fid.write('color gray50, hydro\n')
fid.write('\n')


########################################
for i in range( len(infiles) ):
    infile = prefix+"_superposition_%03d.pdb" % (i+1)

    command = 'rna_test.macosgccrelease  -database ~rhiju/minirosetta_database -print_hbonds -s '+infile

    lines = popen( command ).readlines()

    for line in lines:
        cols = string.split( line )
        if len( cols ) > 0 and cols[0]=='HBOND:':
            res1 = cols[1][1:]
            atom1 = cols[2]
            res2 = cols[4][1:]
            atom2 = cols[5]
            energy = float( cols[7] )
            if (energy > -0.7 ): continue
            fid.write('dist decoy%d_HB, (decoy%d and name %s  and resi  %s  ), (decoy%d and name  %s and resi  %s  )\n' % (i+1,i+1,atom1,res1,i+1,atom2,res2) )
            fid.write('\n')

    fid.write('hide labels, decoy%d_HB\n' % (i+1))
    fid.write('color red, decoy%d_HB\n' % (i+1))
    if (not SHOW_H_BOND):
        fid.write('disable decoy%d_HB\n' % (i+1) )

#    for line in lines:
#        cols = string.split( line )
#        if len( cols ) > 0 and cols[0]=='CHbond:':
#            res1 = cols[2]
#            atom1 = cols[6]
#            res2 = cols[12]
#            atom2 = cols[14]
#            if float( cols[17] ) > -0.75: continue
#            fid.write('dist decoy%d_CHB, (decoy%d and name %s  and resi  %s  ), (decoy%d and name  %s and resi  %s  )\n' % (i+1,i+1, atom1,res1,i+1,atom2,res2) )
#    fid.write('\n')
#    fid.write('hide labels, decoy%d_CHB\n' % (i+1))
#    fid.write('color cyan, decoy%d_CHB\n' % (i+1))
#    if (not SHOW_CH_BOND):
#        fid.write('disable decoy%d_CHB\n' % (i+1) )

##########################
# Figure out chain definitions
lines = open( infiles[0] ).readlines()
chain_ends = []

coords_p = {}
coords_o3star = {}

for line in lines:
    if line[:4] == 'ATOM':
        # need to be very smart about chainbreaks...
        resnum = int( line[22:26] )
        atomname = line[12:16]
        x = float( line[30:38] )
        y = float( line[38:46] )
        z = float( line[46:54] )
        if (atomname == ' O3*'): coords_o3star[resnum] = [x,y,z]
        if (atomname == ' P  '):
            coords_p[resnum] = [x,y,z]
            if (resnum-1) in coords_o3star.keys():
                c1 = coords_o3star[resnum-1]
                c2 = coords_p[ resnum ]
                dist2 = \
                    (c1[0]-c2[0])*(c1[0]-c2[0]) + \
                    (c1[1]-c2[1])*(c1[1]-c2[1]) + \
                    (c1[2]-c2[2])*(c1[2]-c2[2])

                if (dist2 > 9.0 ):
                    chain_ends.append( resnum-1 )

fid.write('color gray80,backbone\n')


if len( chain_ends ) > 0:
    fid.write( 'color gray20, backbone and resi 1-%d\n' % chain_ends[0] )
if len( chain_ends ) > 1:
    fid.write( 'color gray50, backbone and resi %d-%d\n' % (chain_ends[0]+1,chain_ends[1] ) )

fid.write('\n')

fid.close()

#################################################
# Automated png script.
fid = open( 'OUT.pml','w') # output graphics

WIDTH  = 1000
HEIGHT = 800

fid.write('disable all\n')
fid.write('enable decoy1\nenable decoy1_HB\nzoom decoy1\nray %d,%d\nsave %s_native_allatom.png\n' \
              % (WIDTH, HEIGHT, prefix) )
if len( infiles ) > 1:
    fid.write('disable all\n')
    fid.write('enable decoy2\nenable decoy2_HB\nzoom decoy2\nray %d,%d\nsave %s_decoy_allatom.png\n' \
                  % ( WIDTH, HEIGHT, prefix) )
if len( infiles ) > 2:
    for i in range( len(infiles)-2 ):
        fid.write('disable all\\n')
        fid.write('enable decoy%denable decoy%d_HB\n\nray %d,%d\nsave %s_decoy%d_allatom.png\n' % (i+3,i+3,WIDTH,HEIGHT,prefix,i+3) )

fid.close()
