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
fid.write('hide everything, all\nshow cartoon,all\n')
fid.write('set cartoon_round_helices, 1\n')
fid.write('set cartoon_discrete_colors,0\n')
fid.write('set antialias, 1\n')
fid.write('set backface_cull, 0\n')
fid.write('set cartoon_fancy_sheets, 1\n')
fid.write('set cartoon_fancy_helices, 1\n')
fid.write('set cartoon_dumbbell_length, 1.0\n')
fid.write('set cartoon_rect_length, 0.5\n')
fid.write('bg_color white\n')

colors = ['blue','red','green','black']
for i in range( len( infiles) ):
    infile = prefix+"_superposition_%03d.pdb" % (i+1)
    fid.write('color %s,decoy%d\n' % (colors[i],i+1) )
