#!/usr/bin/python

from sys import argv,stdout
from os import system,popen
from math import sqrt
import string

inputfiles = argv[1:]


pdbfiles = []
highlight_residues = []
for i in range( len(inputfiles) ):
    inputfile = inputfiles[i]
    try:
        highlight_residues.append( int( inputfile) )
    except:
        if (not inputfile.find("superposition") > 0):
            pdbfiles.append( inputfile)

#pdbfiles.reverse()
prefix = "TEST"

#superimpose first!
command = "python ~rhiju/python/superimpose.py "
for pdbfile in pdbfiles: command += " "+pdbfile
command += " -R 2.0  > superposition.pdb"
system(command)

#Extract models
command = "python ~rhiju/python/parse_NMR_models.py superposition.pdb"
system(command)

fid = open(prefix+'.pml','w')
#fid = stdout


fid.write('reinitialize\n')
count = 0
for pdbfile in pdbfiles:
    count += 1
    fid.write('load superposition_%03d.pdb,model%d\n' %
              (count, count))


fid.write('\n')
fid.write('hide everything,all\n')

#fid.write('show lines, all\n')
fid.write('color grey, (element c)\n')
fid.write('color red, (element o)\n')
fid.write('color blue, (element n)\n')
fid.write('color yellow, (element s)\n')

fid.write('select backbone, name o+c+ca+n\n')
fid.write('select c_backbone, name c+ca\n')

#fid.write( 'create heavyatoms, not element h\n')
#fid.write('show sticks, not element h\n')
#fid.write('set stick_radius, 0.15, heavyatoms\n')

#fid.write( 'create backbone_obj, backbone\n')
#fid.write('set stick_radius, 0.4, backbone_obj\n')

fid.write('bg_color white\n')
fid.write('show cartoon, all\n')
fid.write('set cartoon_oval_length, 0.5\n')
fid.write('set cartoon_oval_width, 0.5\n')
fid.write('set cartoon_rect_length, 0.5\n')


for count in range( len(pdbfiles) ):
    pdbfile = pdbfiles[ count ]
    fid.write('select model%d_backbone, model%d and c_backbone\n' % (count,count))
    fid.write('cmd.spectrum(selection = "model%d_backbone")\n' % count)

if len( pdbfiles ) == 3:
    colors = ['blue','green','red']
elif len( pdbfiles ) == 4:
    colors = ['blue','green','orange','red']
else:
    colors = ['blue','cyan','green','orange','red']

for count in range( len(pdbfiles) ):
    pdbfile = pdbfiles[ count ]
    fid.write( 'color %s, model%d\n' % ( colors[count], count+1 ) )

########################################
# What's buried?
pdbfile = pdbfiles[ 0 ]
lines = open( pdbfile ).readlines()
oldresnum = ''
CB_position = {}
count = 0
pdb_resnum = {}
for line in lines:

    if (len(line)>54 and  line[0:4] == 'ATOM' ):

        resnum = line[23:26]
        if not resnum == oldresnum:
            count = count + 1
            pdb_resnum[ count] = resnum

        atom_name = line[11:15]
        position = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

        if atom_name == '  CB':
            CB_position[ count ] = position

        oldresnum = resnum



def get_dist( pos1, pos2 ):
    dist2 = 0.0
    for k in range(3):
        dsub = ( pos1[k] - pos2[k] )
        dist2 += dsub*dsub

    return sqrt( dist2 )

nres = count
dist_cutoff = 10.0
buried = []
for  i in CB_position.keys():
    n_neighbors = 0
    for j in CB_position.keys():
        dist = get_dist( CB_position[i], CB_position[j] )
        if ( dist < dist_cutoff ): n_neighbors += 1
    if (n_neighbors > 15): buried.append( string.split(pdb_resnum[i])[0] )

########################################
fid.write( 'select buried, resi %s\n' % string.join(buried,'+') )
fid.write( 'select hydrophobic, resn leu+val+ile+ala+cys+met+phe+trp+tyr+pro\n' )
fid.write( 'show sticks, buried and hydrophobic and not elem h and not c+n+o\n' )

########################################

fid.write('select none\n');

fid.write('\n')
fid.close()


#Output graphics
fid = open(prefix+'2.pml','w')

count = 0
for pdbfile in pdbfiles:
    count += 1
    fid.write('disable all\n')
    fid.write('enable model%d\n' % count)
    fid.write('enable model%d_backbone_obj\n' % count)
    fid.write('enable HBA%d\n' % count)
    fid.write('enable HBD%d\n' % count)
    fid.write('ray 800,800\n')
    fid.write('save '+pdbfile+'.png\n')
