#!/usr/bin/python

from sys import argv,stdout
from os import system,popen


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

#superimpose first!
prefix = 'TEST'#pdbfiles[0].replace( '.pdb','')

if len( pdbfiles ) > 1:
    command = "python ~rhiju/python/superimpose.py "
    for pdbfile in pdbfiles: command += " "+pdbfile
    command += " >  "+ prefix+"_superposition.pdb"
    system(command)

    #Extract models
    command = "python ~rhiju/python/parse_NMR_models.py "+prefix+"_superposition.pdb"
    system(command)

fid = open(prefix+'.pml','w')
#fid = stdout


fid.write('reinitialize\n')
count = 0
if ( len(pdbfiles) > 1 ):
    for pdbfile in pdbfiles:
        count += 1
        fid.write('load %s_superposition_%03d.pdb,model%d\n' %
                  (prefix,count, count))
else:
    fid.write('load %s,model%d\n' %
              ( pdbfiles[0], 1))


fid.write('\n')
fid.write('hide everything,all\n')

fid.write('show lines, all\n')
fid.write('color grey, (element c)\n')
fid.write('color red, (element o)\n')
fid.write('color blue, (element n)\n')
fid.write('color yellow, (element s)\n')

fid.write('select backbone, name o+c+ca+n\n')
fid.write('select c_backbone, name c+ca\n')

#fid.write( 'create heavyatoms, not element h\n')
fid.write('show sticks, not element h\n')
fid.write('set stick_radius, 0.15, heavyatoms\n')

#fid.write( 'create backbone_obj, backbone\n')
#fid.write('set stick_radius, 0.4, backbone_obj\n')

fid.write('bg_color white\n')

count = 0
for pdbfile in pdbfiles:
    count += 1
    fid.write('select model%d_backbone, model%d and c_backbone\n' % (count,count))
    fid.write('cmd.spectrum(selection = "model%d_backbone")\n' % count)


    fid.write('select model%d_don, model%d and (elem n,o and (neighbor hydro))\n' % (count,count) )
    fid.write('select model%d_acc, model%d and (elem o or (elem n and not (neighbor hydro)))\n' % (count,count) )
    fid.write('dist HBA%d, (model%d_acc),(model%d_don), 3.2\n' % (count,count,count) )
    fid.write('dist HBD%d, (model%d_don),(model%d_acc), 3.2\n' % (count,count,count) )
    fid.write('delete model%d_don\n' % count)
    fid.write('delete model%d_acc\n' % count)
    #fid.write('hide (hydro)\n')
    fid.write(' \n')
    fid.write('hide labels,HBA%d\n' % count )
    fid.write('hide labels,HBD%d\n' % count )

    fid.write('create model%d_backbone_obj, model%d and backbone\n' % (count,count))
    fid.write('set stick_radius, 0.6, model%d_backbone_obj\n' % (count) )

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
