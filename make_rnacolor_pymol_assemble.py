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

pdbfiles.reverse()

#superimpose first!
command = "python ~rhiju/python/superimpose.py "
for pdbfile in pdbfiles: command += " "+pdbfile
prefix = pdbfiles[0].replace( '.pdb','')
command += " -R 20.0 >  "+ prefix+"_superposition.pdb"
system(command)

#Extract models
command = "python ~rhiju/python/parse_NMR_models.py "+prefix+"_superposition.pdb"
system(command)

fid = open(prefix+'.pml','w')
#fid = stdout


fid.write('reinitialize\n')
count = 0
for pdbfile in pdbfiles:
    count += 1
    fid.write('load %s_superposition_%03d.pdb,model%d\n' %
              (prefix,count, count))


fid.write('\n')
fid.write('hide everything,all\n')
fid.write('show cartoon,all\n')
fid.write('cartoon oval\n')
fid.write('select a, resn rA\n')
fid.write('select c, resn rC\n')
fid.write('select g, resn rG\n')
fid.write('select u, resn rU\n')
fid.write('\n')
fid.write('select bases, name c2+c4+c5+c6+c8+n1+n2+n3+n4+n6+n7+n9+o2+o4+o6+n1p\n')
fid.write('select backbone, name o1p+o2p+o3p+p+c1*+c2*+c3*+c4*+c5*+o2*+o3*+o4*+o5*\n')
fid.write('select sugar, name c1*+c2*+c3*+c4*+o2*+o4*\n')
fid.write('\n')
fid.write('color lightblue,g\n')
fid.write('color palegreen,c\n')
fid.write('color lightorange,a\n')
fid.write('color tv_red,u\n')

fid.write('color blue,g\n')
fid.write('color green,c\n')
fid.write('color orange,a\n')
fid.write('color red,u\n')

fid.write('show sticks, bases\n')
#fid.write( 'show cartoon, backbone\n')
fid.write( 'set cartoon_ring_mode, 1\n')
fid.write( 'set cartoon_oval_width, 2.0\n')
fid.write( 'set cartoon_oval_length, 2.0\n')
fid.write('bg_color white\n')
fid.write('\n')


CONNECT_BONDS = 0

if CONNECT_BONDS:
    count = 0;
    for pdbfile in pdbfiles:
        count += 1
        lines = popen('~rhiju/python/pdb2fasta.py '+pdbfile).readlines()
        sequence = lines[1]
        for i in range(len(sequence)-1):
            fid.write('bond  (name o3* and resi %d and model%d), (name p and resi %d and model%d)\n' \
                          % (i+1, count,i+2,count) )


count = 0
for pdbfile in pdbfiles:
    count += 1
    fid.write('select model%d_backbone, model%d and backbone\n' % (count,count))
    fid.write('cmd.spectrum(selection = "model%d_backbone")\n' % count)

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
    fid.write('ray 800,800\n')
    fid.write('save '+pdbfile+'.png\n')
