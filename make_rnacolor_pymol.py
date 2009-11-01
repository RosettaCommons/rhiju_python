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
prefix = pdbfiles[0].replace('.pdb','')

if len(pdbfiles) > 1:
    command = "python ~rhiju/python/superimpose.py "
    for pdbfile in pdbfiles: command += " "+pdbfile

#    if len(highlight_residues) > 0:
#        command += " -subset "
#        for i in highlight_residues:
#            command += " %d " % i

    command += " -R 4.0 > "+ prefix+"_superposition.pdb"
#    command += " > "+ prefix+"_superposition.pdb"

    print( command )
    system(command)

#Extract models
    command = "python ~rhiju/python/parse_NMR_models.py "+prefix+"_superposition.pdb"
    system(command)
else:
    command = "cp "+pdbfiles[0]+" "+prefix+"_superposition_001.pdb"
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
fid.write('\n')
fid.write('select a, resn rA\n')
fid.write('select c, resn rC\n')
fid.write('select g, resn rG\n')
fid.write('select u, resn rU\n')
fid.write('\n')
fid.write('select bases, name c2+c4+c5+c6+c8+n1+n2+n3+n4+n6+n7+n9+o2+o4+o6+n1p\n')
fid.write('select backbone, name o1p+o2p+o3p+p+c1*+c2*+c3*+c5*+o2*+o3*+o4*+o5*\n')
#fid.write('select backbone, name c4*')
fid.write('select sugar, name c1*+c2*+c3*+c4*+o2*+o4*\n')
fid.write('\n')

fid.write('set line_width=3.0\n')
fid.write('show lines, bases\n')

#Hey! This is a big pain in the butt!
count = 0
for pdbfile in pdbfiles:
    count += 1
    lines = popen('~rhiju/python/pdb2fasta.py '+pdbfile).readlines()
    sequence = lines[1]
    for i in range(len(sequence)):
        fid.write('select extrabond1, name p and resi %d and model%d\n' \
                      % (i+1,count) )
        if (sequence[i]=='a' or sequence[i]=='g'):
            fid.write('select extrabond2, name n9 and resi %d and model%d\n' \
                          % (i+1,count) )
        else:
            fid.write('select extrabond2, name n1 and resi %d and model%d\n' \
                          % (i+1,count) )
        fid.write('bond extrabond1,extrabond2\n')

fid.write('select extrabondlines, name p\n')
fid.write('show lines, extrabondlines\n')


fid.write('color lightblue,g\n')
fid.write('color palegreen,c\n')
fid.write('color lightorange,a\n')
fid.write('color salmon,u\n')
fid.write('\n')
fid.write('select highlight, resi ')

#for res in highlight_residues: fid.write('%d+' % res)
fid.write('1-1000')
fid.write('\n')


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
fid.write('set stick_radius=0.5\n')
fid.write('show cartoon, backbone\n')
fid.write('cartoon rect, backbone\n')
fid.write('bg_color white\n')
fid.write('\n')


######
#Color backbone
BLACK_AND_WHITE = 0
count = 0
for pdbfile in pdbfiles:
    count += 1
    fid.write('select model%d_backbone, model%d and backbone\n' % (count,count))
    if not BLACK_AND_WHITE:
        fid.write('cmd.spectrum(selection = "model%d_backbone")\n' % count)


if  BLACK_AND_WHITE:
    # New idea ... black/white backbones...
    lines = open( pdbfiles[0] ).readlines()
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

                    if (dist2 > 25.0 ):
                        chain_ends.append( resnum-1 )

    fid.write('color gray80,backbone\n')


    if len( chain_ends ) > 0:
        fid.write( 'color gray20, backbone and resi 1-%d\n' % chain_ends[0] )
    if len( chain_ends ) > 1:
        fid.write( 'color gray50, backbone and resi %d-%d\n' % (chain_ends[0]+1,chain_ends[1] ) )



fid.write('\n')
fid.close()

#Output graphics
fid = open(prefix+'2.pml','w')

count = 0
for pdbfile in pdbfiles:
    count += 1
    fid.write('disable all\n')
    fid.write('enable model%d\n' % count)
    fid.write('ray 1200,1200\n')
#    fid.write('ray 1200,1200\n')
    fid.write('save '+pdbfile+'.png\n')

fid.close()

#Output graphics
fid = open(prefix+'3.pml','w')

count = 0
fid.write('hide everything,not highlight\n')

count = 0
for pdbfile in pdbfiles:
    count += 1
    if (count==1):
        fid.write('enable model%d\n'%count)
        fid.write('color white,model%d\n'%count)
#        fid.write('hide everything, asel+gsel+csel+usel and model%d\n' % count)
#        fid.write('show lines, asel+gsel+csel+usel and model%d\n' % count)
#    elif (count== len(pdbfiles)):
#        fid.write('disable model%d\n' % count)
    else:
        fid.write('enable model%d\n'%count)
fid.write('zoom highlight\n')
fid.write('ray 1200,1200\n')
fid.write('save '+prefix+'_superposition.png \n')
