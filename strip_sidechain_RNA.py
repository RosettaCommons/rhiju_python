#!/usr/bin/python

from sys import argv,stderr
from os import system
import string


def Help():
    print argv[0]," <pdb1> <pdb2> ... "
    print '  Remove sidechains to help in design'
    print '  Note: all pdbs must have the same sequence.'


pdbfiles = argv[1:]
if len(pdbfiles) < 1:
    Help()

backbone_atoms = [' P  ',' O1P',' O2P',' O5*',' C5*',' C4*',' O4*',' C3*',' O3*',' C2*',' O2*',' C1*','2HO*']
backbone_atoms = backbone_atoms + [" OP1"," OP2"," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," O2'"," C1'","2HO'"]

for pdbfile in pdbfiles:

    new_pdbfile = pdbfile[:-4]+'_stripsidechain.pdb'

    fid = open( new_pdbfile,'w')

    lines = open(pdbfile).readlines()
    for line in lines:
        writeout = 0
        if line[:4]=='ATOM':
            old_res = line[17:20]

            #new_res = map_res[ old_res ]
            new_res = old_res

            atom = line[12:16]

            if atom in backbone_atoms:
                writeout = 1
                new_atom = atom


        if writeout:
            fid.write(line[:12]+new_atom+' '+new_res+line[20:])

    fid.close()
    print "Created: ", new_pdbfile
