#!/usr/bin/python

from sys import argv,stderr
from os import system
import string


def Help():
    print argv[0]," <pdb1> <pdb2> ... "
    print '  Remove sidechains to help in design'
    print '  Note: all pdbs must have the same sequence.'


pdbfiles = argv[1:-1]
if len(pdbfiles) < 1:
    Help()

sequence = argv[-1]

resname = { 'a':'RAD','c':'RCY','g':'RGU','u':'URA'}


for pdbfile in pdbfiles:

    new_pdbfile = pdbfile[:-4]+'_newseq.pdb'

    fid = open( new_pdbfile,'w')

    lines = open(pdbfile).readlines()
    count = -1
    old_resnum = '   '
    for line in lines:

        writeout = 0
        if line[:4]=='ATOM':
            resnum = line[22:26]
            if not old_resnum == resnum:
                count = count + 1
            old_resnum = resnum

            old_res = line[17:20]
            new_res = resname[ sequence[count] ]
            writeout = 1



        if writeout:
            fid.write(line[:16]+' '+new_res+line[20:])

    fid.close()

