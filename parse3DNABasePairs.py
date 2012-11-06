#!/usr/bin/python
from sys import argv,stdout
import string
from os.path import exists,basename
from os import system

def Help():
    print
    print argv[0]+ ' <input pdb>'
    print
    print 'R. Das, 2012'

if len (argv )< 2:
    Help()
    exit()

file = argv[1]

rnanum = { 'a':1, 'c':2, 'g':3, 'u':4 }

fid = open( file )

wc_bps = []
all_res = []
wc_bp_steps = []
all_local_base_pair_params = []
all_local_base_pair_step_params = []
all_torsions1 = []
all_torsions2 = []
seq1 = []
seq2 = []

line = fid.readline()
while line:
    if len( line ) > 0:
        cols = line.split()
        if len(cols)>4 and cols[3] == 'explicit': break
    line = fid.readline()


line = fid.readline()

while len( line ) > 4 and line[0] != '#' :


    res1   = int( line[30:34].replace( '.','') )
    chain1 = line[28:29]

    res2   = int( line[54:58].replace( '.','') )
    chain2 = line[60:61]

    print '%s:%d-%s:%d' % (chain1,res1,chain2,res2)

    line = fid.readline()

fid.close()
