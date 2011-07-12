#!/usr/bin/python

import string
from sys import argv

infile = argv[1]
lines = open( infile ).readlines()

pairs = []
count = 0
structure = ''
for line in lines[1:]:
    cols =string.split( line )
    if cols[1]=='ENERGY': break
    pos1 = int( cols[5] )
    pos2 = int( cols[4] )
    if pos2 > 0:
        if pos2 > pos1:
            structure += '('
        else:
            structure += ')'
    else:
        structure +=  '.'

print structure

