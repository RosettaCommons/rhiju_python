#!/usr/bin/python

import sys
import string
from math import sqrt

files = sys.argv[3:]
atomtype = sys.argv[1:3]

atom1_defined = 0
atom2_defined = 0

for file in files:

    lines = open( file, 'r' ).readlines()
    for line in lines:
        cols = string.split(line)
        if len(cols)<3:  continue
        if cols[2] == atomtype[0]:
            N_atom = (float(cols[5]), float(cols[6]), float(cols[7]))
            atom1_defined = 1
        if cols[2] == atomtype[1]:
            CA_atom = (float(cols[5]), float(cols[6]), float(cols[7]))
            atom2_defined = 1
        if atom1_defined and atom2_defined:
            length = sqrt(
                (N_atom[0] - CA_atom[0])*(N_atom[0] - CA_atom[0]) +
                (N_atom[1] - CA_atom[1])*(N_atom[1] - CA_atom[1]) +
                (N_atom[2] - CA_atom[2])*(N_atom[2] - CA_atom[2]))
            print length
            atom1_defined = 0
            atom2_defined = 0
