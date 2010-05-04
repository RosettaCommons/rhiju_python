#!/usr/bin/python

from sys import argv
from generate_constraints import generate_constraints

args = argv
backbone_atoms = [ ' CA ',' C  ',' N  ',' O  ', ' CB ' ]
generate_constraints( args, backbone_atoms, [], 0.0 )
exit( 0 )
