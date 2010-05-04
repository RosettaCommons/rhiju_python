#!/usr/bin/python

from sys import argv, stdout
from generate_constraints import generate_constraints

args = argv
generate_constraints( args, ' CA ', [], 0.0 )
exit( 0 )
