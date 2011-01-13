#!/usr/bin/python

from sys import argv,stdout
import string
from math import sqrt
from parse_options import parse_options
from generate_constraints import generate_constraints

generate_constraints( argv, [ " N1 ", " N2 "," N3 ", " N4 " ], \
                          [ " N3 ", " O2 ", " O4 ", " O6 ", " N7 ", " N1 ", " N6 " ],\
                          3.2 )
