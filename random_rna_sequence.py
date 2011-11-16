#!/usr/bin/python

from sys import argv
import random

try:
    N = int( argv[1] )
except:
    N = 30

add_tail = 0
if argv.count( '-tail'):
    add_tail = 1

rna_letters = 'ACGU'
seq = ''
for i in range( N ):
    seq = seq + rna_letters[ int( 4 * random.random() ) ]

if add_tail:
    seq = 'GGCCAAAACAAC'+seq+'GAAACAAAGAAACAACAACAACAAC'
print seq
