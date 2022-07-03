#!/usr/bin/python

from sys import argv
import string
from string import upper

sequence = string.join(argv[1:],'')

complement = { 'A':'T', 'T':'A', 'U':'A', 'C':'G', 'G':'C', 'N':'N' }
make_DNA = { 'A':'A', 'T':'T', 'U':'T', 'C':'C', 'G':'G', 'N':'N' }

seq1 = ''
seq2 = ''
numchar = len( sequence )

print numchar
scanchar = numchar

MAXCHAR = 6000
if (numchar > MAXCHAR ):
    print 'Truncating down to', MAXCHAR,' !!!!!'
    scanchar = MAXCHAR

for i in range( scanchar ):
    seq1 += make_DNA[ sequence[i].upper() ]
    seq2 += complement[ sequence[numchar - 1 - i].upper() ]


print seq1
print seq2

