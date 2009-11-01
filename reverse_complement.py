#!/usr/bin/python

from sys import argv
import string

sequence = string.join(argv[1:],'')

complement = { 'A':'T', 'T':'A', 'U':'A', 'C':'G', 'G':'C' }
make_DNA = { 'A':'A', 'T':'T', 'U':'T', 'C':'C', 'G':'G' }

seq1 = ''
seq2 = ''
numchar = len( sequence )

print numchar
scanchar = numchar

MAXCHAR = 60
if (numchar > MAXCHAR ):
    print 'Truncating down to', MAXCHAR,' !!!!!'
    scanchar = MAXCHAR

for i in range( scanchar ):
    seq1 += make_DNA[ sequence[i] ]
    seq2 += complement[ sequence[numchar - 1 - i] ]


print seq1
print seq2

