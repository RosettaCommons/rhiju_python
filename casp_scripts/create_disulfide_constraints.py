#!/usr/bin/python

from sys import argv, stdout
from os.path import exists
import string

fasta = argv[ 1 ]
assert( exists( fasta ) )
sequence_lines = open( fasta  ).readlines()[1:]
sequence = string.join(  map( lambda x : x[:-1], sequence_lines) ,  '' )

cys_res = []
for i in range( len ( sequence ) ):
    if ( sequence[i] == 'C' ): cys_res.append( i+1 )
#print cys_res

stdout.write( '[ atompairs ]\n' )
for i in cys_res:
    for j in cys_res:
        if ( i >= j ): continue
        stdout.write( 'SG %3d SG %3d   FADE  -8.0 8.0  4.0  -10.0 0.0\n'  % (i,j) )


