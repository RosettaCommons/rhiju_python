#!/usr/bin/python

from sys import argv,stdout
from os import popen

fastafile = argv[1]
if fastafile[-4:] == '.pdb':
    fasta = popen( '~rhiju/python/pdb2fasta.py '+fastafile).readlines()[1][:-1]
else:
    assert( fastafile[-6:] == '.fasta' )
    fasta = popen( ' cat '+fastafile ).readlines()[1][:-1]


endpoints = map( lambda x:int(x), argv[2:] )

bound_start = -1
bound_end = -1

assert( 1 in endpoints )
assert( len( fasta ) in endpoints )

print
def print_sequence( fid, fasta, bound_start, bound_end ):
    for i in range( len( fasta ) ):
        if ( i+1 < bound_start or i+1 > bound_end):
            fid.write(' ')
        else:
            fid.write(fasta[ i ])
    fid.write('\n')

for i in range( len( endpoints ) - 1 ):
    point_start = endpoints[ i ]
    point_end = endpoints[ i+1 ]

    if ( bound_start == -1 ):
        bound_start = point_start
        bound_end = point_start

    if ( point_end - point_start  ) > 0:

        bound_end_current = bound_end
        for i in range( point_end - bound_end_current ):
            bound_end += 1
            print_sequence( stdout, fasta, bound_start, bound_end )

    else:

        bound_start_current = bound_start
        for i in range( bound_start_current - point_end ):
            bound_start -= 1
            print_sequence( stdout, fasta, bound_start, bound_end )

assert( bound_start == 1 )
assert( bound_end == len( fasta ) )
