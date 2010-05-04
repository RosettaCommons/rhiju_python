#!/usr/bin/python

from sys import argv,stdout
from os import popen

endpoints = map( lambda x:int(x), argv[1:] )

bound_start = -1
bound_end = -1

print 'PATH',

#assert( 1 in endpoints )
#assert( len( fasta ) in endpoints )
#print

#def print_sequence( fid, fasta, bound_start, bound_end ):
#    for i in range( len( fasta ) ):
#        if ( i+1 < bound_start or i+1 > bound_end):
#            fid.write(' ')
#        else:
#            fid.write(fasta[ i ])
#    fid.write('\n')

for i in range( len( endpoints ) - 1 ):
    point_start = endpoints[ i ]
    point_end = endpoints[ i+1 ]

    if ( bound_start == -1 ):
        bound_start = point_start
        bound_end = point_start
        print point_start,

    if ( point_end - point_start  ) > 0:

        bound_end_current = bound_end
        for i in range( point_end - bound_end_current ):
            bound_end += 1
            print bound_end,

    else:

        bound_start_current = bound_start
        for i in range( bound_start_current - point_end ):
            bound_start -= 1
            print bound_start,


print

#assert( bound_start == 1 )
#assert( bound_end == len( fasta ) )

