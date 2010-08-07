#!/usr/bin/python

from sys import argv,stdout

secstruct = argv[1];

stdout.write( "[ atompairs ]\n" )

for n in range( 5, len(secstruct)+1 ):
    # in helix?
    in_helix = True
    for i in range( n-4, n+1 ):
        if ( not secstruct[ i-1 ] == 'H' ):
             in_helix = False
             break
    if not in_helix: continue

    stdout.write('  N %4d    O  %4d   FADE 2.2 3.8 0.5 -10.0 10.0\n' % (i,i-4) )
