#!/usr/bin/python

import string
from sys import argv
from os import system

files = argv[1:]
system( 'mkdir -p original' )

for pdb in files:

    assert( len( pdb ) == 4 )
    command  = 'curl http://www.rcsb.org/pdb/files/%s.pdb.gz > original/%s.pdb.gz' % ( string.upper( pdb ), string.upper( pdb ) )
    print command
    system( command )

    command = 'gunzip original/%s.pdb.gz' % string.upper( pdb )
    print command
    system( command )

    print
    print 'Fetched: original/%s.pdb' % string.upper( pdb )

    command = 'get_pdb.py original/%s.pdb A' % string.upper( pdb )
    print command
    system( command )

    print
    print 'Extracted: %s.pdb' % string.lower( pdb )


