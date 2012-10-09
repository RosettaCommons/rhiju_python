#!/usr/bin/python

import string
from sys import argv
from os import system
from os.path import exists

files = argv[1:]

RNA = False
if files.count( '-rna' ) > 0:
    pos = files.index('-rna')
    del( files[ pos ] )
    RNA = True

if files.count( '-chain' ):
    pos = files.index('-chain')
    del( files[ pos ] )
    chain = argv[pos]
    del( files[ pos ] )

chain = ''

system( 'mkdir -p original' )

for pdb in files:

    assert( len( pdb ) == 4 )

    if not exists( 'original/%s.pdb' % string.upper( pdb ) ):
        command  = 'curl http://www.rcsb.org/pdb/files/%s.pdb.gz > original/%s.pdb.gz' % ( string.upper( pdb ), string.upper( pdb ) )
        print command
        system( command )

        command = 'gunzip original/%s.pdb.gz' % string.upper( pdb )
        print command
        system( command )

    print
    print 'Fetched: original/%s.pdb' % string.upper( pdb )

    if RNA:

        command = 'make_rna_rosetta_ready.py original/%s.pdb %s' % ( string.upper( pdb ), chain )
        print command
        system( command )

    else:
        command = 'get_pdb.py original/%s.pdb %s' % ( string.upper( pdb ), chain )
        print command
        system( command )

        print
        print 'Extracted: %s.pdb' % string.lower( pdb )


