#!/usr/bin/python

from os import system,chdir
from time import sleep
from sys import argv
from glob import glob

NSTRUCT = 10


if len( argv ) > 1:
    pdbfiles = argv[1:]
else:
    pdbfiles = glob( '*-1.pdb' )
    pdbfiles.sort()
    print pdbfiles

system( 'rm -rf damaver/' )
system( 'mkdir damaver' )

for i in range( len( pdbfiles ) ):
    file = pdbfiles[ i ]
    new_file = 'temp_file%02d.pdb' % i
    command = 'cp %s damaver/%s' % (file,new_file)
    print command
    system( command )

chdir( 'damaver' )
system( 'damaver /a' )
chdir( '..' )
system( 'cp damaver/dam* .' )
