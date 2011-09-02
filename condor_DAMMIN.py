#!/usr/bin/python

from os import system
from time import sleep
from sys import argv
NSTRUCT = 10

P2_SYMM = 0
if len( argv) == 3:
    assert( argv[2] == '-p2' )
    del( argv[2] )
    P2_SYMM = 1

assert( len(argv) == 2 )
gnom_out = argv[1]
tag = gnom_out.split( '.' )[0]

for n in range( NSTRUCT ):
    fid = open( 'condorDAMMIN', 'w' )

    fid.write( '+TGProject = TG-MCB090153\n' )
    fid.write( 'universe = vanilla\n' )
    fid.write( 'notification = never\n')

    fid.write( 'executable = /home/rhiju/src/ATSAS-2.3.2-1/bin/dammin\n' )

    p2_symm_tag = ''
    if P2_SYMM: p2_symm_tag = ' /sy P2'

    fid.write( 'arguments  =  %s  /lo %s%02d /mo Slow %s \n' % (gnom_out, tag, n, p2_symm_tag ) )
    fid.write( 'Queue 1\n' )

    fid.close()

    if ( n > 0 ): sleep( 20 )

    system( 'condor_submit condorDAMMIN' )

