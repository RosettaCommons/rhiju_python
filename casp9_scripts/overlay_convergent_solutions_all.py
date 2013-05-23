#!/usr/bin/python

import string
from glob import glob
from sys import argv
from os import system,popen,chdir,getcwd
from os.path import exists,basename
from parse_options import parse_options
from math import sqrt
from make_rhiju_color import make_rhiju_color

globdirs = glob( 'loop*-*/' )
globdirs.sort()
#print globdirs

if len( globdirs ) == 0:
    print 'Looking for directories of the form: loop*-*'
    exit( 0 )

override = parse_options( argv, "override", 0 )

cwd = getcwd()
all_pdbs = []
for dir in globdirs:
    print '------ Entering: ' + dir + '--------'
    chdir( dir )
    command = 'overlay_convergent_solutions.py'
    if override: command += ' -override'
    system( command )
    chdir( cwd )

    pdbs = glob( dir+'pdb/loop*pdb' )
    pdbs.sort()

    for pdb in pdbs: all_pdbs.append( pdb )



if not exists( 'pdb' ):
    system( 'mkdir -p pdb' )

for pdb in all_pdbs:
    system( 'cp %s pdb' % pdb )

#chdir( 'pdb' )
#files = map( lambda x:basename(x), all_pdbs )
#make_rhiju_color( files )


print '##########################'
print ' useful commands'
print '##########################'
for dir in globdirs:
    print 'rs %s/pdb/TEST.script & ' % dir



