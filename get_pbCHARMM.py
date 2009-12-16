#!/usr/bin/python

from sys import argv
from os import putenv,popen
from os.path import expanduser,exists

putenv( 'MMTSBDIR','/home/rhiju/updated_mmtsb/' )
putenv( 'CHARMMEXEC','/home/rhiju/mmtsb/bin/charmm30a1' )
putenv( 'CHARMMDATA','/home/rhiju/mmtsb/data/charmm' )

HOMEDIR = expanduser('~' )

pdbfiles = argv[1:]

ionconcs = ['0.0', '0.1', '1.0', '100.0']

for file in pdbfiles:

    # Must be in CHARMM-ready format!!!
    assert( file.count("min_pdb") > 0 )

    PB_filename = file+'.PB.txt'

    if exists( PB_filename ):
        print "Already calculated", PB_filename
        continue

    energy = {}
    for ionconc in ionconcs:
        command = 'perl /home/rhiju/updated_mmtsb/perl/pbCHARMM.pl  -par pbionconc=%s,nuc5ter=none,nuc3ter=none,deoxy=0 %s ' % (ionconc, file )
        print( command )
        lines = popen( command ).readlines()
        energy[ ionconc ] = float( lines[0][:-1] )

    pbfile = open( PB_filename, 'w' )

    for ionconc in ionconcs:
        line = 'PB_%s %8.3f' % (ionconc,energy[ionconc] )
        print line
        pbfile.write( line+'\n' )

    pbfile.close()
