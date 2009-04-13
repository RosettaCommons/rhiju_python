#!/usr/bin/python

from sys import argv,exit
from os.path import exists,basename,dirname,abspath
from os import system, popen, chdir, getcwd
import string
from random import random
from os import putenv

pdb_file = abspath( argv[1] )

WORKDIR = dirname( pdb_file )+'/'

CWD = getcwd()
chdir( WORKDIR )

putenv( 'MMTSBDIR','/home/rhiju/updated_mmtsb/' )
putenv( 'CHARMMEXEC','/home/rhiju/mmtsb/bin/charmm30a1' )
putenv( 'CHARMMDATA','/home/rhiju/mmtsb/data/charmm' )



SCORE_EXE = '/home/rhiju/updated_mmtsb/perl/enerCHARMM.pl '
command = 'perl '+SCORE_EXE + ' -par gb=gbmva,nuc5ter=none,nuc3ter=none,deoxy=0 -all %s' % ( basename(pdb_file) )
print( command )
system( command )

chdir( CWD )

