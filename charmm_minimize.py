#!/usr/bin/python

from sys import argv,exit
from os.path import exists,basename,dirname,abspath,expanduser
from os import system, popen, chdir, getcwd
import string
from random import random
from os import putenv

pdb_file = abspath( argv[1] )

WORKDIR = dirname( pdb_file )+'/'
HOMEDIR = expanduser( '~' )

CWD = getcwd()
chdir( WORKDIR )

putenv( 'MMTSBDIR',HOMEDIR+'/updated_mmtsb/' )
putenv( 'CHARMMEXEC',HOMEDIR+'/mmtsb/bin/charmm30a1' )
putenv( 'CHARMMDATA',HOMEDIR+'/mmtsb/data/charmm' )


MIN_EXE = HOMEDIR+'/updated_mmtsb/perl/minCHARMM.pl'
if not (exists( basename(pdb_file)+'.min_pdb' ) ):
    command = 'perl '+MIN_EXE+' -par gb=gbmva,nuc5ter=none,nuc3ter=none,deoxy=0 %s > %s.min_pdb' % ( basename(pdb_file), basename(pdb_file) )
    print( command )
    system( command )

    SCORE_EXE = HOMEDIR+'/updated_mmtsb/perl/enerCHARMM.pl '
    command = 'perl '+SCORE_EXE + ' -par gb=gbmva,nuc5ter=none,nuc3ter=none,deoxy=0 -all %s.min_pdb > %s.scores' % ( basename(pdb_file), basename(pdb_file) )
    print( command )
    system( command )

chdir( CWD )

