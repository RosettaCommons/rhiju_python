#!/usr/bin/python

from sys import argv
from os import system, getcwd, chdir
from os.path import exists,abspath,expanduser
from glob import glob

outdirs = argv[1:]

try:
    NUM_JOBS_PER_NODE = int( outdirs[-1] )
    del( outdirs[ -1 ] )
except:
    NUM_JOBS_PER_NODE = 20

HOMEDIR = expanduser('~')
EXE = HOMEDIR+'/python/get_pbCHARMM.py'

bsub_file = open( 'bsubCHARMM','w' )
condor_file = open( 'CHARMM.condor','w' )
condor_file.write('+TGProject = TG-MCB090153\n')
condor_file.write('universe = vanilla\n')
condor_file.write('executable = %s\n' % EXE )

CWD = getcwd()

total_count = 0

for outdir in outdirs:

    print outdir

    assert( exists( outdir ) )

    count = 0
    start = 0
    globfiles = glob( outdir+'/S_*OUT/S*.min_pdb' ) # Must be in CHARMM format, already minimized
    globfiles.sort()
    print len( globfiles ),

    for file in globfiles:

        if ( (count % NUM_JOBS_PER_NODE) == 0):
            if ( start == 1 ):
                bsub_file.write( '\n' )
                condor_file.write( '\nQueue 1\n' )
            else:
                start = 1
            bsub_file.write( '\nbsub -W 16:0 %s ' % EXE )
            condor_file.write( '\narguments = ' )
        count += 1
        bsub_file.write( ' '+file )
        condor_file.write( ' '+file )

    if (not count % NUM_JOBS_PER_NODE == 0):
        bsub_file.write( '\n')
        condor_file.write( '\nQueue 1\n')

    chdir( CWD )

    total_count += count

print 'Total number of PDBs to minimize: ', total_count

