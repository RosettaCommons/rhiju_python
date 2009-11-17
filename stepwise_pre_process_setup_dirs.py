#!/usr/bin/python

from sys import argv
from os import system
from os.path import exists,basename
 
outdir = argv[1]
dir_prev = argv[2]
condor_submit_file = argv[3]

# Need directories to be setup. Could do this with a preprocessing script? That would also allow
# for convenient setup of Queue number in the condor script.
MAX_JOBS = 4000
for q in range( MAX_JOBS ):
    pdbfile = '%s/%s_sample.cluster.out.%d.pdb' % (dir_prev,dir_prev.lower(),q)
    if exists( pdbfile ):
        tag = basename( pdbfile ).replace( '.pdb' ,'' )
        newdir = outdir+'/START_FROM_'+tag.upper()
        if not exists( newdir ):  system( 'mkdir -p '+newdir )
    else:
        break

N_JOBS = q

# Go through condor submission file and update number of jobs...
# May need to be careful about any special cases.
lines = open( condor_submit_file ).readlines()
fid = open( condor_submit_file, 'w' )

for line in lines:
    if len( line ) > 5 and line[:5] == 'Queue':
        fid.write( 'Queue %d\n' % N_JOBS )
    else:
        fid.write( line )

fid.close()


