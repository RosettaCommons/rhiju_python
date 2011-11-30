#!/usr/bin/python

from sys import argv
from os.path import exists
from os import system, chdir, getcwd
from math import sqrt
from glob import glob

job_list = argv[1]
lines = open( job_list ).readlines()

def get_dist( pos1, pos2 ):
    dist2 = 0.0
    for k in range(3):
        dsub = ( pos1[k] - pos2[k] )
        dist2 += dsub*dsub

    return sqrt( dist2 )

def sub( pos1, pos2 ):
    return [ pos1[0]-pos2[0],  pos1[1]-pos2[1],  pos1[2]-pos2[2] ]
def dot( pos1, pos2 ):
    return pos1[0]*pos2[0] + pos1[1]*pos2[1] + pos1[2]*pos2[2]

CWD = getcwd()

print ' pdb [str-stp]:    dist    cos1  cos2  needs N^2'

for line in lines:

    pdb = line[:-1]

    chdir( pdb )
    loop_file = '%s.loop' % pdb
    cols = open( loop_file ).readlines()[0].split()
    loop_start = int(cols[0])
    loop_stop  = int(cols[1])

    pdb_files = glob( 'noloop*pdb')
    if len( pdb_files ) == 0: pdb_files = glob( '*excise*pdb' )
    pdb_file = pdb_files[0]

    pdblines = open( pdb_file ).readlines()
    xyz_N = []
    xyz_C = []
    xyz_CA = []
    for pdbline in pdblines:
        if pdbline[:4] in ['ATOM','HETA']:
            if (pdbline[12:16]==' N  '):
                xyz_N.append( [float(pdbline[30:38]), float(pdbline[38:46]), float(pdbline[46:54])] )
            if (pdbline[12:16]==' C  '):
                xyz_C.append( [float(pdbline[30:38]), float(pdbline[38:46]), float(pdbline[46:54])] )
            if (pdbline[12:16]==' CA '):
                xyz_CA.append( [float(pdbline[30:38]), float(pdbline[38:46]), float(pdbline[46:54])] )

    #d = get_dist( xyz_C[ loop_start-1-1 ], xyz_N[ loop_stop+1-1] )
    dist = get_dist( xyz_CA[ loop_start-1-1 ], xyz_CA[ loop_start-1] )

    try:
        v1 = sub( xyz_CA[ loop_start-1-1], xyz_CA[ loop_start-2-1] )
        v2 = sub( xyz_CA[ loop_start-1], xyz_CA[ loop_start+1-1] )
        r =  sub( xyz_CA[ loop_start-1-1 ], xyz_CA[loop_start-1] )
        #r =  sub( xyz_CA[ loop_start-2-1 ], xyz_CA[loop_start+1-1] )

        c1 = dot( v1, r ) / sqrt( dot( r,r ) * dot( v1,v1) )
        c2 = -dot( v2, r ) / sqrt( dot( r,r ) * dot( v2,v2) )
    except:
        c1 = 0.0
        c2 = 0.0

    asterisk = ''
    costheta_cutoff = 0.4
    costheta_cutoff2 = 0.1
    if ( dist < 16) and \
           ( (c1 > costheta_cutoff) or (c2>costheta_cutoff) or  (c1>costheta_cutoff2 and c2>costheta_cutoff2) ):
        asterisk = '*'
    print '%s [%3d-%3d]: %8.2f  %5.2f %5.2f  %s' % (pdb, loop_start-1, loop_stop+1, dist, c1, c2, asterisk)

    chdir( CWD )
