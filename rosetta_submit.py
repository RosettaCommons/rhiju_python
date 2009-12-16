#!/usr/bin/python

from sys import argv,exit
from os import system
from os.path import basename,dirname
import string

if len( argv ) < 4:
    print argv[0]+' <text file with rosetta command> <outdir> <# jobs>'
    exit()

infile = argv[1]
outdir = argv[2]
try:
    n_jobs = int( argv[3] )
except:
    print 'NEED TO SUPPLY NUMBER OF JOBS'


lines = open(infile).readlines()

bsub_file = 'bsubMINI'
fid = open( bsub_file,'w')

tot_jobs = 0

for line in  lines:

    if len(line) == 0: continue
    if line[0] == '#': continue
    if string.split( line[0]) == []: continue

    for i in range( n_jobs ):
        dir = outdir+'/%03d/' % (tot_jobs)
        system( 'mkdir -p '+ dirname(dir) )
        command_line = line[:-1].replace( '-out:file:silent ', '-out:file:silent '+dir)
        command_line = command_line.replace( '-out::file::silent ', '-out::file::silent '+dir)
        command_line = command_line.replace( 'macosgcc', 'linuxgcc')
        command_line = command_line.replace( 'Users', 'home')

        cols = string.split( command_line )
        if '-total_jobs' in cols:
            pos = cols.index( '-total_jobs' )
            cols[ pos+1 ] = '%d' % n_jobs
            command_line = string.join( cols )
        if '-job_number' in cols:
            pos = cols.index( '-job_number' )
            cols[ pos+1 ] = '%d' % i
            command_line = string.join( cols )


        command =  'bsub -W 16:0 '+command_line
        fid.write( command + '\n')

        tot_jobs += 1

fid.close()

print 'Created bsub submission file ',bsub_file,' with ',tot_jobs, ' jobs queued. To run, type: '
print '>source',bsub_file
