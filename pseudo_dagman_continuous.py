#!/usr/bin/python

from sys import argv,exit,stdout
import string
from os import system,popen
from os.path import basename,exists
from random import randrange
from time import sleep
from glob import glob

from dagman_LSF_continuous import *
#from dagman_condor_continuous import *


# In this implementation, kick off a certain number
#  of long-running scripts -- they will carry out communication
#  this master script based on files that show up in
#  job-specific directories.
N_JOBS = 100
if argv.count( '-j' ):
    pos = argv.index( '-j' )
    del( argv[pos] )
    N_JOBS = int( argv[ pos] )
    del( argv[pos] )

finalize = 1
if argv.count( '-no_finalize' ):
    pos = argv.index( '-no_finalize' )
    del( argv[pos] )
    finalize = 0


#########################################
# Parse dagman file
#########################################
dagman_file = argv[1]
lines = open( dagman_file ).readlines()

jobs = []
condor_submit_file = {}
post_script = {}
pre_script = {}
parents = {}
for line in lines:
    print line[:-1]
    stdout.flush()
    if len( line ) > 4 and line[:3] == "JOB":
        cols = string.split( line )
        job = cols[1]
        jobs.append( job )
        condor_submit_file[ job ] = cols[2]
    elif len( line ) > 6 and line[:6] == "SCRIPT":
        cols = string.split( line )
        job = cols[2]
        if cols[1] == "PRE":
            pre_script[ job ] = string.join( cols[3:] )
        else:
            assert( cols[1] == "POST" )
            post_script[ job ] = string.join( cols[3:] )
    elif len( line ) > 7 and line[:6] == "PARENT":
        cols = string.split( line )
        assert( cols.count( "CHILD" ) )
        child_index = cols.index( "CHILD" )
        parent_jobs =  cols[1:child_index]
        child_jobs =  cols[child_index:]
        for child_job in child_jobs:
            if child_job not in parents.keys():
                parents[ child_job ] = []
            for parent_job in parent_jobs: parents[ child_job ].append( parent_job )

# Kick off jobs!
job_cluster_number = kick_off_slave_jobs( N_JOBS )

done = {}
queued = {}
for job in jobs:
    done[ job ] = 0
    queued[ job ] = 0

#job_tags = {}
output_files = {}
all_done = 0
early_exit = 0

while not all_done:

    ###################################################
    # Find jobs that are ready to go but not done.
    ###################################################
    all_done = 1
    for job in jobs:
        queued_a_job = 0
        if not done[ job ]:
            all_done = 0

            if not queued[ job ]:
                #Consider queuing the job

                ready_to_queue = 1
                if job in parents.keys():
                    for parent_job in parents[job]:
                        if not done[ parent_job ]:
                            ready_to_queue = 0
                            break

                if ready_to_queue:
                    if job in pre_script.keys():
                        pre_command = pre_script[ job ]
                        pre_command_log_file = condor_submit_file[job] +'.pre_script.log'
                        if not exists( pre_command_log_file  ):
                            command =  pre_command + ' > '+pre_command_log_file
                            system( command )

                    ( output_files[ job ], actually_queued ) = condor_submit( condor_submit_file[ job ] )

                    if len( output_files[ job ] ) == 0 and not actually_queued:
                        print "problem?"
                        #all_done = 1
                        #early_exit = 1

                    queued[ job ] = actually_queued
                    if actually_queued:
                        queued_a_job = 1
                    else:
                        done[ job ] = 1

            #if ( queued_a_job ): break # some other jobs may have finished while we were queuing.

    if not all_done:
        sleep(1)

    if len( glob( 'core.*' ) ) > 0 : early_exit = 1 # totally uncool to have cores.

    if early_exit: break

    ###################################################
    # Check for any jobs that are done
    ###################################################

    for job in jobs:
        if not done[ job ] and queued[ job ]:
            still_running = check_output_files( output_files[ job ] )
            print "Jobs still running: ",output_files[job][-1],still_running
            stdout.flush()

            if not still_running:
                if job in post_script.keys():
                    command = post_script[ job ]
                    print( command )
                    system( command )
                done[ job ] = 1
                queued[ job ] = 0


if finalize: finish_jobs()

