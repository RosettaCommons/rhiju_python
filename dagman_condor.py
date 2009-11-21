#!/usr/bin/python


from sys import argv,exit
import string
from os import system,popen
from os.path import basename,exists
from time import sleep

#######################################
# This could actually be a class!
#######################################

##############################################
def check_output_files( output_files_ ):

    still_running = 0
    done_running = 0

    for output_file in output_files_:

        if not exists( output_file ):
            # Need to mark that we are not done!
            still_running += 999
            continue

        lines = open( output_file ).readlines()

        jobs_running = []
        for line in lines:
            if line.find( 'Job executing' ) > 0:
                job_tag = string.split( line, '(')[1].split(')')[0]
                jobs_running.append( job_tag )
                print 'EXECUTED', job_tag

        if len( jobs_running) == 0: # logfile not filled out yet.
            still_running += 999
            continue

        terminated = {}
        for job in jobs_running: terminated[ job ] = 0

        for line in lines:
            if line.find( 'Job terminated' ) > 0:
                job_tag = string.split( line, '(')[1].split(')')[0]
                terminated[ job_tag ] = 1

        for job in jobs_running:
            still_running += 1
            if  terminated[ job ]: done_running += 1

    still_running = still_running - done_running
    print 'Checkaroo: ',still_running, done_running

    command = 'condor_reschedule > /dev/null 2> /dev/null '
    #print( command )
    system( command )

    return still_running

##############################################
def wait_for_pending_jobs_to_clear():
    too_many_pending_jobs = 1
    PEND_LIMIT = 199
    while too_many_pending_jobs:
        too_many_pending_jobs = 0
        lines = popen( "condor_q " ).readlines()

        pos = 0 # Look for "STATUS" column
        for count in range( len( lines ) ):
            if len( line ) > 3 and line[:3] == " ID":
                pos = line.find('ST')
                break
        num_idle = 0

        for count in range( len(lines ) ):
            if len( line ) > pos:
                if ( line[pos:pos+1] == 'I' ):
                    num_idle += 1

        if num_idle >= PEND_LIMIT:
            too_many_pending_jobs = 1
            print "Waiting for pending job number ", len(lines), " to fall below ", PEND_LIMIT
            sleep( 10 )

##############################################
def condor_submit( condor_submit_file_ ):
    print 'Run: ', condor_submit_file_
    lines = open( condor_submit_file_ ).readlines()
    log = ""
    output = ""
    err = "/dev/null"
    exe = ""
    args = ""
    queue_num = 0
    for line in lines:
        if len( line ) > 2:
            cols = string.split( line )
            if cols[0] ==  "executable":
                assert( cols[1] == "=" )
                exe = cols[2]
            elif cols[0] == "arguments":
                assert( cols[1] == "=" )
                args = string.join(cols[2:])
            elif cols[0] == "log":
                assert( cols[1] == "=" )
                log = cols[2]
            elif cols[0] == "output":
                assert( cols[1] == "=" )
                output = cols[2]
            elif cols[0] == "error":
                assert( cols[1] == "=" )
                err = cols[2]
            elif (cols[0]).lower() == "queue":
                if len( cols ) > 1: queue_num = int( cols[1] )

    if ( exe == "" ) or  ( args == "" ) or (output == "") or (queue_num == 0) or ( log == "" ):
        print "PROBLEM!!!"
        exit()

    job_tags_ = []
    output_files_ = []
    actually_queued = 0

    if exists( log ) and (not check_output_files( [ log ] )):
        print "Already done", log
    else:
        command = "condor_submit "+condor_submit_file_
        print( command )
        system( command )
        actually_queued = 1


    output_files_.append( log  ) # NOTE: on condor systems, check the ".log" files for queuing and finishing.

    return ( output_files_, actually_queued )
