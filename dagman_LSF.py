#!/usr/bin/python


from sys import argv,exit
import string
from os import system,popen
from os.path import basename,exists
from random import randrange
from time import sleep


##############################################
def check_output_files( output_files_ ):

    still_running = 0
    done_running = 0

    for output_file in output_files_:
        still_running += 1
        #print "CHECKING ", output_file

        if not exists( output_file ): continue
        lines = popen( "tail -n 20 "+output_file ).readlines()

        for line in lines:
            if ( len( line ) > 22 and  line[:23] == "Successfully completed." ):
                #print "completed!"
                done_running += 1
                break
    #print "checkaroo: ", still_running, done_running
    still_running = still_running - done_running
    return still_running

##############################################
def wait_for_pending_jobs_to_clear():
    too_many_pending_jobs = 1
    PEND_LIMIT = 199
    while too_many_pending_jobs:
        too_many_pending_jobs = 0
        lines = popen( "bjobs -w | grep PEND " ).readlines()
        if len( lines ) >= PEND_LIMIT:
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

    output_files_ = []
    commands = []
    for q in range( queue_num ):

        q_string = '%d' % q
        args_new = args.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)
        err_new = err.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)
        output_new = output.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)
        output_files_.append( output_new )

        if exists( output_new ) and not check_output_files( [ output_new ] ):
            print "Already done", output_new
            already_done += 1
        else:
            command = "bsub -W 4:0 -o %s -e %s  %s %s " % \
                      (output_new,err_new, exe,args_new)

            commands.append( command )

    actually_queued = 0
    for command in commands: # There may no be anything to queue...
        wait_for_pending_jobs_to_clear()
        print command
        system( command )
        actually_queued = 1

    return ( output_files_, actually_queued )
