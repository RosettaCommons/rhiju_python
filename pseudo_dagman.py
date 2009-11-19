#!/usr/bin/python

from sys import argv,exit
import string
from os import system,popen
from os.path import basename,exists
from random import randrange
from time import sleep

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


##############################################
# This doesn't seem to be very robust...
# Some issues with file server synchronization with bjobs queue listing
def check_for_job_tags( job_tags_ ):
    lines = popen( "bjobs -w" ).readlines()
    still_running = 0
    for job_tag in job_tags_:
        for line in lines:
            cols = string.split( line )
            if ( cols.count( job_tag ) ):
                still_running += 1
                break
    return still_running


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
def biox2_condor_submit( condor_submit_file_ ):
    print 'Run: ', condor_submit_file_
    lines = open( condor_submit_file_ ).readlines()
    log = "/dev/null"
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

    if ( exe == "" ) or  ( args == "" ) or (output == "") or (queue_num == 0):
        print "PROBLEM!!!"
        exit()

    job_tags_ = []
    output_files_ = []
    for q in range( queue_num ):
        job_tag_ = basename( condor_submit_file_.replace('.condor','') ).upper()
        job_tag_ += "_%d_%d" % ( q, randrange(0,1000000) )

        q_string = '%d' % q
        args_new = args.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)
        err_new = err.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)
        output_new = output.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)

        command = "bsub -W 4:0 -o %s -e %s -J %s    %s %s " % \
                  (output_new,err_new,job_tag_,   exe,args_new)

        wait_for_pending_jobs_to_clear()

        print command
        system( command )

        job_tags_.append( job_tag_ )
        output_files_.append( output_new )

        #Need to wait until jobs actually appear in the queue...
        #check_for_job_tags( job_tags_ )

    #print "Sleeping to let queue update... "
    #sleep( 5 )
    #return log_files_
    return output_files_

done = {}
queued = {}
for job in jobs:
    done[ job ] = 0
    queued[ job ] = 0

#job_tags = {}
output_files = {}
all_done = 0

while not all_done:

    ###################################################
    # Find jobs that are ready to go but not done.
    ###################################################
    all_done = 1
    for job in jobs:
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
                        command = pre_script[ job ]
                        print( command )
                        system( command )

                    output_files[ job ] = biox2_condor_submit( condor_submit_file[ job ] )
                    queued[ job ] = 1

    if not all_done:
        sleep(1)

    ###################################################
    # Check for any jobs that are done
    ###################################################

    for job in jobs:
        if not done[ job ] and queued[ job ]:
            still_running = check_output_files( output_files[ job ] )
            print "Jobs still running: ",output_files[job][-1],still_running

            if not still_running:
                if job in post_script.keys():
                    command = post_script[ job ]
                    print( command )
                    system( command )
                done[ job ] = 1
                queued[ job ] = 0


