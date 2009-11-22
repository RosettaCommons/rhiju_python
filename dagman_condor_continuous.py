#!/usr/bin/python


from sys import argv,exit
import string
from glob import glob
from os import system,popen,getcwd
from os.path import basename,exists,expanduser,dirname
from time import sleep

#######################################
# This could actually be a class!
#######################################

HOMEDIR = expanduser('~')

PYDIR = HOMEDIR+'/python'
assert( exists( PYDIR ) )

SLAVE_EXE = PYDIR + '/dagman_slave.py'

def kick_off_slave_jobs( N_JOBS ):

    for n in range( N_JOBS ):
        JOBDIR = 'SLAVE_JOBS/%d' % n
        if not exists( JOBDIR) :
            command = 'mkdir -p '+JOBDIR
            print( command )
            system( command )

    condor_submit_file = 'SLAVE_JOBS/slave_jobs.condor'
    fid = open( condor_submit_file, 'w' )
    fid.write('universe = vanilla\n')
    fid.write('executable = %s\n' % SLAVE_EXE )
    fid.write('arguments = SLAVE_JOBS/$(Process)\n' )
    fid.write('output = SLAVE_JOBS/$(Process)/slave_jobs.out\n' )
    fid.write('log = SLAVE_JOBS/$(Process)/slave_jobs.log\n' )
    fid.write('error = SLAVE_JOBS/$(Process)/slave_jobs.err\n' )
    fid.write('notification = never\n')
    fid.write('Queue %d\n' % N_JOBS )
    fid.close()

    job_log_file = 'slave_jobs.log'
    command = ( 'condor_submit '+condor_submit_file+' > ' + job_log_file )
    print( command )
    system( command )

    # This could also keep track of start time...
    for n in range( N_JOBS ):
        JOBDIR = 'SLAVE_JOBS/%d' % n
        running_file = JOBDIR+'/running.txt'
        command = 'echo "running" > '+running_file
        print( command )
        system( command )

    lines =open( job_log_file ).readlines()

    cols = string.split( lines[-1] )
    assert( cols[2] == 'submitted' )
    job_cluster_number = int( cols[-1][:-2] )

    return job_cluster_number

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


def find_a_slave( qsub_script ):

    #jobdirs = glob( 'SLAVE_JOBS/*/running.txt' )
    #jobdirs.sort()
    #jobdirs = map( lambda x:dirname( x ), jobdirs )

    job_log_file = 'slave_jobs.log'
    lines =open( job_log_file ).readlines()
    cols = string.split( lines[-1] )
    assert( cols[2] == 'submitted' )
    job_cluster_number =  cols[-1][:-1]


    # Could also look at condor_queue or whatever to make sure this job is still running...
    num_slave = -1
    while num_slave < 0:

        jobdirs = []

        # Running jobs...
        lines = popen( 'condor_q rhiju | grep " R " ' ).readlines()
        for line in lines:
            cols = string.split( line, '.')
            if len(line) > 10 and cols[0] == job_cluster_number:
                jobdir = 'SLAVE_JOBS/' + string.split( cols[1], ' ')[0]
                jobdirs.append( jobdir )

        for jobdir in jobdirs:
            command_file_name = jobdir + '/run_this_script.txt'
            if not exists( command_file_name ):
                fid = open( command_file_name, 'w' )
                fid.write( '%s\n' % qsub_script )
                fid.close()
                num_slave = int( basename( jobdir ) )
                break

        if num_slave >= 0:
            print 'assigning  %s to slave %d' % ( qsub_script, num_slave )
        else:
            # Did not find anything ... wait a couple seconds
            print "waiting for a slave to free up."
            sleep( 2 )

    return num_slave

def finish_jobs():
    jobdirs = glob( 'SLAVE_JOBS/*/slave_jobs.log' )
    jobdirs.sort()
    jobdirs = map( lambda x:dirname( x ), jobdirs )
    for jobdir in jobdirs:
        fid = open( jobdir + '/finished.txt' )
        fid.write( ' Finished!! ' )
        fid.close()

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

    # Stolen from script for PBS.
    output_files_ = []
    qsub_scripts = []
    for q in range( queue_num ):

        q_string = '%d' % q
        args_new = args.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)
        err_new = err.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)
        output_new = output.replace('$(Process)',q_string).replace('$(PROCESS)',q_string)
        output_files_.append( output_new )

        already_done = 0
        if exists( output_new ) and not check_output_files( [ output_new ] ):
            print "Already done", output_new
            already_done += 1
        else:
            # Need to prepare a script for qsub
            qsub_script = err_new.replace('.err','.qsub')
            fid = open( qsub_script, 'w' )
            fid.write( '#!/bin/bash\n' )
            fid.write( 'cd %s\n\n' % getcwd() )
            command = "%s %s > %s 2> %s\n\n" % \
                      (exe, args_new,output_new,err_new)
            fid.write( command )
            command = "echo 'Successfully completed.\n' >> %s " % \
                      ( output_new )
            fid.write( command )

            qsub_scripts.append( qsub_script )

    actually_queued = 0
    if not already_done:
        for qsub_script in qsub_scripts: # There may not be anything to queue...
            num_slave = find_a_slave( qsub_script )
            actually_queued = 1

    return ( output_files_, actually_queued )


