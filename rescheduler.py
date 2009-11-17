#!/usr/bin/python

from time import sleep
from os import system,popen,getlogin
import string

username = getlogin()

condor_running = 1
while condor_running:

    lines = popen( 'condor_q '+username ).readlines()
    condor_running = 0
    for line in lines:
        cols =  string.split( line )
        if cols.count( username ):
            condor_running += 1

    print 'Jobs running:', condor_running

    if condor_running:
        command = 'condor_reschedule'
        print( command )
        system( command  )
        sleep( 1 )
