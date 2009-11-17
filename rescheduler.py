#!/usr/bin/python

from time import sleep
from os import system,popen
import string

condor_running = 1
while condor_running:

    lines = popen( 'condor_q rhiju' ).readlines()
    condor_running = 0
    for line in lines:
        cols =  string.split( line )
        if cols.count('rhiju'):
            condor_running += 1

    print 'Jobs running:', condor_running

    if condor_running:
        command = 'condor_reschedule'
        print( command )
        system( command  )
        sleep( 1 )
