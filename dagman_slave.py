#!/usr/bin/python

from os import system
from os.path import exists
from sys import argv
from time import sleep

jobdir = argv[ 1 ]

assert( exists( jobdir ) )

finished_file = jobdir + '/finished.txt'

while not exists( finished_file ):

    command_file_name = jobdir + '/run_this_script.txt'

    ran_a_job = 0
    if exists( command_file_name ):

        sleep( 2 ) # stupid file locking...

        lines = open( command_file_name ).readlines()
        command = 'source '+lines[0]
        print command
        system( command  )

        # After done, remove the script file...
        # A more robust signal might be to create a "done" file.
        system( 'rm -rf ' + command_file_name )

        ran_a_job = 1

    sleep( 2 )



