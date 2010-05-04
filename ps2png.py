#!/usr/bin/python

from sys import argv
from os import system

for file in argv[1:]:

    file_new = file.replace('.eps','.png').replace('.ps','.png')

    command =  'convert -rotate 90 '+file+' '+file_new
    print( command )
    system( command )
