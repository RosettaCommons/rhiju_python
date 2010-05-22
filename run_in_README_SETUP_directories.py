#!/usr/bin/python

from sys import argv
from os import popen,system,chdir,getcwd
from os.path import dirname
import string


command = string.join( argv[1:] )

lines = popen( "find ./ -name 'README_SETUP'" ).readlines()

cwd = getcwd()
for line in lines:

    chdir( dirname( line ) )
    system( command )
    chdir( cwd )
