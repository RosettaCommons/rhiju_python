#!/usr/bin/python
from os import system
from sys import argv
from os.path import exists

for file in argv[1:]:

    assert( file[-4:] == '.out' )

    scorefile = file.replace('.out','.sc')
    if not exists( scorefile ):
        command = 'grep SCORE '+file+' > '+scorefile
        print( command )
        system( command )

