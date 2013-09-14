#!/usr/bin/python
from sys import argv
from os import system
from parse_options import parse_options
from os.path import exists, basename, dirname, expanduser
import string


executable_dir = parse_options( argv, 'exe_dir', '' )
database = parse_options( argv, 'd', '' )

submit_script = argv[1]
lines = open( submit_script ).readlines()

executable = ''
arguments = ''
for line in lines:
    if line[:10] == 'executable' and len( executable ) == 0:
        executable = string.split( line, '=' )[1]
        executable = executable.replace( ' ','')[:-1]
    if line[:9] == 'arguments':
        arguments = string.split( line, '=' )[1]
        arguments = arguments.replace( '$(Process)', '0' )

assert( len( executable ) > 0)
assert( len( arguments ) > 0 )

cols = arguments.split()
pos = cols.index( '-out:file:silent' )
outfile = cols[pos+1]
if not exists( dirname( outfile ) ):
    print 'Making directory: ', dirname( outfile )
    system( 'mkdir -p '+dirname(outfile) )

pos = cols.index( '-database' )
database_original = cols[ pos + 1 ]
if len( database ) == 0: database = database_original

if len( executable_dir ) > 0:
    executable = executable_dir + '/'+ basename( executable )

if not exists( executable ):
    executable = executable.replace( 'linux','macos' ).replace( '/home','/Users')
if not exists( executable ):
    executable = executable.replace( '~',expanduser( '~' ) )
if not exists( database ):
    database = database.replace( '/home', '/Users' )
if not exists( database ):
    database = database.replace( '~',expanduser( '~' ) )

cols[ pos+1 ] = database
arguments = string.join( cols )

print executable
assert( exists( executable ) )
#print database
assert( exists( database ) )

command = executable +' '+arguments
print command
system( command )
