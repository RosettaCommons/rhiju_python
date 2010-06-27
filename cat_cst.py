#!/usr/bin/python

from sys import argv,stdout

cst_files = argv[1:]

coord_lines = []
atompair_lines = []

for file in cst_files:
    lines = open( file ).readlines()
    in_coord_cst = 0
    in_atompair_cst = 0
    for i in range( len( lines ) ):
        line = lines[ i ]
        if ( len( line ) >= 15 and line[:15] == "[ coordinates ]" ):
            in_coord_cst = 1
            in_atompair_cst = 0
            continue
        if ( len( line ) >= 13 and line[:13] == "[ atompairs ]" ):
            in_coord_cst = 0
            in_atompair_cst = 1
            continue
        if (in_coord_cst):
            coord_lines.append( line )
        elif (in_atompair_cst):
            atompair_lines.append( line )
        else:
            assert( in_coord_cst or in_atompair_cst )

if len( coord_lines ) > 0:
    stdout.write( "[ coordinates ]\n" )
    for line in coord_lines: stdout.write( line )

if len( atompair_lines ) > 0:
    stdout.write( "[ atompairs ]\n" )
    for line in atompair_lines: stdout.write( line )
