#!/usr/bin/python

from sys import argv
import string

outfiles = argv[1:]

def output_file( outfile, count, header_lines, current_lines ):

    outfile_new = outfile.replace( '.out', '_%05d.out' % count )

    fid = open( outfile_new, 'w' )
    for l in header_lines: fid.write( l )
    for l in current_lines: fid.write( l )
    fid.close()


for outfile in outfiles:

    lines = open( outfile ).readlines()

    NUM_SCORE_LINES = 0
    current_lines = []
    header_lines = []
    count = 0

    for line in lines:

        if ( string.split( line )[0] == 'SCORE:' ):
            NUM_SCORE_LINES += 1

            if ( NUM_SCORE_LINES == 2 ):
                header_lines = current_lines
                current_lines = []

            if NUM_SCORE_LINES > 2:
                output_file( outfile, count, header_lines, current_lines )
                current_lines = []
                count += 1

        if ( string.split( line )[0] == 'SEQUENCE:' and NUM_SCORE_LINES > 0):
            NUM_SCORE_LINES = 0
            current_lines = []

        current_lines.append( line )

    output_file( outfile, count, header_lines, current_lines )
