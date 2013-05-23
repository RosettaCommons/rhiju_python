#!/usr/bin/python

from sys import argv
from os import system,chdir
from os.path import basename, exists
from parse_options import parse_options
import string

decoy_pdb = argv[1]
new_working_directory = argv[2]
example_working_directory = argv[3]

assert( exists( decoy_pdb ) )

##########################################################################
# Set up the directory -- grab files from previous directory that worked.
##########################################################################
if not exists( new_working_directory ):
    system( 'mkdir -p '+new_working_directory )

assert( exists( example_working_directory ) )
command = 'rsync -avz %s/ %s --exclude="region*"' % ( example_working_directory, new_working_directory )
system( command )

command = 'rsync -avz %s %s --exclude="region*"' % ( decoy_pdb, new_working_directory )
system( command )


chdir( new_working_directory )

##########################################################################
# Check the "README_SETUP" script...
##########################################################################
assert( exists( "README_SETUP" ) )
lines = open( "README_SETUP" ).readlines()

endpoints = parse_options( argv, "endpoints", [-1] )

outfile = "README_SETUP"
fid = open( outfile, 'w' )
for line in lines:
    if line.count( "grinder_dagman.py" ):

        # Add start_pdb information
        line += " -start_pdb %s " % basename( decoy_pdb )

        # Add endpoint information
        line += " -endpoints"
        if len( endpoints ) > 0:
            for m in endpoints: line += ' %d' % m
        else:
            cols = string.split( line )
            assert( cols.count( "-fasta" ) )
            pos = cols.index( "-fasta" )
            fasta_file = cols[ pos+1 ]

            sequence_lines = open( fasta_file  ).readlines()[1:]
            sequence = string.join(  map( lambda x : x[:-1], sequence_lines) ,  '' )

            for m in [1,len(sequence)]: line += ' %d' % m

        # Remove any virtual res tags!!!
        cols = string.split( line )
        if cols.count( "-virtual_res") > 0:
            pos = cols.index( "-virtual_res" )
            del( cols[ pos ] )
            while pos < len( cols ) and cols[ pos ][0] != '-': del( cols[ pos ] )
            line = string.join( cols )

        line += '\n'

        # If there's a new pathway to follow, create pathway_file, and
        #  specify this in the line

    fid.write( line )
fid.close()
print "Edited: ",outfile


outfile = "README_SUB"
lines = open( outfile ).readlines()
fid = open( outfile, 'w' )
for line in lines:
    if line.count( "pseudo_dagman" ):

        cols = string.split( line )
        if cols.count( "-j") > 0:
            pos = cols.index( "-j" )
            cols[ pos+1 ] = "20"
            line = string.join( cols )
            line += '\n'

    fid.write( line )
fid.close()
print "Edited: ",outfile



print
print "Set up working directory : ", new_working_directory


