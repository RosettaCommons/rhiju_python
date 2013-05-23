#!/usr/bin/python

import string
from sys import argv
from os import getcwd
from os.path import abspath
from parse_options import parse_options
from get_sequence import get_sequence

author_code = '7006-8004-3833'

if len( argv ) < 3:
    print argv[0]," <list of 5 pdbs> <CASP target number> "
    print "   Additional possible options: "
    print "     -stdev <float>         [error expected over residues]"
    print "     -stdev_loose <float>   [looser error expected over some residues]"
    print "     -loose_res <residues>  [the loose residues]"
    print "     -no_cst_res <residues> [crazy residues --> error = 12 A]"

list_file = argv[1]
target_num = int( argv[2] )

current_path = abspath( getcwd() )
if ( string.find( current_path, argv[2] ) < 0 ):
    print
    print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    print " Could not find target number ", argv[2], " in current path!", current_path
    print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    print

################################################
# What are the files?
################################################
files = map( lambda x:x[:-1], open( list_file ).readlines() )
assert( len( files ) > 0 )
assert( len( files ) <= 5 )

sequence = get_sequence( files[0] )
totres = len( sequence )

################################################
# To specify errors... same format as used to
#  generate constraints.
#
# -stdev (general error, by default 2 A)
#
# -stdev_loose (loose error, by default 4 A)
# -loose_res (residues that are given stdev_loose errors)
#
# -no_cst_res (residues with very poor errors, assumed 12 A)
#
# OR:  "per_res" file.
################################################
per_res_file = parse_options( argv, "per_res", "" )
per_res_scale = parse_options( argv, "per_res_scale", 1.5 )

if len( per_res_file ):
    lines = open( per_res_file ).readlines()
    expected_errors = []
    for line in lines:
        expected_errors.append( per_res_scale * float( string.split( line[:-1] )[1] ) )
else:
    STDEV = parse_options( argv, "stdev", 2.0 )
    STDEV_LOOSE = parse_options( argv, "stdev_loose", 4.0 )
    loose_res = parse_options( argv, "loose_res", [-1] )
    cst_res = parse_options( argv, "cst_res", [-1] )
    no_cst_res = parse_options( argv, "no_cst_res", [-1] )
    BIG_ERROR = 8.0

    if len( cst_res ) == 0 and len( no_cst_res ) == 0 and len( loose_res ) == 0:
        print "No -cst_res, -no_cst_res, -loose_res specified ..."
        print "  assuming that these are denovo submissions!"
        print "  All 'errors' assumed to be 8 A."
        print
        STDEV = 8.0

    if len( cst_res ) == 0:
        for m in range( 1, totres+1):
            if m not in no_cst_res:
                cst_res.append( m )

    free_res = []
    for i in range( 1, totres+1 ):
        if ( i not in cst_res ) and ( i not in loose_res ):
            free_res.append( i )

    print "Assuming default error over residues of: %6.2f A" %  STDEV
    if len( loose_res ) > 0:
        print "Assuming error over loose residues   of: %6.2f A" %  STDEV_LOOSE
    if len( free_res ) > 0:
        print "Assuming error over no_cst residues  of: %6.2f A" %  BIG_ERROR


    expected_errors = []
    for i in range( totres ):
        res = i+1
        if res in free_res:
            expected_errors.append( BIG_ERROR )
        elif res in loose_res:
            expected_errors.append( STDEV_LOOSE )
        else:
            expected_errors.append( STDEV )

################################################
# Create CASP format file
################################################

for i in range( len( files ) ):
    file = files[ i ]
    lines = open( file ).readlines()

    casp_file = 'TS%3s_%s.casp' % (target_num,i+1)
    fid = open( casp_file, 'w' )
    fid.write( 'PFRMAT TS\n')
    fid.write( 'TARGET T%04d\n' % target_num )
    fid.write( 'AUTHOR %s\n' % author_code )
    fid.write( 'METHOD Refinement and building of de novo 3D structures in small steps, while sampling under a high resolution energy function. Unpublished and largely untested.\n' )
    fid.write( 'MODEL %d\n' % (i+1) )
    fid.write( 'PARENT N/A\n' )

    for line in lines:
        if line.count('ATOM'):
            resnum = int( line[22:26] )
            new_line = '%60s%6.2f\n' % ( line[:60],expected_errors[ resnum-1 ] )
            fid.write( new_line )

    fid.write( 'TER\n' )
    fid.write( 'END\n' )
    fid.close()
    print "Made: ", casp_file
