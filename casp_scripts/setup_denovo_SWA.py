#!/usr/bin/python

from sys import argv
from parse_options import parse_options
from os.path import exists, expanduser
from os import popen,system
from glob import glob
from make_tag import make_tag
import string

native = parse_options( argv, "native", "" )
fasta = parse_options( argv, "fasta", "" )
cluster_radius = parse_options( argv, "cluster_radius", 0.5 )
final_number = parse_options( argv, "final_number", 50 )
cst_file = parse_options( argv, "cst_file", "" )
optimal_length = parse_options( argv, "optimal_length", 10 )
pathway = parse_options( argv, "pathway", [-1] )

PYDIR = expanduser( "~rhiju" ) + "/python"
assert( exists( PYDIR ) )

assert( len( argv ) == 1 )

if len( fasta ) == 0:
    globfiles = glob( "*.fasta" )
    assert( len( globfiles ) == 1 )
    fasta = globfiles[ 0 ]

globfiles = glob( "*aa*v1_3*" )
assert( len( globfiles ) > 0 )
frag_files = globfiles

frag_lengths = []
for frag_file in frag_files:
    pos = frag_file.find("_05.200_v1_3" )
    assert( pos > 0 )
    frag_length = int( frag_file[pos-2 : pos] )
    frag_lengths.append( frag_length )

frag_lengths.sort()
#min_length = frag_lengths[ 0 ]

#
# Currently assumes a spacing bridged by 9mers.
# Later allow for arbitrary spacing and optimize
# stepping stones based on secondary structure.
#
ASSUMED_SPACING = 9
assert( ASSUMED_SPACING in frag_lengths )

fid = open( "README_SETUP", "w" )
fid.write(  "rm -rf REG* *~ CONDOR core.* SLAVE* \n" )

native_tag = ""
if len( native ) > 0:
    assert( exists( native ) )
    native_tag = " -native %s" % native

assert( exists( fasta ) )
sequence_lines = open( fasta  ).readlines()[1:]
sequence = string.join(  map( lambda x : x[:-1], sequence_lines) ,  '' )

#endpoints = []
#for i in range( 1 + (len( sequence ) - 1 )/ (ASSUMED_SPACING-1) ):
#    endpoints.append( 1 + i * (ASSUMED_SPACING-1) )
#print endpoints
#print len( sequence )
#if ( len( sequence ) not in endpoints): endpoints.append( len( sequence ) )
#print "Last frag spacing: ", (endpoints[-1] - endpoints[-2] + 1)

if len( pathway ) > 0:
    pathway_file = "guess.pathway"
    command = PYDIR+"/generate_pathway.py "
    for m in pathway: command += ' %d' % m
    command += " > "+pathway_file
    system( command )
    assert( exists( pathway_file ) )

psipred_file = fasta.replace( '.fasta', '.psipred_ss2' )
optimal_length_tag = ' -optimal_length %d ' % optimal_length
endpoint_tag = popen( 'endpoint_picker.py ' + optimal_length_tag + ' '+ psipred_file ).readlines()[-1][:-1]

command = "grinder_dagman.py %s -fasta %s -cluster_radius %6.2f  -final_number %d  -frag_files %s %s" % \
    ( native_tag, fasta, cluster_radius, final_number, string.join( frag_files ), endpoint_tag )

if len( cst_file ) > 0 : command += " -cst_file " + cst_file
if len( pathway ) > 0: command += " -pathway_file " + pathway_file

fid.write( command + '\n' )
fid.close()


fid= open( 'README_SUB','w' )
fid.write( 'rm -rf blah.*\n' )
fid.write( 'bsub -W 48:0 -o blah.out -e blah.err SWA_pseudo_dagman_continuous.py -j 50 protein_build.dag\n' )
fid.close()
