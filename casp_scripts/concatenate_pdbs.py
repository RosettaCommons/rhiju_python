#!/usr/bin/python

from sys import argv
import string
from glob import glob
from os import system
from os.path import dirname, expanduser, exists,basename
from parse_options import parse_options
from get_sequence import get_sequence
from make_tag import make_tag
from random import random


EXE = expanduser('~rhiju')+'/src/mini/bin/stepwise_protein_test.macosgccrelease'
if not( exists( EXE ) ):
    EXE = expanduser('~rhiju')+'/src/mini/bin/stepwise_protein_test.linuxgccrelease'
assert( exists( EXE ) )

EXTRACT = expanduser('~rhiju')+'/python/extract_lowscore_decoys.py'
assert( exists( EXTRACT ) )

fasta_file = parse_options( argv, "fasta", "" )
start_files = parse_options( argv, "s", [""] )
assert( len( argv ) == 1 )

sequence_lines = open( fasta_file  ).readlines()[1:]
full_sequence = string.join(  map( lambda x : x[:-1], sequence_lines) ,  '' )
NRES = len( full_sequence )

file_bounds = []
########################################################
for file in start_files:
    seq = get_sequence( file )
    pos = string.find( full_sequence, seq )
    assert( pos > -1 )
    file_bounds.append( [pos+1, pos+len(seq)] )
    print file, file_bounds[-1]

########################################################
# Need to go through and see if there are points
# where we need to do fragment insertions.
numfiles =  len( start_files )
current_start = file_bounds[0][0]
current_end   = file_bounds[0][1]

MIN_FRAG_LENGTH = 3
MAX_FRAG_LENGTH = 20

all_files = []
all_file_bounds = []
frag_bounds = []
all_current_start = []
all_current_end   = []

def  fill_with_frag_moves( frag_start, frag_end, current_start, current_end,   \
                           all_files, all_file_bounds, frag_bounds, all_current_start, all_current_end ):

    frag_length = frag_end - frag_start + 1

    if ( frag_length > MAX_FRAG_LENGTH ):
        num_frag_steps = frag_length / MAX_FRAG_LENGTH

        for k in range( num_frag_steps ):
            all_files.append( "" )
            all_file_bounds.append( [] )
            frag_bounds.append( [ frag_start, frag_start+MAX_FRAG_LENGTH-1 ] )
            all_current_start.append( current_start )
            all_current_end.append( current_end )

            current_end = frag_start + MAX_FRAG_LENGTH - 1
            print "frag ", [ frag_start, current_end ]

            frag_start = current_end + 1

    if ( frag_end - frag_start + 1 ) < MIN_FRAG_LENGTH:
        frag_start = current_end
        frag_end = frag_start + MIN_FRAG_LENGTH - 1

    return   ( frag_start, frag_end, current_end )


final_pdb_file = basename(start_files[0]).replace( '.pdb','')

for i in range( 1, numfiles ):

    final_pdb_file += '_' + basename( start_files[i] ).replace('.pdb','')

    # These have to be in order.
    assert( file_bounds[i-1][0] < file_bounds[i][0] )

    frag_start = current_end + 1
    frag_end   = file_bounds[i][0] - 1
    ( frag_start, frag_end, current_end ) = fill_with_frag_moves( frag_start, frag_end, current_start, current_end,   \
                                                                  all_files, all_file_bounds, frag_bounds, all_current_start, all_current_end )

    print "file ", [ frag_start, frag_end ]

    all_files.append( start_files[ i ] )
    all_file_bounds.append( file_bounds[ i ] )
    frag_bounds.append( [frag_start, frag_end ] )
    all_current_start.append( current_start )
    all_current_end.append( current_end )

    current_end = file_bounds[ i ][-1]

########################################################
# Add termini with this script too.
########################################################
# C-terminus
if ( current_end < NRES ):

    final_pdb_file += '_addCterm'

    frag_start = current_end + 1
    frag_end = NRES

    ( frag_start, frag_end, current_end ) = fill_with_frag_moves( frag_start, frag_end, current_start, current_end,   \
                                                                  all_files, all_file_bounds, frag_bounds, all_current_start, all_current_end )

    print "frag ", [ frag_start, frag_end ]

    all_files.append( "" )
    all_file_bounds.append( [] )
    frag_bounds.append( [ frag_start, NRES ] )
    all_current_start.append( current_start )
    all_current_end.append( current_end )
    current_end = NRES


# N-terminus (not coded up yet) .. this actually needs to creep backwards
if ( current_start > 1 ):
    final_pdb_file += '_addNterm'

    frag_start = 1
    frag_end = current_start-1

    # This does not work! Need to allow reversal of polarity.
    #( frag_start, frag_end, current_end ) = fill_with_frag_moves( frag_start, frag_end, current_end,   \
    #                                                              all_files, all_file_bounds, frag_bounds, all_current_end )

    print "frag ", [ frag_start, frag_end ]

    all_files.append( "" )
    all_file_bounds.append( [] )
    frag_bounds.append( [ frag_start,frag_end ] )
    all_current_start.append( current_start )
    all_current_end.append( current_end )
    current_start = 1



final_pdb_file += '.pdb'
print final_pdb_file

########################################################
# Go through one by one and concatenate.
prev_file = start_files[0]
random_num = int( random() * 1000000 )
for i in range( len(all_file_bounds) ):

    command = '%s -database ~rhiju/minirosetta_database ' % EXE
    command += ' -fasta %s' % fasta_file

    fixed_res = []

    current_start = all_current_start[ i ]
    current_end = all_current_end[ i ]
    input_res1 = range( current_start, current_end+1)
    command += ' -s1 %s -input_res1 %s' % ( prev_file, make_tag( input_res1 ) )
    for m in input_res1:
        if m not in fixed_res: fixed_res.append( m )

    added_file   = all_files[ i ]
    if len( added_file ) > 0:
        file_bounds = all_file_bounds[ i ]
        input_res2 = range( file_bounds[0], file_bounds[1]+1)
        command += ' -s2 %s -input_res2 %s' % ( added_file, make_tag( input_res2 ) )
        for m in input_res2:
            if m not in fixed_res: fixed_res.append( m )

    frag_start = frag_bounds[ i ][ 0 ]
    frag_end   = frag_bounds[ i ][ 1 ]
    frag_length = frag_end - frag_start + 1

    assert( frag_length >= MIN_FRAG_LENGTH )
    assert( frag_length <= MAX_FRAG_LENGTH )

    fasta_dir = dirname( fasta_file )
    globfiles = glob( fasta_dir + '/*aa*%02d_05.200_v1_3*' % frag_length )
    assert( len( globfiles ) == 1 )
    frag_file = globfiles[ 0 ]

    for m in range( frag_start, frag_end+1):
        if m in fixed_res:
            pos = fixed_res.index( m )
            del( fixed_res[ pos ] )

    command += ' -in:file:frag_files %s' % ( frag_file )
    command += ' -sample_res %s' % ( make_tag( range( frag_start, frag_end+1) ) )

    command += ' -nstruct 1'

    #command += ' -centroid' # hmmm...

    outfile = 'temp%d_%d.out' % (random_num, i)
    command += ' -out:file:silent %s' % outfile

    command += ' -fixed_res %s' % make_tag( fixed_res )

    command += ' -pack_weights pack_no_hb_env_dep_strongrg.wts'
    command += ' -score:weights score12_no_hb_env_dep.wts'

    print 'CURRENT_END: ', current_end, '   FRAG: ',frag_start, frag_end

    print command
    system( command )

    ###########
    command = EXTRACT + ' %s 1' % outfile
    print command
    system( command )

    prev_file = '%s.1.pdb' %  outfile

command =  'mv ' + prev_file + ' '+final_pdb_file
print command
system( command )

system( 'rm -rf temp%d* '  % random_num )

