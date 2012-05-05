#!/usr/bin/python

from sys import argv
import string
from os import popen, system
from os.path import exists
from glob import glob

# Let's get going!
start_model = 'S_0'
if len( argv ) > 1 :
    start_model = argv[1]
    assert( start_model[:2] == 'S_' )

def initialize():
    fasta_file = glob( '????.fasta' )[0]
    sequence = open( fasta_file ).readlines()[1][:-1]
    nres = len( sequence )

    cols = open( 'README_SETUP' ).readlines()[1].split()
    pos = cols.index( '`seq' )
    loop_start = int( cols[pos+1] )
    loop_end   = int( cols[pos+2][:-1] )
    return (sequence, loop_start, loop_end)

# some useful functions
def find_parent_tag_and_source( outfile, model ):
    parent_tag = None
    source = None

    lines = popen( 'grep -A3 SCORE '+outfile ).readlines()

    found_model = 0
    for line in lines:
        cols = string.split( line[:-1] )
        if cols[-1] == model:
            found_model = 1
        if found_model and cols[0] == 'REMARK' and cols[1] == 'PARENT_TAG':
            parent_tag = cols[-1]
        if found_model and cols[0] == 'REMARK' and cols[1] == 'SOURCE':
            source = cols[-1]
            break
    return (parent_tag, source )

def get_outfile_ready( outfile ):
    if not exists( outfile ):
        assert( exists( outfile +'.gz' ) )
        system( 'gunzip ' + outfile )


def extract_model( outfile, model, start_model_name ):
    extract_model_name = '%s.backtrack_%s.pdb' % (outfile.replace('.out',''), start_model_name )
    if not exists( extract_model_name ):
        print 'EXTRACTING: ', extract_model_name
        command = '/home/rhiju/src/rosetta_TRUNK/rosetta_source/bin/extract_pdbs.linuxgccrelease -in:file:silent  %s  -in:file:silent_struct_type binary  -in:file:tags %s -database ~/rosetta_database_TRUNK/' % (outfile, model )
        print command
        system( command )
        assert( exists( model+'.pdb' ) )
        system( 'mv %s.pdb %s' % ( model, extract_model_name ) )
    return extract_model_name

def get_previous_model_and_outfile( parent_tag, source ):
    prev_outfile = source.split( '/' )[1]
    cols = prev_outfile.split( '_' )
    outfile = string.join( cols[2:5], '_' )+'_sample.cluster.out'
    get_outfile_ready( outfile )

    nmodel = 0
    cols = parent_tag.split( '_' )
    if len( cols ) > 2:
        nmodel = int(cols[1])

    return (nmodel,outfile)


def figure_out_actual_model( nmodel, outfile ):
    # rescore -- this serves two purposes
    # first, it gives back the rmsd to the desired model -- this provides a consistency check
    # second, it 'resorts' the models in the order in which they would be sorted by glob --
    #       S_0, S_1, S_10, S_100, etc.
    # and that lets us figure out which decoy to track back. Its roundabout but it works.
    outfile_calcRMSD = outfile.replace( '.out', '.backtrack_%s.out' % start_model )
    if not exists( outfile_calcRMSD ):
        i = int( outfile.split( '_' ) [2] )
        j = int( outfile.split( '_' ) [1] )
        command = ' ~/src/rosetta_TRUNK/rosetta_source/bin/stepwise_protein_test.linuxgccrelease -in:file:silent %s  -calc_rms -out:file:silent %s -native %s  -database ~/rosetta_database_TRUNK/ -calc_rms_res %d-%d -working_res 1-%d %d-%d' % ( outfile, outfile_calcRMSD, start_model_pdb_file, loop_start, loop_end, i, j, len( sequence ) )
        print( command )
        system( command )

    scorefile_calcRMSD = outfile_calcRMSD.replace('.out','.sc')
    if not exists( scorefile_calcRMSD ):
        command =  'grep SCORE %s | head -n %d > %s' % ( outfile_calcRMSD, nmodel+2, scorefile_calcRMSD )
        print command
        system( command )

    lines = open( scorefile_calcRMSD ).readlines()
    cols = lines[0].split()
    rms_col = cols.index( 'rms' )
    cols = lines[-1][:-1].split()
    rms = float( cols[ rms_col ] )
    model = cols[-1]

    print outfile_calcRMSD, model, rms
    extract_model( outfile, model, start_model ) # extract for viewing.

    return (model,rms)

# basic initialization -- what's the sequence? the loop boundaries?
(sequence, loop_start, loop_end) = initialize()
outfiles_backtrack = []

# starting outfile is special -- a loop closure
outfile = 'region_FINAL.out'
outfiles_backtrack.append(  (outfile, start_model, 0.0)  )
( parent_tag, source ) = find_parent_tag_and_source(  outfile, start_model )
start_model_pdb_file = extract_model( outfile, start_model, start_model )

# MAIN LOOP
count = 0
while parent_tag: # and len( outfiles_backtrack ) < 7:
    (nmodel, outfile ) = get_previous_model_and_outfile( parent_tag, source )
    (model,rms) = figure_out_actual_model( nmodel, outfile ) # extracts also.
    outfiles_backtrack.append( (outfile, model, rms ) )
    # OK, let's go to the next one.
    ( parent_tag, source ) = find_parent_tag_and_source( outfile, model )
    print parent_tag, source

# What we found in the backtrack
print
print 'BACKTRACK:'
for n in range( len( outfiles_backtrack )): print '%s %s %6.2f' % outfiles_backtrack[n]

