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
    fasta_file = glob( '*.fasta' )[0]
    sequence = open( fasta_file ).readlines()[1][:-1]
    nres = len( sequence )

    cols = open( 'README_SETUP' ).readlines()[1].split()
    pos = cols.index( '`seq' )
    loop_start = int( cols[pos+1] )
    loop_end   = int( cols[pos+2][:-1] )

    rosetta_path = '/home/rhiju/src/rosetta_TRUNK/'
    new_cat_outfile_convention = 0

    if sequence.count( 'a' ) > 0:
        rosetta_path = '/home/rhiju/src/rosetta_protein_rna/'
        #new_cat_outfile_convention = 1

    return (sequence, loop_start, loop_end, rosetta_path, new_cat_outfile_convention )

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
        command = '%s/rosetta_source/bin/extract_pdbs.linuxgccrelease -in:file:silent  %s  -in:file:silent_struct_type binary  -in:file:tags %s -database %s/rosetta_database/' % (rosetta_path, outfile, model , rosetta_path)
        print command
        system( command )
        assert( exists( model+'.pdb' ) )
        system( 'mv %s.pdb %s' % ( model, extract_model_name ) )

    command = 'renumber_pdb_in_place.py %s ' % extract_model_name
    try:
        i = int( outfile.split( '_' )[2] )
        j = int( outfile.split( '_' )[1] )
        for n in range(1,i+1): command += ' %d' % n
        for n in range(j+1,len(sequence)+1): command += ' %d' % n
        #print command
        system( command )
    except:
        print outfile

    return extract_model_name

def get_previous_model_and_outfile( parent_tag, source, new_cat_outfile_convention ):
    prev_outfile = source.split( '/' )[1]
    cols = prev_outfile.split( '_' )
    outfile = string.join( cols[2:5], '_' )+'_sample.cluster.out'
    get_outfile_ready( outfile )

    nmodel = 0
    cols = parent_tag.split( '_' )
    if len( cols ) > 2:
        nmodel = int(cols[1])
        if new_cat_outfile_convention: nmodel = int( cols[2] )

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
        command = '%s/rosetta_source/bin/stepwise_protein_test.linuxgccrelease -in:file:silent %s  -calc_rms -out:file:silent %s -native %s  -database %s/rosetta_database/ -calc_rms_res %d-%d -working_res 1-%d %d-%d' % (rosetta_path, outfile, outfile_calcRMSD, start_model_pdb_file, rosetta_path, loop_start, loop_end, i, j, len( sequence ) )
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
(sequence, loop_start, loop_end, rosetta_path, new_cat_outfile_convention) = initialize()
outfiles_backtrack = []

# starting outfile is special -- a loop closure
outfile = 'region_FINAL.out'
outfiles_backtrack.append(  (outfile, start_model, 0.0 )  )
( parent_tag, source ) = find_parent_tag_and_source(  outfile, start_model )
start_model_pdb_file = extract_model( outfile, start_model, start_model )

# kind of weird -- need to track down model inside chain-closure file.
chain_closure_outfile = source.split( '/' )[0].lower()+'_sample.cluster.out'
outfiles_backtrack.append( (chain_closure_outfile, 'ND', 0.0 ) )

# MAIN LOOP
count = 0
while source:# and len( outfiles_backtrack ) < 2:
    (nmodel, outfile ) = get_previous_model_and_outfile( parent_tag, source, new_cat_outfile_convention )
    (model,rms) = figure_out_actual_model( nmodel, outfile )
    outfiles_backtrack.append( (outfile, model, rms ) )
    # OK, let's go to the next one.
    ( parent_tag, source ) = find_parent_tag_and_source( outfile, model )
    print parent_tag, source

# What we found in the backtrack
summary_file = 'backtrack_summary_%s.txt' % start_model
fid = open( summary_file, 'w' )
for n in range( len( outfiles_backtrack )): fid.write( '%s %s %6.2f\n' % outfiles_backtrack[n] )
fid.close()

print
print 'BACKTRACK:'
system( 'cat '+summary_file )
