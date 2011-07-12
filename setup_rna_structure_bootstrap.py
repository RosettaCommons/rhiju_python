#!/usr/bin/python

from sys import argv
from os import system
from os.path import exists,dirname,basename,expanduser
from random import randint
from parse_options import parse_options
import string

HOMEDIR = expanduser('~rhiju')

EXE = HOMEDIR+'/projects/rdat/external/RNAstructure_Rochester/RNAstructure/exe/Fold'
#EXE = HOMEDIR+'/projects/rdat/external/RNAstructure_DASLAB/exe/Fold'
assert( exists(EXE ) )

PYDIR = HOMEDIR+'/python'
assert( exists( PYDIR ) )

seq_file = parse_options( argv, 'seq_file', '' )
outdir = parse_options( argv, 'outdir','RUNS')
min_res = parse_options( argv, 'min_res', -1 )
max_res = parse_options( argv, 'max_res', -1 )
sm = parse_options( argv, 'sm', -99.9 )
si = parse_options( argv, 'si', -99.9 )
maxdistance = parse_options( argv, 'maxdistance', -1 )
shape_file = parse_options( argv, 'shape_file', '' )
force_unpair_file = parse_options( argv, 'force_unpair_file', '' )
nboot = parse_options( argv, 'nboot', 0 )
jobfile = parse_options( argv, 'jobfile', 'bsubRNAstructure')

assert( len( argv ) == 1 )

force_unpair = []
if len( force_unpair_file ) > 0:
    lines = open( force_unpair_file ).readlines()
    for line in lines:
        for m in string.split( line ):
            force_unpair.append( int(m) )

if (min_res != 1 or max_res != 1 ):
    system( '%s/slice_seq.py %s %d %d' % (PYDIR,seq_file,min_res,max_res) )
    seq_file = 'region_%d_%d_%s' % (min_res,max_res,seq_file)
    if len( shape_file ) > 0:
        system( '%s/slice_shape.py %s %d %d' % (PYDIR,shape_file,min_res,max_res) )
        shape_file = 'region_%d_%d_%s' % (min_res,max_res,shape_file)

    if len( force_unpair) > 0:
        force_unpair_original = force_unpair
        force_unpair = []
        for m in force_unpair_original:
            if ( m >= min_res and m <= max_res ):
                force_unpair.append( m - min_res + 1)

if ( maxdistance > 0 ): EXE += ' -md %d ' % maxdistance

assert( len( seq_file ) > 0 )
if nboot > 0: assert( len( shape_file ) > 0 )

if not exists( outdir ): system( 'mkdir -p '+outdir )

def outputcommand( seq_file, ct_file, shape_file, force_unpair_constraints_file, fid ):
    shape_tag = ''
    if len( shape_file ) > 0: shape_tag = ' -sh %s' % shape_file
    command_line = '%s %s %s %s' % ( EXE, seq_file, ct_file, shape_tag )

    if len( force_unpair_constraints_file ) > 0:
        command_line += ' -c %s ' % force_unpair_constraints_file

    if ( sm != -99.9 ): command_line += ' -sm %5.2f ' % sm
    if ( si != -99.9 ): command_line += ' -si %5.2f ' % si

    outfile = '/dev/null'
    errfile = '/dev/null'
    bsub_command =  'bsub -W 16:0 -o %s -e %s %s \n' % (outfile, errfile, command_line )
    fid.write( bsub_command )

force_unpair_constraints_file = ''
if len( force_unpair ) > 0:
    force_unpair_constraints_file = 'region_%d_%d_%s' % (min_res, max_res, 'force_unpair.con' )
    fid = open( force_unpair_constraints_file, 'w' )
    fid.write( 'DS:\n')
    fid.write( '-1\n')
    fid.write( 'SS:\n')
    for m in force_unpair: fid.write( '%d\n' % m)
    fid.write( '-1\n')
    fid.write( 'Mod:\n')
    fid.write( '-1\n')
    fid.write( 'Pairs:\n')
    fid.write( '-1 -1\n')
    fid.write( 'FMN:\n')
    fid.write( '-1\n')
    fid.write( 'Forbids:\n')
    fid.write( '-1 -1\n')


fid = open( jobfile, 'w' )
ct_file = outdir+'/'+seq_file+'.ct'
if not exists( ct_file ):
    outputcommand( seq_file, ct_file, shape_file, force_unpair_constraints_file, fid )
else:
    print 'Job done already ==> ', ct_file

if len( shape_file ) > 0 and nboot > 0:
    lines = open( shape_file ).readlines()
    numlines = len( lines )

    nboot_actual = 0
    for n in range( nboot ):

        boot_ct_file = ct_file+'.boot%05d' % n

        #print boot_ct_file

        if not exists( boot_ct_file ):
            rand_index = []
            for m in range( numlines ): rand_index.append(  randint(0,numlines-1) )
            rand_index.sort()

            boot_shape_file = outdir + '/'+shape_file+'.boot%05d' % n
            fid2 = open( boot_shape_file, 'w' )
            for m in range( numlines ):
                fid2.write( lines[ rand_index[m] ] )
            fid2.close()

            outputcommand( seq_file, boot_ct_file, boot_shape_file, force_unpair_constraints_file, fid )
            nboot_actual += 1

    print 'Set up %d bootstraps for %s' % (nboot_actual,seq_file)
