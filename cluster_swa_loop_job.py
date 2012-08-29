#!/usr/bin/python

from sys import argv
from parse_options import parse_options
from os import system, popen
from os.path import basename,exists
from make_tag import make_tag
from glob import glob
import string

native_pdb = parse_options( argv, "native", "" )
nstruct = parse_options( argv, "nstruct", 5 )

assert( len( argv) == 4 )
loop_silent_file =  argv[1]
loop_start = int( argv[2] )
loop_stop  =  int( argv[3] )


in_loop = []

EXE = '/Users/rhiju/src/rosetta_TRUNK/rosetta_source/bin/stepwise_protein_test.macosgccrelease'
DB =  '/Users/rhiju/src/rosetta_TRUNK/rosetta_database'
SCORE_DIFF_CUT = 10.0
NSTRUCT = 100

if len(native_pdb) > 0:

    rescore_against_native_silent_file = loop_silent_file.replace( '.out', '.calcRMSD.out' )
    if not exists( native_pdb ):
        print 'Could not find ', native_pdb
        exit( 0 )

    if not exists( rescore_against_native_silent_file ):
        command = '%s  -calc_rms -in:file:silent %s  -native %s  -calc_rms_res %s  -out:file:silent %s -database %s' % ( EXE, loop_silent_file, native_pdb, make_tag( range( loop_start, loop_stop+1) ),rescore_against_native_silent_file, DB )
        print command
        system( command )

    loop_silent_file = rescore_against_native_silent_file


cluster_silent_file = basename(loop_silent_file).replace( '.out', '.cluster1.0A.out' )

if not exists( cluster_silent_file ):

    command =  '%s  -database %s -in:file:silent %s -calc_rms_res %d-%d -cluster_test -cluster:radius 1.0 -out:file:silent %s -score_diff_cut %8.3f -in:file:silent_struct_type binary -nstruct %d' % ( EXE, DB, loop_silent_file, loop_start, loop_stop, cluster_silent_file, SCORE_DIFF_CUT, NSTRUCT )
    print command
    system( command )


loop_pdb = cluster_silent_file + '.1.pdb'
if not exists( loop_pdb ):
    command = 'extract_lowscore_decoys.py %s %d' % (cluster_silent_file,nstruct)
    print command
    system( command )


all_score_gap = []
all_best_rms = []
all_best_cluster_num = []
all_top_score_rms = []
all_tag = []
all_top_score = []

print
print '==================================='
plines = popen( 'grep SCORE '+cluster_silent_file).readlines()

n_less_than_score_cut =  len( plines )-1
n_less_than_score_cut2 = 0
score_min  = 0
TIGHT_SCORE_CUT = 2.5
scores = []
for pline in plines[1:]:
    score = float( string.split( pline )[1] )
    scores.append( score )
    if score_min == 0: score_min = score
    if ( score <= score_min + TIGHT_SCORE_CUT ): n_less_than_score_cut2 += 1

scores.sort()
score_gap = 999
if len( scores ) > 1: score_gap = scores[1] - scores[0]

loop_tag = 'Loop %d-%d' % (loop_start,loop_stop)
print '%s   n<%5.1f: %2d    n<%5.1f: %2d   score_gap:%5.2f' % ( loop_tag, SCORE_DIFF_CUT, n_less_than_score_cut , TIGHT_SCORE_CUT, n_less_than_score_cut2, score_gap )
print '==================================='

cols = plines[0].split()
col_idx = []
try:
    fields = ['score','all_rms','backbone_rms','rms']
    for field in fields: col_idx.append( cols.index( field ) )
except:
    fields = ['score','looprms','loopcarms']
    for field in fields: col_idx.append( cols.index( field ) )

best_rms = 9999
best_cluster_num = 0
top_score_rms = 0
for c in range( nstruct ):
    if c+1 >= len( plines ): continue
    pline = plines[ c+1 ]
    cols = pline.split()
    for i in col_idx:
        print '%12s' % cols[i],
    print

    try:
        rms =  float( cols[ col_idx[-1] ] )
    except:
        continue

    if rms < best_rms:
        best_rms = rms
        best_cluster_num = c+1
    if ( top_score_rms == 0 ):
        top_score_rms = rms

all_score_gap.append( score_gap )
all_best_rms.append( best_rms )
all_best_cluster_num.append( best_cluster_num )
all_top_score_rms.append( top_score_rms )
all_top_score.append( scores[0] )
all_tag.append( loop_tag )


