#!/usr/bin/python

from sys import argv
from os.path import exists, basename
from os import system, chdir, getcwd
from glob import glob
from os import popen
import string

job_list = argv[1]
lines = open( job_list ).readlines()

assert( job_list.find( '.txt' ) > 0 )
loop_model_dir = job_list.replace( '.txt', '_clusters' )
if not exists( loop_model_dir ): system( 'mkdir '+loop_model_dir )

def make_tag( int_vector ):
    tag = ''
    for m in int_vector: tag += ' %d' % m
    return tag

CWD = getcwd()

print 'Working directory:', CWD
print 'Script command', string.join( argv )


# Changed later -- originally used SCORE_DIFF_CUT = 5;
SCORE_DIFF_CUT = 40;
NSTRUCT = 20;

for line in lines:
    loop_tag = line[:-1]

    chdir( loop_tag )

    loop_silent_file = 'region_FINAL.out'

    sequence = open( '%s.fasta' % loop_tag ).readlines()[-1][:-1]
    in_loop = []
    for m in range( len( sequence ) ): in_loop.append( 0 )

    loop_file = '%s.loop'  % loop_tag
    cols = open( loop_file ).readlines()[0].split()
    loop_start = int( cols[0] )
    loop_stop = int( cols[1] )

    native_pdb = '%s_min.pdb' % loop_tag
    assert( exists( native_pdb ) )

    cluster_silent_file = basename(loop_silent_file).replace( '.out', '.cluster1.0A.out' )

    if not exists( cluster_silent_file ):

        input_res = []
        for m in range( len( sequence ) ):
            if not in_loop[ m ] or (m+1 >= loop_start) or (m+1 <= loop_stop ): input_res.append( m+1 )

        EXE = '/home/rhiju/src/mini/bin/stepwise_protein_test.linuxgccrelease'
        #assert( exists(EXE) )
        if not exists( EXE ): EXE =  '/Users/rhiju/src/mini/bin/stepwise_protein_test.macosgccrelease'
        command =  '%s  -database ~/minirosetta_database -in:file:silent %s -calc_rms_res %d-%d -cluster_test -cluster:radius 1.0 -out:file:silent %s -score_diff_cut %8.3f -in:file:silent_struct_type binary -nstruct %d' % ( EXE, loop_silent_file, loop_start, loop_stop, cluster_silent_file, SCORE_DIFF_CUT, NSTRUCT )
        print command
        system( command )

    chdir(CWD )

    loop_cluster_file = '%s_FINAL.cluster1.0A.out' % loop_tag
    loop_silent_file_copy = '%s/%s' % ( loop_model_dir, loop_cluster_file )
    if not exists( loop_silent_file_copy ):
        command = 'rsync -avz %s/%s  %s' % (loop_tag, cluster_silent_file, loop_silent_file_copy )
        print command
        system( command )

    chdir( loop_model_dir )
    loop_pdb = loop_cluster_file + '.1.pdb'
    if not exists( loop_pdb ):
        command = 'extract_lowscore_decoys.py %s 5' % loop_cluster_file
        print command
        system( command )

    chdir( CWD )


all_score_gap = []
all_best_rms = []
all_top_score_rms = []
all_tag = []

chdir( loop_model_dir )
for line in lines:
    loop_tag = line[:-1]
    loop_cluster_file = '%s_FINAL.cluster1.0A.out' % loop_tag
    print
    print '==================================='
    plines = popen( 'grep SCORE '+loop_cluster_file).readlines()

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
    score_gap = scores[1] - scores[0]

    print '%s   n<%5.1f: %2d    n<%5.1f: %2d   score_gap:%5.2f' % ( loop_tag, SCORE_DIFF_CUT, n_less_than_score_cut , TIGHT_SCORE_CUT, n_less_than_score_cut2, score_gap )
    print '==================================='

    fields = ['score','all_rms','backbone_rms','rms']
    cols = plines[0].split()
    col_idx = []
    for field in fields: col_idx.append( cols.index( field ) )

    best_rms = 9999
    top_score_rms = 0
    for pline in plines[:6]:
        cols = pline.split()
        for i in col_idx:
            print '%12s' % cols[i],
        print


        try:
            rms =  float( cols[ col_idx[-1] ] )
        except:
            continue

        if rms < best_rms: best_rms = rms
        if ( top_score_rms == 0 ): top_score_rms = rms

    all_score_gap.append( score_gap )
    all_best_rms.append( best_rms )
    all_top_score_rms.append( top_score_rms )
    all_tag.append( loop_tag )

chdir( CWD )


cluster_summary_file = job_list.replace( '.txt', '_cluster_summary.txt' )
print
print 'Making ... ', cluster_summary_file
fid = open( cluster_summary_file, 'w' )



fid.write( '%4s  %6s  %6s  %6s\n' % ( 'ID','rms1','rms5','Egap'));
for i in range( len( all_tag ) ):
    fid.write( '%4s  %6.2f  %6.2f  %6.2f\n' % ( all_tag[i], all_top_score_rms[i], all_best_rms[i], all_score_gap[i]) );

fid.close()
