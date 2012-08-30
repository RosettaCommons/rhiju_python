#!/usr/bin/python
from sys import argv,stdout
import string
from os.path import exists,basename
from os import system
from read_pdb import read_pdb

def Help():
    print
    print argv[0]+ ' <input pdb OR list of input pdbs>'
    print
    print ' Runs 3DNA utilities find_pair & analyze, and then goes '
    print '  through the outfile to pull out useful parameters'
    print '  in file with name  *base_pair_step_stats.txt and *base_pair_stats.txt.  '
    print
    print 'R. Das, 2012'

if len (argv )< 2:
    Help()
    exit()

file_list = argv[1]

if file_list[-4:] == '.pdb':
    X3DNA_outfiles = [ file_list ]
else:
    X3DNA_outfiles = map( lambda x:x[:-1], open( file_list ).readlines() )

file_out1 = file_list.replace('.txt','.list').replace('.pdb','.list').replace('.list','_base_pair_stats.txt' )
file_out2 = file_out1.replace( 'base_pair_stats', 'base_pair_step_stats' )

fid_out1 = open( file_out1, 'w' )
fid_out2 = open( file_out2, 'w' )

rnanum = { 'a':1, 'c':2, 'g':3, 'u':4 }

def advance_to_tag( fid, tag ):
    line = fid.readline()
    while line:
        # Let's get to base pair definition
        if len( line ) > 1:
            cols = string.split( line )
            if len( cols ) > 1 and cols[0] == tag: break
        line = fid.readline()
    line = fid.readline()
    return line

all_loop_info = {}
for file in X3DNA_outfiles:

    outfile = file
    PDB_READIN = 0

    if file[-4:] == '.pdb':
        PDB_READIN = 1

        inp_file = file.replace( '.pdb','.inp' )
        if not exists( inp_file):
            command = 'find_pair %s %s' % (file, inp_file )
            print command
            system( command )
        assert( exists( inp_file ) )

        outfile  = basename( file.replace( '.pdb','.out' ) )
        if not exists( outfile ):
            command = 'analyze %s'  % inp_file
            print command
            system( command )

    if not exists( outfile ):continue
    [ coords, pdb_lines, sequence ] = read_pdb( file )

    print outfile
    fid = open( outfile )

    wc_bps = []
    all_res = []
    all_chain = []
    wc_bp_steps = []
    seq1 = []
    seq2 = []

    line = advance_to_tag( fid, 'Strand' )

    while len( line ) > 4 and line[0] != '*' :
        n = int( line[:4] ) # base pair number

        bp_tag = line[34:39]
        if ( bp_tag == '-----'):  wc_bps.append( n )
        if ( bp_tag == 'x----'):  wc_bps.append( n )
        if ( bp_tag == '----x'):  wc_bps.append( n )

        seq1.append( line[33].lower() )
        seq2.append( line[39].lower() )

        res1 = int( line[22:26].replace( '.','') )
        res2 = int( line[46:50].replace( '.','') )
        chain1 = line[20]
        chain2 = line[52]

        if res2 > res1 :
            all_res.append( [res1,res2] )
            all_chain.append( [chain1,chain2] )
        else:
            all_res.append( [res2,res1] )
            all_chain.append( [chain2,chain1] )

        line = fid.readline()


    # find loops -- can't have anything nested inside. and must be part of a single chain.
    tag = file.replace( ".pdb","")
    for n in range( len( all_res ) ):

        if not (n+1) in wc_bps: continue

        if ( all_chain[n][1] != all_chain[n][0] ): continue
        #print all_res[n]

        nested_base_pair_exists = False
        for m in range( len( all_res ) ): # This is totally inefficient...
            if not (m+1) in wc_bps: continue
            if all_chain[m][0] != all_chain[n][0]: continue
            if all_chain[m][1] != all_chain[n][1]: continue
            if (all_res[m][0] > all_res[n][0]) and (all_res[m][1] < all_res[n][1]):
                nested_base_pair_exists = True
                break
        if nested_base_pair_exists: continue

        loop_length = all_res[n][1] - all_res[n][0] - 1
        if ( all_chain[n][1] == all_chain[n][0] ):
            selection_string = "%s and chain %s and resi %d+%d"  % (tag,  all_chain[n][0], all_res[n][0],all_res[n][1] )
        else:
            selection_string = "%s and (chain %s and resi %d) or (chain %s and resi %d)"  % (tag, all_chain[n][0], all_res[n][0], all_chain[n][1], all_res[n][0] )

        if loop_length not in all_loop_info.keys(): all_loop_info[ loop_length ] = []

        sequence_string = ""
        for m in range( all_res[n][0], all_res[n][1]+1 ):
            sequence_string += sequence[ all_chain[n][0] ][ m ]

        all_loop_info[ loop_length].append( [ selection_string, sequence_string ] )

        print "Loop length %d  %s %s" % (loop_length, selection_string, sequence_string )

for loop_length in all_loop_info.keys():
    loop_file = "loops%d.txt" % loop_length
    print "Writing to: ", loop_file
    fid = open( loop_file, 'w' )
    for loop_info in all_loop_info[ loop_length ]:
        fid.write( "%40s %s\n" % ( loop_info[0], loop_info[1] ) )
    fid.close()
