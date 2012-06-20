#!/usr/bin/python
from sys import argv,stdout
import string
from os.path import exists,basename
from os import system

def Help():
    print
    print argv[0]+ ' <input pdb OR list of input pdbs>'
    print
    print ' Runs 3DNA utilities find_pair & analyze, and then goes '
    print '  through the outfile to pull out useful parameters'
    print '  in file with name  *base_pair_step_stats.txt and *base_pair_stats.txt.  '
    print
    print ' File format of base_pair_step_stats has one line for each base pair step (i,j) to '
    print '      (i+1,j-1) with the following parameters:'
    print
    print '   type(i) type(i+1) type(j) type(j-1)    [a,c,g,u = 1,2,3,4]'
    print '   resnum(i) resnum(i+1) resnum(j) resnum(j-1) '
    print '   shift slide rise tilt roll twist'
    print '   delta(i) epsilon(i) zeta(i) alpha(i+1) beta(i+1) gamma(i+1)   chi(i) chi(i+1)'
    print '   delta(j-1) epsilon(j-1) zeta(j-1) alpha(j) beta(j) gamma(j)   chi(j-1) chi(j)'
    print '   shear(i,j) stretch(i,j) stagger(i,j) buckle(i,j) propeller(i,j) opening(i,j)'
    print '   shear(i+1,j-1) stretch(i+1,j-1) stagger(i+1,j-1) buckle(i+1,j-1) propeller(i+1,j-1) opening(i+1,j-1)'
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

    print outfile
    fid = open( outfile )

    wc_bps = []
    all_res = []
    wc_bp_steps = []
    all_local_base_pair_params = []
    all_local_base_pair_step_params = []
    all_torsions1 = []
    all_torsions2 = []
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

        all_res.append( [res1,res2] )

        line = fid.readline()

    line = advance_to_tag( fid, 'bp' )
    line = advance_to_tag( fid, 'bp' )
    while len( line ) > 12 and line[11] != '~':
        cols = string.split( line )
        all_local_base_pair_params.append( cols[2:] )
        line = fid.readline()

    line = advance_to_tag( fid, 'step' )
    while len( line ) > 12 and line[11] != '~':
        cols = string.split( line )
        all_local_base_pair_step_params.append( cols[2:] )
        line = fid.readline()

    line = advance_to_tag( fid, 'base' )
    while len( line ) > 12:
        cols = string.split( line )
        all_torsions1.append( cols[2:] )
        line = fid.readline()

    line = advance_to_tag( fid, 'base' )
    while len( line ) > 12:
        cols = string.split( line )
        all_torsions2.append( cols[2:] )
        line = fid.readline()


    # Now go through and pick out WC/WC base steps
    for n in range( len( all_res ) - 1):

        if not (n+1) in wc_bps: continue
        if not (n+2) in wc_bps: continue

        if not (all_res[n][0] == all_res[n+1][0]-1): continue
        if not (all_res[n][1] == all_res[n+1][1]+1): continue

        try:
            for k in range(6): check_if_parameter_defined = float( all_local_base_pair_step_params[n][k] )
            for k in range(3,7): check_if_parameter_defined = float( all_torsions1[n][k] )
            for k in range(4): check_if_parameter_defined = float( all_torsions2[n][k] )
        except:
            continue

        fid_out2.write( '%d '    % rnanum[ seq1[n]   ] )
        fid_out2.write( '%d    ' % rnanum[ seq1[n+1] ] )
        fid_out2.write( '%d '    % rnanum[ seq2[n+1] ] )
        fid_out2.write( '%d    ' % rnanum[ seq2[n]   ] )

        fid_out2.write( '%4d '    % all_res[n][0]    )
        fid_out2.write( '%4d    ' % all_res[n+1][0]  )
        fid_out2.write( '%4d '    % all_res[n+1][1]  )
        fid_out2.write( '%4d    ' % all_res[n][1]    )

        fid_out2.write( '    ' )

        # roll, twist, etc.
        for k in range( 6 ):
            fid_out2.write( ' %6s ' % all_local_base_pair_step_params[n][k] )

        fid_out2.write( '    ' )

        # suite torsions. Strand 1
        for k in [3,4,5]:   fid_out2.write( ' %6s ' % all_torsions1[n]  [k] )
        for k in [0,1,2,3]: fid_out2.write( ' %6s ' % all_torsions1[n+1][k] )

        fid_out2.write( '    ' )

        # suite torsions. Strand 2
        for k in [3,4,5]:   fid_out2.write( ' %6s ' % all_torsions2[n+1]  [k] )
        for k in [0,1,2,3]: fid_out2.write( ' %6s ' % all_torsions2[n][k] )

        fid_out2.write( '    ' )

        # chi angles
        fid_out2.write( ' %6s ' % all_torsions1[n]  [6] )
        fid_out2.write( ' %6s ' % all_torsions1[n+1][6] )
        fid_out2.write( '   ' )
        fid_out2.write( ' %6s ' % all_torsions2[n+1][6] )
        fid_out2.write( ' %6s ' % all_torsions2[n][6] )


        fid_out2.write( '      ' )

        # buckle, etc.
        for k in range( 6 ):
            fid_out2.write( ' %6s ' % all_local_base_pair_params[n][k] )
        fid_out2.write( '   ' )
        for k in range( 6 ):
            fid_out2.write( ' %6s ' % all_local_base_pair_params[n+1][k] )

        fid_out2.write( '\n' )

    for n in range( len( all_res ) - 1):
        if not (n+1) in wc_bps: continue

        fid_out1.write( '%d '    % rnanum[ seq1[n]   ] )
        fid_out1.write( '%d    ' % rnanum[ seq2[n] ] )

        fid_out1.write( '%4d '    % all_res[n][0]    )
        fid_out1.write( '%4d    ' % all_res[n][1]  )

        # buckle, etc.
        for k in range( 6 ):
            fid_out1.write( ' %6s ' % all_local_base_pair_params[n][k] )

        fid_out1.write( '\n' )

    if PDB_READIN: system( 'rm -rf '+outfile )


fid_out1.close()
fid_out2.close()

if PDB_READIN: system( 'rm -rf auxiliary.par bestpairs.pdb bp_helical.par bp_order.dat bp_step.par cf_7methods.par col_chains.scr col_helices.scr hel_regions.pdb hstacking.pdb poc_haxis.r3d ref_frames.dat stacking.pdb'  )

print 'Outputted: ', file_out1
print 'Outputted: ', file_out2


system( 'cat ' + file_out2 )
