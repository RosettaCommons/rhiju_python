#!/usr/bin/python
from sys import argv,stdout
import string
from os.path import exists,basename
from os import system
from numpy import array, cross, mat
from numpy.linalg import norm
from read_pdb import read_pdb

file_list = argv[1]
outfile_out = file_list + '.out'
fid_out = open( outfile_out, 'w' )

if file_list[-4:] == '.pdb':
    files = [ file_list ]
else:
    files = map( lambda x:x[:-1], open( file_list ).readlines() )

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

for file in files:

    assert( file[-4:] == '.pdb' )

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
    seq1 = []
    seq2 = []
    res1 = []
    res2 = []
    chain1 = []
    chain2 = []
    origins = []
    normals = []

    line = advance_to_tag( fid, 'Strand' )

    while len( line ) > 4 and line[0] != '*' :
        n = int( line[:4] ) # base pair number

        bp_tag = line[34:39]
        if ( bp_tag == '-----'):  wc_bps.append( n )
        if ( bp_tag == 'x----'):  wc_bps.append( n )
        if ( bp_tag == '----x'):  wc_bps.append( n )

        seq1.append( line[33].lower() )
        seq2.append( line[39].lower() )

        chain1.append( line[20] )
        chain2.append( line[52] )

        res1.append( int( line[22:26].replace( '.','') ) )
        res2.append( int( line[46:50].replace( '.','') ) )

        line = fid.readline()

    line = advance_to_tag( fid, 'bp' )
    while len( line ) > 12 and line[11] != '*':
        cols = string.split( line )
        origins.append( array( [ float(cols[2]), float(cols[3]), float(cols[4]) ] ))
        normals.append( array( [ float(cols[5]), float(cols[6]), float(cols[7]) ] ))
        line = fid.readline()

    # need to write this reader!   coords[chain][res][atom_name] = [ x, y, z]
    coords = read_pdb( file )

    for m in range( len(wc_bps)-1 ):

        # offset by one, 3DNA indexes by 1.
        i = wc_bps[m]-1
        j = wc_bps[m+1]-1

        # check that these are on the same chain.
        if ( chain1[i]  != chain1[j] ): continue
        if ( chain2[i]  != chain2[j] ): continue

        # figure out sequence separations.
        seqsep1 = res1[j] - res1[i]
        seqsep2 = res2[i]   - res2[j]

        # define coordinate system of this base pair
        # I've got the normal from 3DNA
        z = normals[ i ]
        # Now get a vector from C1* to C1*
        x = array( coords[ chain2[i] ][ res2[i] ][ " C1*" ] ) -  array( coords[ chain1[i] ][ res1[i] ][ " C1*" ] )

        # do cross products to get x,y
        y = cross( z, x )
        y /= norm(y)
        x = cross( y, z)
        x /= norm(x)

        # orthonormal matrix
        m = mat( [x, y, z] )

        # get C1* of next base pair
        pos1 = array( coords[ chain1[j] ][ res1[j] ][ " C1*" ]  )
        pos_local1 = ( pos1 - origins[i] ) * m.H

        # get other C1* of next base pair
        pos2 = array( coords[ chain2[j] ][ res2[j] ][ " C1*" ] )
        pos_local2 = ( pos2 - origins[i] ) * m.H

        # output to a file for matlab visualization. Damn, probably should just use my pylab plugin. Anyway.
        fid_out.write( "%4d   %4d %4d  %4d %4d    %4d   %4d  %8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f\n"  % \
                       (i, res1[i],res2[i],  res1[i+1],res2[i+1],  seqsep1, seqsep2,  \
                            pos_local1[0,0],pos_local1[0,1],pos_local1[0,2], \
                            pos_local2[0,0],pos_local2[0,1],pos_local2[0,2]  ) )


fid_out.close()
print 'Put stats in ', outfile_out
