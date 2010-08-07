#!/usr/bin/python

from sys import argv
from os.path import exists, expanduser
from os import system
from get_sequence import get_sequence
from math import exp,sqrt

PYDIR = expanduser('~rhiju')+'/python/'
assert( exists( PYDIR) )

file1 = argv[1]
other_files = argv[2:]

def get_dist2( v1, v2 ):
    dx = v1[0] - v2[0]
    dy = v1[1] - v2[1]
    dz = v1[2] - v2[2]
    return ( dx*dx + dy*dy + dz*dz )

def get_score( dist2 ):
    CUTOFF2 = 6.0 * 6.0
    # bonus for distances less than cutoff; penalties for longer distances.
    score = exp( -1.0 * (dist2/CUTOFF2) ) #- exp( -1.0 )
    #score = ( dist2 < CUTOFF2 )
    return score


seq1 = get_sequence( file1 )

for file2 in other_files:
    seq2 = get_sequence( file2 )

    tmp_sup = 'tmp.pdb'
    command = PYDIR+'/superimpose.py '+file1+' '+file2+' > '+tmp_sup
    print command
    system( command )

    lines = open( tmp_sup ).readlines()
    count = 0
    rescount = 0
    prev_resnum = ''
    model_xyzs = [ [], [] ]
    totres = []
    for line in lines:
        if len( line ) > 5 and line[:6] == 'ENDMDL':
            totres.append( rescount )
            count += 1
            rescount = 0

        if  len(line) < 40: continue

        resnum = line[22:26]
        if not resnum == prev_resnum:
            rescount += 1

        if ( line[12:16]==' CA ' ): model_xyzs[ count ].append( [float(line[30:38]), float(line[38:46]), float(line[46:54])] )

        prev_resnum = resnum

    totres.append( rescount )

    print 'Lengths: ', totres[0], ' and ', totres[1]
    assert( len( model_xyzs[ 0 ] ) == totres[ 0 ] )
    assert( len( model_xyzs[ 1 ] ) == totres[ 1 ] )
    assert( len( seq1 ) == totres[ 0 ] )
    assert( len( seq2 ) == totres[ 1 ] )

    print 'Calculating distances ... '

    dist2 = []
    for i in range( totres[0] ):
        dist2.append( [] )
        for j in range( totres[1] ):
            dist2[ i ].append( get_dist2( model_xyzs[0][i], model_xyzs[1][j] ) )

    # consistency check
    for i in range( totres[0] ):
        for j in range( totres[1] ):
            if ( j == 0 or dist2[i][j] < mindist ):
                mindist = dist2[i][j]
                best_j = j
        #print i+1, best_j+1, sqrt( mindist)

    #print map( lambda x:int(x), dist2[ 42 ] )

    # Dynamic programming matrix
    DP_score = []
    # Boundary conditions
    for i in range( totres[0]+1):
        DP_score.append([])
        for j in range( totres[1]+1):
            DP_score[i].append( 0.0 )


    # Fill in matrix, left to right, top to bottom.

    #initialize
    choice = []
    for i in range( totres[0]+1 ):
        choice.append( [] )
        for j in range( totres[1]+1 ):
            choice[i].append( [] )

    for i in range( 1, totres[0] + 1 ):
        choice[i][0] = [i-1,0]

    for j in range( 1, totres[1] + 1 ):
        choice[0][j] = [0,j-1]

    for i in range( 1, totres[0]+1 ):

        for j in range( 1, totres[1]+1 ):

            alternatives = []

            # stupid off by one's in dist2 matrix...
            score = DP_score[i-1][j-1] + get_score( dist2[i-1][j-1] )
            alternatives.append( [ score, [i-1,j-1] ] )

            score = DP_score[i-1][j] # currently no gap penalty
            alternatives.append( [ score, [i-1,j]] )

            score = DP_score[i][j-1] # currently no gap penalty
            alternatives.append( [ score, [i,j-1]] )

            alternatives.sort()
            alternative = alternatives[ -1 ] # Last element, highest score

            DP_score[i][j] = alternative[ 0 ]
            choice[i][j]   = alternative[ 1 ]



    # Backtrack
    i = totres[0]
    j = totres[1]
    align_seq1 = ''
    align_seq2 = ''
    corresponding_pairs = []
    while ( i > 0 or j > 0 ):

        #print [i,j]

        i_prev = choice[i][j][0]
        j_prev = choice[i][j][1]

        if ( i_prev == i-1  and j_prev == j-1 ):
            # These will be reversed later
            align_seq1 += seq1[ i-1 ]
            align_seq2 += seq2[ j-1 ]
            corresponding_pairs.append( [i, j] )
        elif ( i_prev == i-1 ):
            align_seq1 += seq1[ i-1 ]
            align_seq2 += '-'
        else:
            align_seq1 += '-'
            align_seq2 += seq2[ j-1 ]

        i = i_prev
        j = j_prev

    # Reverse order
    align_seq1_new = ''
    align_seq2_new = ''
    for i in range( len(align_seq1) ): align_seq1_new += align_seq1[ -1 - i ]
    for i in range( len(align_seq2) ): align_seq2_new += align_seq2[ -1 - i ]

    print align_seq1_new
    print align_seq2_new


    if file2.find( '.pdb' ) > 0:
        seqfile = file2.replace('.pdb','.mapping')
    else:
        seqfile = file2 + '.mapping'

    print 'Outputting this alignment to: ', seqfile

    fid = open( seqfile, 'w')
    fid.write( align_seq1_new+'\n' )
    fid.write( align_seq2_new+'\n' )
    fid.close()

    pymol_file = seqfile.replace('.mapping','.pml')
    print 'Display alignment through pymol script: ', pymol_file
    fid2 = open( pymol_file, 'w' )
    fid2.write( 'reinitialize\n')

    command = PYDIR+'/parse_NMR_models.py '+tmp_sup
    system( command )
    file2_sup = file2+'.sup'
    command = 'mv '+tmp_sup.replace('.pdb','_002.pdb') + ' ' + file2_sup
    system( command )

    fid2.write( 'load %s\n' % file1 )
    fid2.write( 'load %s\n' % file2_sup )
    fid2.write( 'hide everything\n' )
    fid2.write( 'show cartoon\n' )
    fid2.write( 'set cartoon_rect_length, 0.5 \n' )
    fid2.write( 'set cartoon_rect_width,  0.2 \n' )
    fid2.write( 'set cartoon_discrete_colors, 1\n' )
    fid2.write( 'color white\n')
    fid2.write( 'bg_color white\n')

    colors = [ 'blue','red','green' ]
    count = 0
    for pair in corresponding_pairs:
        fid2.write( 'color %s, %s and resi %d\n' % ( colors[count], file1.replace('.pdb',''), pair[0])  )
        fid2.write( 'color %s, %s and resi %d\n' % ( colors[count], file2_sup, pair[1])  )
        count += 1
        count = count % len( colors )
    fid2.close()

command = 'rm -rf tmp*'
system( command )
