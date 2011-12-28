#!/usr/bin/python

from math import sqrt

def get_positions( lines, atom_names_in ):

    count = 0
    positions = []
    which_res = []
    atom_names = []
    oldresnum = ''
    for line in lines:

        if (len(line)>54 and  line[0:4] == 'ATOM' ):

            resnum = line[23:26]
            if not resnum == oldresnum:
                count = count + 1

            atom_name = line[12:16]
            position = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

            if atom_name in atom_names_in:
                positions.append( position )
                which_res.append( count )
                atom_names.append( atom_name )

            oldresnum = resnum

    return ( positions, which_res, atom_names, count )

def get_dist( pos1, pos2 ):
    dist2 = 0.0
    for k in range(3):
        dsub = ( pos1[k] - pos2[k] )
        dist2 += dsub*dsub

    return sqrt( dist2 )


def get_disulf( pdb_file ):

    lines = open( pdb_file ).readlines()
    ( positions1, which_res1, atom_names1, totres) = get_positions( lines, [ " SG " ] )
    ( positions2, which_res2, atom_names2, totres) = get_positions( lines, [ " SG " ] )
    SEQ_SEP_CUTOFF = 1
    DIST_CUT = 2.9

    disulf_pairs = []

    for i in range( len( which_res1) ):

        res1 = which_res1[ i ]
        atom1 = atom_names1[ i ]

        for j in range( i, len( which_res2) ):

            res2 = which_res2[ j ]
            atom2 = atom_names2[ j ]

            if ( abs( res1 - res2 ) <= SEQ_SEP_CUTOFF ): continue

            dist = get_dist( positions1[ i ], positions2[ j ] )

            if ( dist < DIST_CUT ):
                disulf_pairs.append( [res1, res2] )

    return disulf_pairs
