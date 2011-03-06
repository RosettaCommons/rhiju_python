#!/usr/bin/python

from sys import argv
from os.path import exists
from os import system

EXE = '/Applications/RNAVIEW/bin/RNAVIEW'

native_pdb = argv[ 1]
decoy_pdbs  = argv[ 2: ]

def is_base( c ):
    if c in [ '+','W','S','H' ]:
        return 1
    return 0

def get_bp( pdb ):
    outfile =  pdb+'.out'

    if not exists(outfile ):
        system( EXE + ' '+pdb+' > /dev/null 2> /dev/null ' )
    assert( exists( outfile ) )
    lines = open( outfile ).readlines()

    base_pair_list = []
    start_read = 0
    for line in lines:
        if line.count('BEGIN_base-pair')>0 :
            start_read = 1
            continue

        if start_read and len( line ) > 40:
            bp1 = line[ 33 ]
            bp2 = line[ 35 ]

            res1 = int( line[17:19] )
            res2 = int( line[27:29] )

            orientation = line[ 37:41 ]
            info = line[49:-1]

            if not is_base( bp1 ): continue
            if not is_base( bp2 ): continue

            if info.count( 'O1P' ) : continue

            if info.count( '!1H(b_b)' ):  # These cases are ambiguous
                bp1 = 'W'
                bp2 = 'W'

            #print bp1, bp2, res1, res2, orientation, info
            base_pair_list.append( [res1, res2, bp1, bp2 ] )

        if line.count('END_base-pair')>0 :
            break

    #print base_pair_list

    return base_pair_list

#####################################
bp_native = get_bp( native_pdb )

for decoy_pdb in decoy_pdbs:
    bp_decoy = get_bp( decoy_pdb )

    N_BP = 0
    N_BP_IN_DECOY = 0
    for bp in bp_native:
        N_BP += 1
        if bp in bp_decoy:
            N_BP_IN_DECOY += 1
        else:
            continue
            print bp

    f_N_BP = (1.0*N_BP_IN_DECOY) / N_BP

    print decoy_pdb, f_N_BP
