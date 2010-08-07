#!/usr/bin/python

import string
from sys import argv
from os.path import basename,exists
from get_sequence import get_sequence
from amino_acids import short_to_long
from blast import NBAlign


native_pdb = argv[1]
ref_pdb = argv[2]

native_sequence = get_sequence( native_pdb )
ref_sequence = get_sequence( ref_pdb )

al = NBAlign( ref_sequence, native_sequence )

new_file =  native_pdb.replace('.pdb','_MATCH.pdb' )
outid = open( new_file, 'w' )

lines = open( native_pdb ).readlines()
oldresnum = ''

all_res_lines = []
res_lines = []
for line in lines:
    line_edit = line
    if line[0:3] == 'TER':
        continue

    if line_edit[0:4] == 'ATOM' or line_edit[0:6] == 'HETATM':

        if not (line[16]==' ' or line[16]=='A'): continue

        resnum = line_edit[23:26]
        if not resnum == oldresnum:
            oldresnum = resnum
            if len( res_lines) > 0: all_res_lines.append( res_lines )
            res_lines = []

        res_lines.append( line )

if len( res_lines) > 0: all_res_lines.append( res_lines )

atomnum = 0
for i in range( len( ref_sequence ) ):
    if i in al.keys():
        #print i, al[ i ]
        for line in all_res_lines[ al[ i ] ]:
            atomnum += 1
            newnum = '%4d' %  (i+1)
            line_edit = '%s%5d%s%s%s' % (line[0:6],atomnum,line[11:22], newnum, line[26:] )
            outid.write( line_edit )
    else:
        atomnum += 1
        outid.write( 'ATOM   %4d  CA  %3s   %3d     -17.149  26.321  14.330  1.00 96.92           C\n' % ( atomnum, short_to_long[  ref_sequence[ i ] ], i+1 ) )
        outid.write( 'ATOM   %4d  C   %3s   %3d     -17.149  26.321  14.330  1.00 96.92           C\n' % ( atomnum, short_to_long[  ref_sequence[ i ] ], i+1 ) )
        outid.write( 'ATOM   %4d  N   %3s   %3d     -17.149  26.321  14.330  1.00 96.92           C\n' % ( atomnum, short_to_long[  ref_sequence[ i ] ], i+1 ) )

outid.close()




