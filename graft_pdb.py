#!/usr/bin/python

from read_pdb import read_pdb
from sys import argv

# Note, this assumes that the pdb's are already aligned!

main_pdb = argv[1]
scratch_pdb = argv[2]

[ coords_main, lines_main ] = read_pdb( main_pdb )
[ coords_scratch, lines_scratch ] = read_pdb( scratch_pdb )

# remove atoms in main_pdb that are being replaced by scratch_pdb:
for chain in lines_scratch.keys():
    if chain not in lines_main.keys(): continue

    for resi in lines_scratch[ chain ]:
        if resi not in lines_main[ chain ].keys(): continue

        for atom in lines_scratch[ chain ][ resi ]:
            if atom not in lines_main[ chain ][ resi ].keys(): continue

            del( lines_main[ chain ][ resi ][ atom ] )


for chain in lines_scratch.keys():
    if not chain in lines_main.keys(): lines_main[ chain ] = {}

    for resi in lines_scratch[ chain ]:
        if not resi in lines_main[ chain ].keys(): lines_main[ chain ][ resi ] = {}

        for atom in lines_scratch[ chain ][ resi ]:
            lines_main[ chain ][ resi ][ atom ] = lines_scratch[ chain ][ resi ][ atom ]

chains = lines_main.keys()
chains.sort()
for chain in chains:

    residues = lines_main[ chain ].keys()
    residues.sort()

    for resi in residues:
        for atom in lines_main[ chain ][ resi ]:
            print lines_main[ chain ][ resi ][ atom ]

