#!/usr/bin/python

from sys import argv
from get_disulf import get_disulf

pdb_file = argv[ 1 ]

disulf_pairs = get_disulf( pdb_file )

for i in range( len( disulf_pairs ) ):
    print disulf_pairs[i][0], disulf_pairs[i][1]
