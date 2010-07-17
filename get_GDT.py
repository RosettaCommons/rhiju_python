#!/usr/bin/python

from os import system
from os.path import basename
from subprocess import *
import string
from sys import argv
from get_sequence import get_sequence
from blast import NBAlign
from make_tag import make_tag

native = argv[1]
predictions = argv[2:]

#sequence matcher...
native_sequence = get_sequence( native )
tot_length = len( native_sequence )

for prediction in predictions:

    prediction_sequence = get_sequence( prediction )
    pdb_to_superimpose = 'tmp_'+basename(prediction)

    if native_sequence == prediction_sequence:
        system( 'cp %s %s' % ( prediction, pdb_to_superimpose ) )
    else:
        al = NBAlign( prediction_sequence, native_sequence )

        slice_res = []
        for i in range( len( prediction_sequence ) ):
            if i in al.keys(): slice_res.append( i+1 )

        command = 'pdbslice.py %s -subset %s tmp_' % ( prediction, make_tag( slice_res ) )
        system( command )

    if not ( get_sequence( pdb_to_superimpose ) == native_sequence ): continue

    maxsubs = []
    #cluster_threshold = [  2,   3, 5, 8, 16 ]
    cluster_threshold = [  0.5, 1, 2, 4,  8 ]
    rmsd_threshold     = [ 0.5, 1, 2, 4,  8 ]
    assert( len( cluster_threshold) == len( rmsd_threshold ) )

    for i in range( len( cluster_threshold) ):
        command = 'superimpose.py %s  %s  -R %5.2f -per_res > q 2> per_res.txt' % \
            (native, pdb_to_superimpose, cluster_threshold[ i ] )
        #command = 'superimpose.py %s  %s  -per_res > q 2> per_res.txt' % \
        #    (native, pdb_to_superimpose )
        system( command )

        lines = open( 'per_res.txt' ).readlines()
        nali = 0
        #print lines[1]
        for line in lines:
            cols = string.split( line )
            if len(cols) == 2 and float( cols[1] ) <= rmsd_threshold[ i ]:
                nali += 1
        maxsubs.append( nali )
    #print maxsubs

    system( 'rm -f '+pdb_to_superimpose )

    gdt_ha = (0.25 *(maxsubs[0]+maxsubs[1]+maxsubs[2]+maxsubs[3]))/tot_length
    gdt_ts = (0.25 *(maxsubs[1]+maxsubs[2]+maxsubs[3]+maxsubs[4]))/tot_length

    print 'GDT_HA: %8.3f  GDT_TS: %8.3f ' % (gdt_ha, gdt_ts ),
    print ' MM0.5:%8.3f  MM1:%8.3f'  % ( 1.0*maxsubs[0]/tot_length, 1.0*maxsubs[1]/tot_length),
    print '%s' % ( prediction)


