#!/usr/bin/python

import string
from sys import argv,stdout,stderr
from os import system

############################
subset_residues = []
if argv.count('-segments'):
    use_subset = 1
    pos = argv.index('-segments')
    del argv[pos]

    segment_residues = []
    stderr.write( 'SLICE using a subset of residues: '  )
    goodint = 1
    while goodint:
        try:
            segment_residue = int(argv[pos])
            segment_residues.append( segment_residue )
            del argv[pos]
            stderr.write('%d ' % segment_residue )
        except:
            goodint = 0

    stderr.write( '\n'  )

    pdbfiles = argv[1:-1]

    prefix = argv[-1]
    startseq = 1
    endseq = 10000000000

    segment_ends = []
    for i in range( len(segment_residues)/2):
        segment_ends.append( segment_residues[2*i+1] )
        for j in range( segment_residues[2*i], segment_residues[2*i+1]+1 ):
            subset_residues.append( j )

    print subset_residues

prefix = argv[-1]

for file in argv[1:-1]:

    if len(file ) > 6 and file[-6:] == '.fasta':
        lines = open( file ).readlines()
        fid = open( prefix + file, 'w' )
        name_line = lines[0]
        fid.write( '>' + prefix+ name_line[1:] )

        seq_line = lines[1]
        for k in subset_residues:
            fid.write( seq_line[k-1] )
        fid.write('\n')
        fid.close()
    elif len(file ) > 4 and file[ -4:] == '.pdb':
        command = 'pdbslice.py '+file+' -segments '
        for k in segment_residues:
            command +=' %d' % k
        command += ' '+prefix
        print( command )
        system( command )
    elif len(file ) > 4 and file[-4:] == '.cst':
        lines = open( file ).readlines()
        fid = open( prefix + file, 'w' )
        for line in lines:
            if len(line) > 1 and line[0] == '[':
                fid.write( line )
                continue

            cols = string.split( line )
            pos1 = int( cols[1] );
            pos2 = int( cols[3] );
            if ( pos1 in subset_residues) and \
                    ( pos2 in subset_residues):
                pos1_map = subset_residues.index( pos1 )+1
                pos2_map = subset_residues.index( pos2 )+1
                fid.write('%s %d %s %d %s\n' %  \
                             (cols[0],pos1_map,cols[2],pos2_map,string.join( cols[4:] ) ) )
        fid.close()
    elif len(file ) > 5 and file[-5:] == '.data':
        lines = open( file ).readlines()
        pos_map = []
        for line in lines:
            if len(line) > 15 and line[:15] == 'BACKBONE_BURIAL':
                cols = string.split( line )
                for k in cols[1:]:
                    if int(k) in subset_residues:
                        print k, len( subset_residues ), subset_residues.index( int(k) ) + 1
                        pos_map.append( subset_residues.index( int(k) ) + 1 )
        if len( pos_map ) > 0:
            fid = open( prefix + file, 'w' )
            fid.write( 'BACKBONE_BURIAL ' )
            for k in pos_map:
                fid.write(' %d' % k );
            fid.write('\n')
            fid.close()
    elif len(file ) > 7 and file[-7:] == '.params':
        lines = open( file ).readlines()
        fid = open( prefix + file, 'w' )
        for line in lines:
            cutpoints = []

            #for k in segment_ends[:-1]:
            #    cutpoints.append( subset_residues.index( k ) + 1 )

            if len(line) > 14 and line[:13] == 'CUTPOINT_OPEN':
                cols = string.split( line )
                for k in cols[1:]:
                    if ( int( k ) in subset_residues ):
                        cutpoints.append( subset_residues.index( int(k) ) + 1 )
                if len( cutpoints ) > 0 :
                    fid.write( 'CUTPOINT_OPEN ' )
                    for k in cutpoints: fid.write( ' %d' % k )
                    fid.write( '\n' )
            if len(line) > 14 and line[:15] == 'CUTPOINT_CLOSED':
                cols = string.split( line )
                for k in cols[1:]:
                    if ( int( k ) in subset_residues ):
                        cutpoints.append( subset_residues.index( int(k) ) + 1 )
                if len( cutpoints ) > 0 :
                    fid.write( 'CUTPOINT_CLOSED ' )
                    for k in cutpoints: fid.write( ' %d' % k )
                    fid.write( '\n' )
            elif len( line ) > 8 and line[:8] == 'OBLIGATE':
                cols = string.split( line )
                pos1 = int( cols[2] )
                pos2 = int( cols[3] )
                if ( pos1 in subset_residues) and \
                        ( pos2 in subset_residues):
                    pos1_map = subset_residues.index( pos1 )+1
                    pos2_map = subset_residues.index( pos2 )+1
                    fid.write('%s %d %d %s %s %s\n' %  \
                                  ('OBLIGATE  PAIR',pos1_map,pos2_map,cols[4],cols[5],cols[6] ) )
            elif len( line ) > 4 and line[:4] == 'STEM':
                cols_all = string.split( line )
                pairs = []
                for i in range( len( cols_all ) / 6 ):
                    cols = ['blah']
                    for k in range(6): cols.append( cols_all[ 1+6*i+k ] )
                    pos1 = int( cols[2] )
                    pos2 = int( cols[3] )
                    if ( pos1 in subset_residues) and \
                            ( pos2 in subset_residues):
                        pos1_map = subset_residues.index( pos1 )+1
                        pos2_map = subset_residues.index( pos2 )+1
                        pairs.append(' %s %d %d %s %s %s ' %  \
                                      ('PAIR',pos1_map,pos2_map,cols[4],cols[5],cols[6] ))
                if len( pairs ) > 0:
                    fid.write( 'STEM ' )
                    for k in pairs: fid.write( k )
                    fid.write( '\n' )
        fid.close()



