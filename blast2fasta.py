#!/usr/bin/python

from sys import argv

blast_file = argv[1]

lines = open( blast_file ).readlines()

in_titles = False
finished_titles = False
titles = []
block_count = 0
seq_count = 0
seqs = []
tags = []
for i in range( len( lines ) ):
    line = lines[i]
    if len( line ) > 9 and line[:9] == 'Sequences':
        in_titles = True
        tags.append( 'Query' )
        titles.append( 'Query' )
        continue
    if len( line ) == 1 and len( titles ) > 1:
        in_titles = False
        if not finished_titles:
            finished_titles = True
        else:
            block_count += 1
            seq_count = 0
        continue
    if len( line ) == 1: continue
    if line[:10] ==  'ALIGNMENTS': continue
    if in_titles:
        tags.append( line.split( '|' )[1].split( '.' )[0] )
        titles.append( line[:65] )
        continue

    if len( line ) > 12 and line[:11] == '  Database:': break

    if finished_titles:
        pos = 22
        while pos < len(line ) and line[pos] not in '0123456789': pos += 1
        pos -=2
        seq_tract = line[22:pos].replace( ' ','-' )
        seq_count += 1
        tag = line.split()[0]
        while ( tag != tags[seq_count-1] ):
            dummy_seq_tract = '-'*60
            if len( seqs ) < seq_count :  seqs.append( '' )
            seqs[ seq_count-1 ] += dummy_seq_tract
            seq_count += 1

        if len( seqs ) < seq_count :  seqs.append( '' )
        seqs[ seq_count-1 ] += seq_tract



for i in range( len( titles ) ):
    print '>',titles[i].split( '|' )[-1].replace( ' ','_' )
    print seqs[i]
    print


