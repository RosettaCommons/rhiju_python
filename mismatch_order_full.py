#!/usr/bin/python

from sys import argv,stdout,exit
import string

def Help():
    print
    print argv[0]+ ' <name> <offset> [<library_num>] <sequence> <pos1> <pos2> ... '
    print
    print ' A little script to generate oligo sequences and 96-well plate assignments'
    print '  for mutate/map experiments.'
    print
    print ' name = arbitrary tag for sequence'
    print ' offset = position in sequence that is numbered "1"'
    print ' sequence = string of A,G,C, T (or U), including T7 promoter'
    print ' library_num = 1, 2, or 3, for the different mutation possibilities.'
    print
    print ' pos1, pos2, ... = positions at which to introduce mismatch'
    print
    exit()

if ( len( argv ) < 6 ): Help()
name = argv[1]
offset = int( argv[2] )

sequence = argv[3]
positions = [ -999 ] + map( lambda x:int(x), argv[4:] )

make_DNA = { 'A':'A', 'T':'T', 'U':'T', 'C':'C', 'G':'G' }

complement = { 'A':'T', 'T':'A', 'U':'A', 'C':'G', 'G':'C' }
alterna1   = { 'A':'C', 'T':'C', 'U':'C', 'C':'A', 'G':'A' }
alterna2   = { 'A':'G', 'T':'G', 'U':'G', 'C':'T', 'G':'T' }
mismatch = [ complement, alterna1, alterna2 ]

numchar = len( sequence )
scanchar = numchar

MAXCHAR = 60
if (numchar > MAXCHAR ):
    #print 'Truncating down to', MAXCHAR,' !!!!!'
    scanchar = MAXCHAR

seq1s = []
seq2s = []
new_names = []

for library_num in range( 1, len( mismatch )+1):
    for position in positions:

        sequence_new = ''
        for j in range( len( sequence ) ):
            if j == offset + position - 1 and not position == -999:
                seqchar = sequence[j]
                seqchar_new = mismatch[library_num-1][ seqchar ]
                sequence_new += seqchar_new
            else:
                sequence_new += sequence[j]

        seq1 = ''
        seq2 = ''

        for i in range( scanchar ):
            seq1 += make_DNA[ sequence_new[i] ]
            seq2 += complement[ sequence_new[numchar - 1 - i] ]

        if position == -999:
            seq1_WT = seq1
            seq2_WT = seq2

        seq1s.append( seq1 )
        seq2s.append( seq2 )

        if position == -999:
            new_name = '%s-WT' % name
        else:
            new_name = '%s-%s%02d%s' % (name, seqchar, position, seqchar_new )

        new_names.append( new_name )

tot_constructs = len( new_names )
ncols = ( (tot_constructs-1) / 8 ) + 1

letters = 'ABCDEFGH'

for j in range( len( new_names ) ):
    well = '%s%d' % ( letters[ j % 8 ],  (j/8)+1 )
    if not( seq1s[j]  == seq1_WT ) or new_names[j][-2:] == 'WT':
        stdout.write('%s;%s;%s;\n' % (well, new_names[j]+'-F', seq1s[j]) )
print
print

for j in range( len( new_names ) ):
    well = '%s%d' % ( letters[ j % 8 ],  (j/8)+1 )
    if not( seq2s[j]  == seq2_WT ) or new_names[j][-2:] == 'WT':
        stdout.write('%s;%s;%s;\n' % (well, new_names[j]+'-R', seq2s[j] ) )


############

keys_file = name+'_keys.txt'
print 'Putting keys in ', keys_file
fid = open( keys_file , 'w' )
for j in range( len( new_names ) ):
    fid.write( new_names[j].replace( name+'-','')  + '\n' )
fid.close()


############
make_RNA = { 'A':'A', 'U':'U', 'T':'U', 'C':'C', 'G':'G' }

sequence_truncate = sequence;
if sequence_truncate[:20] == 'TTCTAATACGACTCACTATA':
    sequence_truncate = sequence_truncate[20:]

sequence_file = name+'_sequence.txt'
print 'Putting RNA sequence in ', sequence_file
fid = open( sequence_file, 'w' )
for j in range( len( sequence_truncate) ):
    fid.write( make_RNA[ sequence_truncate[j] ] )
fid.write( '\n')
fid.close()

