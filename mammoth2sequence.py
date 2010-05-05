#!/usr/bin/python

from sys import argv
import string

file = argv[1]
lines = open( file ).readlines()

seq1 = ''
seq2 = ''

for i in  range( len(lines)-4):
    if len( lines[i] ) < 10: continue
    if len( lines[i+1] ) < 10: continue
    if len( lines[i+3] ) < 10: continue
    if len( lines[i+4] ) < 10: continue
    if ( not lines[i][:10]  == 'Prediction' ): continue
    if ( not lines[i+1][:10]  == 'Prediction' ): continue
    if ( not lines[i+3][:10]  == 'Experiment' ): continue
    if ( not lines[i+4][:10]  == 'Experiment' ): continue

    seq_chunks = string.split( lines[i] )[1:]
    seq1 += string.join( seq_chunks, '' )

    seq_chunks = string.split( lines[i+4] )[1:]
    seq2 += string.join( seq_chunks, '' )

seq1 = seq1.replace('.','-')
seq2 = seq2.replace('.','-')

# This is really silly, but mammoth does not align last residue in sequence.

