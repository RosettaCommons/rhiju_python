#!/usr/bin/python

from sys import argv,stdout
import string
from make_tag import make_tag
from parse_options import parse_options

OPTIMAL_LENGTH = parse_options( argv, 'optimal_length', 10 )
LENGTH_STRENGTH = parse_options( argv, 'length_strength', 0.1 )

file = argv[1]
assert( file.count( 'psipred_ss2' ) )

lines = open( file ).readlines()

ss_prob = []
ss = []
for line in lines[2:]:
    cols = string.split( line )
    ss_prob.append( 1.0 - float( cols[3] ) )
    ss.append( cols[2] )

nres = len( ss_prob )

##################################
def plot_stuff( ss_prob, ss, endpoints):

    MAX_BAR_SIZE = 40
    EXTRA_SIZE = 10;

    for i in range( nres ):
        stdout.write( '%3d %s ' % (i+1, ss[i]) )
        bar_size =  int( ss_prob[i]*MAX_BAR_SIZE + 0.5 )
        for j in range( bar_size ):
            stdout.write( 'X' )
        if (i+1) in endpoints:
            for j in range( MAX_BAR_SIZE+EXTRA_SIZE - bar_size ):
                stdout.write( '-' )

        stdout.write( '\n' )


# Pick endpoints by dynamic programming.
MIN_FRAG_LENGTH = 3
MAX_FRAG_LENGTH = 20

# Main loop -- fill in scores recursively
scores = [ 0.0 ]
prev_pos_choice = [ -1 ]

for i in range( 1, MIN_FRAG_LENGTH ):
    scores.append( 99999.99999 )
    prev_pos_choice.append( -1 )

for i in range( MIN_FRAG_LENGTH, nres ):

    local_scores = []
    for j in range( MIN_FRAG_LENGTH, MAX_FRAG_LENGTH+1 ):
        prev_pos = i - j + 1
        if ( prev_pos < 0 ): continue

        score =  scores[ prev_pos ]
        score += ( 1.0 - ss_prob[ i ] )
        score += LENGTH_STRENGTH * abs( j - OPTIMAL_LENGTH )

        local_scores.append( [score, prev_pos ] )

    local_scores.sort()
    scores.append( local_scores[ 0 ][ 0 ] )
    prev_pos_choice.append( local_scores[ 0 ][ 1 ] )

print nres
pos = nres - 1
endpoints = []
while pos > -1:
    endpoints.append( pos+1 )
    pos = prev_pos_choice[ pos ]
endpoints.sort()

plot_stuff( ss_prob, ss, endpoints )

print
print ' -endpoints', make_tag( endpoints )


