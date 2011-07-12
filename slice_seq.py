#!/usr/bin/python

import sys

infile = sys.argv[1]
startpos = int( sys.argv[2] )
endpos   = int( sys.argv[3] )
lines = open(  infile ).readlines()

newfile = 'region_%d_%d_%s' % (startpos,endpos,infile)
fid = open( newfile, 'w' )

fid.write( lines[0] )
fid.write( lines[1] )
sequence = lines[2][:-2]
fid.write( sequence[startpos-1:endpos] )
fid.write( '1\n')


