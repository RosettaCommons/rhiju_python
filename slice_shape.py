#!/usr/bin/python

import sys
import string

infile = sys.argv[1]
startpos = int( sys.argv[2] )
endpos   = int( sys.argv[3] )
lines = open(  infile ).readlines()

newfile = 'region_%d_%d_%s' % (startpos,endpos,infile)
fid = open( newfile, 'w' )

for line in lines:
    cols = string.split(line)
    pos = int( cols[0] )
    if pos >= startpos and pos <= endpos:
        fid.write( '%d %s\n' % (pos-startpos+1, cols[1] ) )

fid.close()

