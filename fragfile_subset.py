#!/usr/bin/python

from os import popen,system
import sys
import string

fragfilefile = sys.argv[1]
outfile = sys.argv[2]
subset_residues = map( lambda x:int(x),  sys.argv[3:] )
subset_residues.sort()

gzipped = 0

if outfile[-2:] == 'gz':
    outfile = outfile[:-3]
    gzipped = 1
outid = open(outfile,'w')

if fragfilefile[-2:] == 'gz':
    lines = popen('gzcat '+fragfilefile).readlines()
    fragfilefile = fragfilefile[:-3]
else:
    lines = open(fragfilefile).readlines()

startoutput = 0
fragmentsize = int( fragfilefile[-14:-12])
for line in lines:
    cols = string.split(line)
    if cols.count('position:'):
        i = int(cols[1])
        if i in subset_residues:
            startoutput = 1
        else:
            startoutput = 0
            continue
        newresnum = '%4d' % ( subset_residues.index(i)+1 )
        line = line[:19] + newresnum + line[23:]

    if startoutput:
        outid.write(line)

outid.close()

if gzipped:
    command = 'gzip -f '+outfile
    system(command)
