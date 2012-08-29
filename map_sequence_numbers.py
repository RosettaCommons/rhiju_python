#!/usr/bin/python

from sys import argv, exit, stdout
import string
from os.path import dirname, basename, exists,abspath



alignfile = argv[1]
startnum = int(argv[2])
endnum   = int(argv[3])

newprefix='temp_'
if len(argv)>4:
    newprefix = argv[4]

alignfile = abspath(alignfile)

lines = open( alignfile, 'r').readlines()

fid = open( dirname(alignfile)+'/'+newprefix+basename(alignfile),'w')

def figureoutnewnum(targetsequence,sequence,startnum):
    targetseqpos = 0
    seqpos = 0
    for k in range( len(targetsequence)):
        alignpos = k+1
        if targetsequence[k] != '-':
            targetseqpos += 1
        if sequence[k] != '-':
            seqpos += 1
        if (targetseqpos >= startnum):
            break

    if seqpos == 0:
        seqpos = 1
    return (seqpos, alignpos)

line = lines[0]
targetsequence = string.split(line)[1]

for i in range( len(lines)):
    line = lines[i]
    sequence = string.split(line)[1]
    (startnum_map, startnum_in_alignfile) = figureoutnewnum(targetsequence,sequence,startnum)
    (endnum_map  , endnum_in_alignfile)   = figureoutnewnum(targetsequence,sequence,endnum)

    tag = string.split(line)[2]
    print tag,startnum_map,endnum_map
    fid.write('ALIGN %s %s\n' % (sequence[startnum_in_alignfile-1 : endnum_in_alignfile], newprefix+tag))

fid.close()
