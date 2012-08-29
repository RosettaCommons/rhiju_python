#!/usr/bin/python

from sys import argv, exit, stdout
import string
from os import system, popen
from os.path import basename,dirname, exists



alignfile = argv[1]
startnum = int(argv[2])
endnum   = int(argv[3])
newprefix = argv[4]

indirs = argv[5:]

command = '/users/rhiju/python/map_sequence_numbers.py %s %d %d %s' % (alignfile, startnum,endnum,newprefix)
print(command)
lines = popen(command).readlines()

print len(lines)
print len(indirs)

assert( len(lines) == len(indirs))

for i in range(len(lines)):
    line  = lines[i]
    startnum = int( string.split(line)[1])
    endnum   = int( string.split(line)[2])
    indir = indirs[i]

    command = '/users/rhiju/python/truncate_rosetta_files.py %s %s %d %d %s' % (indir,
                                          newprefix+basename(indir),startnum,endnum,newprefix)
    print(command)
    system(command)
