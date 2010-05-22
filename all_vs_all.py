#!/usr/bin/python

from sys import argv,stdout,stderr
from os import system
import string
from parse_options import parse_options

rmsd_threshold = parse_options( argv, 'R', 8.0);
subset_res = parse_options( argv,'subset',[-1])

infilelist = argv[1]

lines = open(infilelist,'r').readlines()
lines = [x[:-1] for x in lines]

subset_tag = ''
if len( subset_res ) > 0:
    subset_tag = ' -subset'
    for m in subset_res: subset_tag+= ' %d' % m

fit_threshold_save = {}
maxsub_save = {}
for line in lines:
    maxsub_save[line] = {}
    fit_threshold_save[line] = {}


for i in range(len(lines)):
    line1 = lines[i]
    for j in range(i, len(lines)):
        line2 = lines[j]

        command = '~rhiju/python/superimpose.py %s %s  %s -R %d > q 2> blah.err' % (line1,line2,subset_tag,rmsd_threshold)
#        print(command)
        system(command)

        superimposeline = open( 'blah.err', 'r').readlines()[-1]
 #       print superimposeline

        fit_threshold = float( string.split(superimposeline)[3] )
        maxsub = int( string.split(superimposeline)[5] )
        if fit_threshold > rmsd_threshold:
            maxsub = 0 # Failure!
            fit_threshold = -1

        stderr.write('%s %s %4.2f %d\n' % (line1, line2, fit_threshold, maxsub))

        maxsub_save[line1][line2] = maxsub
        maxsub_save[line2][line1] = maxsub

        fit_threshold_save[line1][line2] = fit_threshold
        fit_threshold_save[line2][line1] = fit_threshold

print
print


maxlen = max( [len(x) for x in lines] )
blanks = ' '*200
for line1 in lines:
    print '%s' % line1+blanks[len(line1):maxlen],
    for line2 in lines:
        print '%03d' % maxsub_save[line1][line2],
    print

print
for line1 in lines:
    print '%s' % line1+blanks[len(line1):maxlen],
    for line2 in lines:
        print '%4.2f' % fit_threshold_save[line1][line2],
    print

command = 'rm -rf maxsub*pdb blah*pdb'
system(command)
