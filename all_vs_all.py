#!/usr/bin/python

from sys import argv,stdout,stderr
from os import system
import string
from parse_options import parse_options

rmsd_threshold = parse_options( argv, 'R', 8.0);
subset_res = parse_options( argv,'subset',[-1])

EXE = 'superimpose.py'

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

        if ( rmsd_threshold == 8.0 ):
            command = EXE+' %s %s  %s > q 2> blah.err' % (line1,line2,subset_tag)
        else:
            command = EXE+' %s %s  %s -R %6.2f > q 2> blah.err' % (line1,line2,subset_tag,rmsd_threshold)
        print(command)
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

print


mean_maxsubs = []
for line1 in lines:
    mean_maxsub = 0.0
    for line2 in lines:
        if line1 == line2: continue
        mean_maxsub += maxsub_save[ line1 ][ line2 ]
    mean_maxsub /= (len( lines ) - 1 )
    mean_maxsubs.append( [mean_maxsub, line1] )

mean_maxsubs.sort()
mean_maxsubs.reverse()
for i in range( len(mean_maxsubs) ):
    line1 = mean_maxsubs[i][1]
    mean_maxsub = mean_maxsubs[i][0]
    print '%s' % line1+blanks[len(line1):maxlen],
    print '%4.2f' % mean_maxsub
print


mean_rmsds = []
for line1 in lines:
    mean_rmsd = 0.0
    for line2 in lines:
        mean_rmsd += fit_threshold_save[ line1 ][ line2 ]
    mean_rmsd /= (len( lines ) - 1 )
    mean_rmsds.append( [mean_rmsd, line1] )

mean_rmsds.sort()
for i in range( len(mean_rmsds) ):
    line1 = mean_rmsds[i][1]
    mean_rmsd = mean_rmsds[i][0]
    print '%s' % line1+blanks[len(line1):maxlen],
    print '%4.2f' % mean_rmsd

command = 'rm -rf maxsub*pdb blah*pdb'
system(command)
