#!/usr/bin/python

from sys import argv,exit,stderr
from os import popen, system
from os.path import basename
import string



def Help():
    print
    print 'Usage: '+argv[0]+' <silent out file 1> < silent file 2> ... <N> '
    print '  Will extract N decoys with lowest score from each silent file.'
    print '  If you want to select based on another column, say 12 (Rg), the'
    print '    last arguments should be -12 <N>  (for lowest Rg) or +12 <N>'
    print '    (for highest Rg).'
    print

    exit()


if len(argv)<2:
    Help()


try:
    NSTRUCT_IN = float(argv[-1])
    del(argv[-1])
except:
    NSTRUCT_IN = 2

scorecol_defined = 0
try:
    scorecol = int(argv[-1])
    del(argv[-1])
    scorecol_defined = 1
except:
    scorecol = -1

REVERSE = ''
if scorecol > 0:
    REVERSE = ' --reverse '

#Another possibility... user supplies -rms or +rms
scorecol_name_defined = 0
if not scorecol_defined:
    scorecol_name = argv[-1]
    if scorecol_name[0] == '-':
        scorecol_name_defined = 1
        scorecol_name = scorecol_name[1:]
        del( argv[-1] )
        REVERSE = ''
    if scorecol_name[0] == '+':
        scorecol_name_defined = 1
        scorecol_name = scorecol_name[1:]
        REVERSE = '-r'
        del( argv[-1] )


infiles = argv[1:]

for infile in infiles:
    tags = []

    firstlines = popen('head -n 3 '+infile).readlines()
    scoretags = string.split( firstlines[1] )
    scoretag=''
    if scorecol_defined:
        scoretag = scoretags[ abs(scorecol) ]

    if scorecol_name_defined:
        scorecol_names = string.split( scorecol_name,',' )
        scorecols = []
        for s in scorecol_names:
            assert( scoretags.count( s ))
            scorecol = scoretags.index( s )
            scorecols.append( scorecol )
        scoretag = scorecol_name
    else:
        scorecols  = [scorecol]

    #   print 'grep SCORE '+infile+' |  sort -k %d -n %s | head -n %d' % (abs(SCORECOL)+1, REVERSE, NSTRUCT+1)

    NSTRUCT = NSTRUCT_IN
    if (NSTRUCT < 1.0 ):
        NUMDECOYS = int( string.split(popen('grep SCORE '+infile+' | wc').readlines()[0])[0] ) - 1
        NSTRUCT = round( NSTRUCT_IN * NUMDECOYS )

    #lines = popen('grep SCORE '+infile+' | grep -v NATIVE | grep -v RT | sort -k %d -n %s | head -n %d' % (abs(SCORECOL)+1, REVERSE, NSTRUCT+1) ).readlines()

    # Make the list of decoys to extract
    command = 'grep SCORE '+infile+' | grep -v NATIVE'
    lines = popen( command ).readlines()
    #stderr.write( '%s %d\n' % ( command, len( lines ) ) )

    score_plus_lines = []
    for line in lines:
        cols = string.split( line )
        score = 0.0
        try:
            for scorecol in scorecols: score += float( cols[ abs(scorecol) ] )
        except:
            continue
        if REVERSE: score *= -1
        score_plus_lines.append( ( score, line ))

    score_plus_lines.sort()

    #stderr.write( '%d' % len(score_plus_lines[:NSTRUCT]) )

    #stderr.write( '%d %d %d \n' % (NSTRUCT,int( NSTRUCT),len(lines) ) )

    lines = map( lambda x:x[-1], score_plus_lines[:int(NSTRUCT)] )


    #lines = popen('grep SCORE: '+infile+' | grep -v NATIVE | grep -v rms | sort -k %d -n %s | head -n %d' % (abs(SCORECOL)+1, REVERSE, NSTRUCT+1) ).readlines()

    templist_name = 'temp.%s.list'% basename(infile)

    fid = open(templist_name,'w')
    count = 0
    for line in lines:
        cols = string.split(line)
        tag = cols[-1]
        if tag.find('desc') < 0:
            fid.write(tag+'\n')
            tags.append(tag)
            count = count+1
        if count >= NSTRUCT:
            break
    outfilename = infile
    fid.close()


    if (firstlines[2][:6] == 'REMARK' ):
        command = 'head -n 3 '+infile
        system(command)
    else:
        command = 'head -n 2 '+infile
        system(command)


    count = 1
    fid = open( infile )
    line = fid.readline()

    writeout = 0
    while line:
        cols = string.split(line)
	if (len(cols)>1 and cols[0]=='SCORE:'):
	    if tags.count(cols[-1]) > 0:
		writeout = 1
	    else:
		writeout = 0
        if writeout:
            print line[:-1]
        line = fid.readline()

    command = 'rm '+templist_name
    #    print(command)
    system(command)


