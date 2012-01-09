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
score_plus_lines = []
firstlines = []
IS_OUTFILE = 1

for infile in infiles:

    if len( firstlines ) == 0:
        firstlines = popen('head -n 3 '+infile).readlines()
        scoretags = string.split( firstlines[1] )
        if firstlines[0].find( "SEQUENCE" ) < 0  and   firstlines[0].find("SCORE:") >= 0:
            IS_OUTFILE = 0
            scoretags = string.split( firstlines[0] )

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

    # Make the list of decoys to extract
    command = 'grep SCORE '+infile+' | grep -v NATIVE'
    lines = popen( command ).readlines()

    for line in lines:
        cols = string.split( line )
        score = 0.0
        try:
            for scorecol in scorecols: score += float( cols[ abs(scorecol) ] )
        except:
            continue
        if REVERSE: score *= -1
        score_plus_lines.append( ( score, line, infile ))

score_plus_lines.sort()

NSTRUCT = NSTRUCT_IN
if (NSTRUCT < 1.0 ):
    NUMDECOYS = len( score_plus_lines )
    NSTRUCT = round( NSTRUCT_IN * NUMDECOYS )
else:
    NSTRUCT = int( NSTRUCT )

tags_for_infile = {}
for infile in infiles: tags_for_infile[ infile ] = []

for score_plus_line in score_plus_lines[:NSTRUCT]:
    line = score_plus_line[1]
    cols = string.split(line)
    tag = cols[-1]

    infile = score_plus_line[2]
    tags_for_infile[ infile ].append( tag )

if not IS_OUTFILE:
    command = 'head -n 1 '+infile
    system(command)
elif (firstlines[2][:6] == 'REMARK' ):
    command = 'head -n 3 '+infile
    system(command)
else:
    command = 'head -n 2 '+infile
    system(command)

# following basically stolen from cat_outfiles.py
n = -1
for infile in infiles:

    n += 1

    ok_tags = tags_for_infile[ infile ]
    if len( ok_tags ) == 0: continue

    data = open(infile,'r')
    line = data.readline() # Skip first two lines
    line = data.readline()

    writeout = 0
    while line:
        line = data.readline()[:-1]
        if len( line ) < 2: continue
        if line[:9] == 'SEQUENCE:': continue # Should not be any more sequence lines!
        cols = string.split( line )
        tag = cols[-1]

        if cols[0][:5] == "SCORE":
            if ok_tags.count( tag ) > 0: writeout = 1
            else: writeout = 0

        if not writeout: continue

        if (n>0):
            tagcols = string.split( tag, '_')
            try:
                tagnum = int( tagcols[-1] )
                tagcols[-1] = '%04d_%06d' %  ( n, tagnum )
                newtag = string.join( tagcols,'_')
                description_index = line.find( tag )
                line = line[:description_index] + newtag
            except:
                continue

        print line

    data.close()
