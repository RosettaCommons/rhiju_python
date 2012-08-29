#!/usr/bin/python

from sys import argv
from os import popen,system
from os.path import exists
import string
from glob import glob

submitfile = argv[1]
batchid = 0
if len(argv)>2:
    batchid = int( argv[2])

data = open(submitfile,'r')

if not exists('filters'): system('mkdir filters')


keeptesting = 1
count = 1
dirnames = []
line = data.readline()
while line and keeptesting:
    while line and line[:4] != 'name': line = data.readline()
    if not line: break

    WUname = string.split(line,'=')[-1][1:-1]

    #Mirror the results from ralph
    dirname = '/work/boinc/results_ralph/'+WUname[:8]
    if dirname not in dirnames:
        dirnames.append(dirname)
        command = '/users/rhiju/python/update_results.py '+dirname
        print(command)
        system(command)

    relax_score_filter1 = []
    relax_score_filter2 = {}

    # Find the scorefile -- assume its the first scorefile available with the WU name.
    if batchid > 0:
        scorefiles = glob('%s/%s*%d*sc.bz2' % (dirname,WUname,batchid))
    else:
        scorefiles = glob('%s/%s*sc.bz2' % (dirname,WUname))
    assert( len(scorefiles)>0)
    scorefile = scorefiles[-1]

    # Figure out scorefilters.
    lines = popen('bzcat '+scorefile+ '| grep SCORE ').readlines()

    line = lines[0]
    cols = string.split(line)
    if not cols.count('rlxfilt1'):
        print 'Hey! Outfiles must have columns with rlxfilt1 and rlxfilt2'

    index1 = cols.index('rlxfilt1')
    index2 = cols.index('rlxfilt2')

    for line in lines[1:]:
        cols = string.split(line)
        try:
            relax_score_filter1.append( float(cols[index1]))
            relax_score_filter2[ float(cols[index1]) ] = float(cols[index2])
        except:
            print 'BAD SCORE? ', cols[index1], cols[index2]
            continue

    relax_score_filter1.sort()
    numdecoys = len( relax_score_filter1 )

    set_score_filter1 = relax_score_filter1[ int(numdecoys/ 2) ]


    relax_score_filter2_cull = []
    for x in relax_score_filter1[ :int(numdecoys/2)]:
        if x in relax_score_filter2.keys():
            relax_score_filter2_cull.append( relax_score_filter2[x] )

    relax_score_filter2_cull.sort()

    set_score_filter2 = relax_score_filter2_cull[ int(numdecoys/ 4) ]


    #Save a filter file.
    while line and line[:9] != 'arguments': line = data.readline()
    if not line: break

    arguments = string.split(line,'=')[-1][1:-1]
    cols = string.split(arguments)
    fivelettercode = cols[1]+cols[2]

    pos = cols.index('-protein_name_prefix')
    prefix = cols[pos+1]


    fid = open( 'filters/'+prefix+fivelettercode+'filters.txt','w')
    fid.write('%6.2f %6.2f\n' % (set_score_filter1, set_score_filter2) );
    print set_score_filter1, set_score_filter2
    fid.close()

