#!/usr/bin/python

from sys import stdout,argv
from os import system

def extractchain(actualpdbname, out, chain_to_extract):

    if actualpdbname[-3:] =='.gz':
        lines = popen( 'zcat '+actualpdbname).readlines()
    else:
        lines = open(actualpdbname,'r').readlines()

#    out = open(actualpdbname_chain_to_extract,'w')
    for i in range( len(lines)):
        line = lines[i]
        if line.count('ATOM') and (line[21:22] == chain_to_extract ):
            line = line[0:21]+chain_to_extract+line[22:]
            out.write(line)
    out.close()


actualpdbnames = argv[1:-1]
chain_to_extract = argv[-1]
assert( len( chain_to_extract) == 1)

for actualpdbname in actualpdbnames:
    newpdbfile = actualpdbname.replace('.pdb',chain_to_extract+'.pdb')

    out = open( newpdbfile, 'w' )

    print 'Extracting to ',newpdbfile,'...'

    extractchain(actualpdbname, out, chain_to_extract)

