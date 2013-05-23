#!/usr/bin/python

from os import popen,system
from os.path import exists,dirname,basename,expanduser
import sys
import string
from glob import glob

indir = sys.argv[1]
outdir = sys.argv[2]

PYDIR = expanduser('~rhiju')+'/python/'
assert( exists( PYDIR ) )

inputres = 0
if len(sys.argv)>4:
    startseq = int(sys.argv[3])
    endseq = int(sys.argv[4])
    inputres = 1


newprefix = 'truncate_termini_'
if len(sys.argv)>5:
    newprefix = sys.argv[5]

command = 'mkdir '+outdir
print(command)
system(command)

if not inputres:
    secstructprobfile = glob(indir+'/*.secstructprob')
    outfile = outdir+'/truncate_sequence.txt'
    assert(len(secstructprobfile)>0)
    command = PYDIR+'/decide_termini_truncate.py '+secstructprobfile[0]+ ' ' + outfile
    print(command)
    system(command)

    assert( exists( outfile))
    line = open(outfile).readlines()
    cols = string.split(line[0])
    startseq = int(cols[0])
    endseq   = int(cols[1])

print
print 'Using start and end residues: ',startseq,endseq
print


infile = glob(indir+'/*.pdb')
if(len(infile)>0): # PDB file is optional.
    infile = infile[0]
    outfile = outdir + '/'+newprefix+basename(infile)
    command = PYDIR+'/termini_truncate_pdb.py %s %d %d %s' % \
              (infile,startseq,endseq,outfile)
    print(command)
    system(command)
else:
    print 'COULD NOT FIND PDB FILE BUT THAT IS OK IF YOU ARE DOING CASP.'

infile = glob(indir+'/*.fasta*')
assert(len(infile)>0)
infile = infile[0]
outfile = outdir + '/'+newprefix+basename(infile)
command = PYDIR+'/termini_truncate_fasta.py %s %d %d %s' % \
          (infile,startseq,endseq,outfile)
print(command)
system(command)

infile = glob(indir+'/*.psipred_ss2*')
assert(len(infile)>0)
infile = infile[0]
outfile = outdir + '/'+newprefix+basename(infile)
command = PYDIR+'/termini_truncate_psipred_ss2.py %s %d %d %s' % \
          (infile,startseq,endseq,outfile)
print(command)
system(command)


infiles = glob(indir+'/*v1_3*')
assert(len(infiles)>1)
for infile in infiles:
    outfile = outdir + '/'+newprefix+basename(infile)
    if basename(infile)[:6] == 'boinc_': # A special case.
        outfile = outdir + '/boinc_'+newprefix + basename(infile)[6:]

    command = PYDIR+'/termini_truncate_fragfile.py %s %d %d %s' % \
              (infile,startseq,endseq,outfile)
    print(command)
    system(command)



