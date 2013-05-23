#!/usr/bin/python

from sys import argv,exit
from os import system,chdir
from os.path import exists
from glob import glob

target_num = int( argv[1] )

workdir = "/home/rhiju/projects/casp9/t%04d/rosetta_frags/" % target_num

if exists( workdir ):
    print 'You already have something in: ', workdir
    chdir( workdir )
else:
    system( "mkdir -p "+workdir )

    chdir( workdir )

    kyledir = '/home/kyleb/casp9/frags/T%04d/seq/' %  target_num
    if exists( kyledir ):
        print 'Found fragment directory in kyle''s stuff!'
        system( 'rsync -avzL %s/ %s' % (kyledir, '.' ) )

fasta_file = "t%3d_.fasta" % target_num
if not exists( fasta_file ):
    fasta_file = "T%04d.fasta" % target_num

if not exists( fasta_file ):
    command = "wget 'http://predictioncenter.gc.ucdavis.edu/casp9/target.cgi?target=T%04d&view=sequence'" % target_num
    system( command )

    command = "mv 'target.cgi?target=T%04d&view=sequence' %s" % (target_num, fasta_file )
    system( command )
assert( exists( fasta_file ) )


frag_file3 = 'aat%3d_03_05.200_v1_3'  % target_num
if not exists( frag_file3)  and not exists( frag_file3 + '.gz' ):
    frag_file3 = 'aaT%04d03_05.200_v1_3'  % target_num

if not exists( frag_file3)  and not exists( frag_file3 + '.gz' ):
    FRAG_EXE = '/home/kyleb/fragments/nnmake_new/make_fragments.pl'
    assert( exists( FRAG_EXE ) )
    print "Launching frag picker..."
    command  = "nohup %s  -verbose %s > make_frags.log 2> make_frags.err &" % (FRAG_EXE, fasta_file )
    print command
    system( command )

globfiles = glob(  'aa*200_v1_3' )
globfiles.sort()
for globfile in globfiles:
    command = 'gzip '+globfile
    print command
    system( command )


