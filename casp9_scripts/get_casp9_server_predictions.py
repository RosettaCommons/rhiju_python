#!/usr/bin/python

from sys import argv
from os.path import exists,expanduser
from os import system,chdir

target_num = int( argv[1] )
HOMEDIR = expanduser('~')

if ( target_num <= 386 ):   casp_num = 7
elif ( target_num <= 514 ): casp_num = 8
else:                       casp_num = 9

prefix = 'http://predictioncenter.org/download_area/CASP%d/server_predictions/' % casp_num
savedir = HOMEDIR +'/projects/casp9/t%04d/' % (target_num)
if not exists( savedir ):
    print 'Making directory: ' + savedir
    system( 'mkdir -p %s' % savedir )

tarfile = 'T%04d.3D.srv.tar.gz' % target_num
command  =  'curl %s/%s > %s/%s' % (prefix,tarfile,savedir,tarfile)
print command
system( command )

chdir( savedir)
command = 'tar xvfz %s ' % tarfile
print command
system( command )

command = 'rm  %s ' % tarfile
print command
system( command )

command = 'mv %s %s' % (tarfile.replace('.3D.srv.tar.gz',''), 'server_predictions' )
print command
system( command )

print
print 'Server predictions for T%04d are in %s/server_predictions/' % (target_num, savedir )


