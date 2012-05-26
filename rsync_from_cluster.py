#!/usr/bin/python

from sys import argv,exit
import string
from os import system
from os.path import basename,dirname,abspath,exists
from cluster_info import *

def Help():
    print
    print argv[0]+' <cluster> <any extra rsync flags>'
    print
    exit()

if len(argv)<2:
    Help()

cluster_in = argv[1]
(cluster,remotedir) = cluster_check( cluster_in )

if cluster == 'unknown':
    Help()

extra_args = argv[2:]

dir = '.'
clusterdir = abspath(dir).replace('/Users/rhiju/','')
clusterdir = abspath(dir).replace('/scratch/users/rhiju/','')
clusterdir = clusterdir.replace('/work/rhiju/','')

clusterdir = remotedir+clusterdir

#if cluster[:3]=='syd':
#    n = cluster[3]
#    cluster = 'syd'
#    clusterdir = 'work'+n+'/'+clusterdir

#command = 'ssh ' + cluster + ' mkdir -p '+clusterdir
#print(command)
#system(command)


cluster_prefix = cluster+':'
if len(cluster) == 0: cluster_prefix = ''

command = 'rsync -avzL '+cluster_prefix+clusterdir+'/'+string.join(extra_args)+' '+dir+' --exclude="condor*log"'
print(command)
system(command)

