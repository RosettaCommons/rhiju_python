#!/usr/bin/python

from sys import argv,exit
import string
from os import system
from os.path import basename,dirname,abspath,exists
from cluster_info import *

def Help():
    print
    print argv[0]+' <cluster> <command you want to run>'
    print
    exit()


cluster_in = argv[1]
(cluster,remotedir) = cluster_check( cluster_in )
if len( cluster ) == 0:
    Help()

extra_args = argv[2:]

dir = '.'
clusterdir = abspath(dir).replace('/Users/rhiju/','')
clusterdir = clusterdir.replace('/work/rhiju/','')

clusterdir = remotedir+clusterdir


command_to_run =  string.join(extra_args)

# This is really dumb. I cannot seem to get ssh to apply aliases from the .bashrc. Weird.
if command_to_run == 'cleanup':
    command_to_run = 'find ./ -name ''*REGI*'' | xargs rm -rf; find ./ -name ''*SLAVE*'' | xargs rm -rf;  find ./ -name ''*CONDOR*'' | xargs rm -rf'



command = 'ssh '+cluster+' "cd '+clusterdir+'; '+command_to_run+'"'
print(command)
system(command)
