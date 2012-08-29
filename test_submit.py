#!/usr/bin/python

from sys import argv
from os import system,chdir,getcwd
from os.path import exists
import string

file = argv[1]

do_all = 0
if argv.count('-all'):
    do_all = 1

fast = 0
if argv.count('-fast'):
    fast = 1

data = open(file,'r')

command = 'rm -rf test; mkdir test'
print command
system(command)

chdir('test')

line = data.readline()

keeptesting = 1
count = 1
while line and keeptesting:

    # Input files first
    while line and line[:10] != 'inputfiles': line = data.readline()
    if not line: break
    inputfiles_string = string.split(line,'=')[-1][1:-1]
    inputfiles = string.split(inputfiles_string,';')

    for inputfile in inputfiles:
        command = 'cp '+inputfile+' .'
        print(command)
        print inputfile
        assert( exists(inputfile) )
        system(command)


    line = data.readline()
    while line and line[:9] != 'arguments': line = data.readline()
    if not line: break

    args = string.split(line,'=')[-1][:-1]

    if do_all: args=args.replace(' xx',' %02d' % count)
    if fast:
        args=args.replace('-abrelax','')
        args=args.replace('-farlx','')

    args = args.replace('-output_silent_gz','')
    args = args.replace('-close_chainbreaks','')

    cols = string.split(args)
    pos = cols.index('-increase_cycles')
    cols[pos+1] = '0.1'
    args = string.join(cols,' ')

    args = args.replace('nstruct 30','nstruct 1')
    args += ' -farlx_cycle_ratio 0.01'
#    args += ' -benchmark'
#    args += ' -use_fold_constraints_no_minimize'
    args += ' -paths /users/rhiju/paths.txt'


    print 'CURRENT DIRECTORY ',getcwd()

    command = '/users/rhiju/rosetta++/rosetta.gcc '+args
#    command = '/users/bqian/rosetta.gcc '+args
    print(command)
    system(command)

    if not do_all:
        keeptesting = 0
    count=count+1

