#!/usr/bin/python

from sys import argv
from os import system,chdir,getcwd
from os.path import basename,dirname,abspath
import string

filename = argv[ 1 ]
if filename[-5:] == 'fasta':
    tmp_file_name = '/tmp/'+basename( filename )
    system( 'cp '+filename+' '+tmp_file_name)
else:
    tmp_file_name = '/tmp/blah.fasta'
    sequence = string.join( argv[1:] )
    fid = open( tmp_file_name , 'w' )
    fid.write( sequence+'\n')
    fid.close()

CWD = getcwd()
chdir( '/tmp')
command = 'UNAFold.pl  '+basename( tmp_file_name)
system( command )

SIR_GRAPH_EXE = '/Users/rhiju/src/mfold_util-3.3/src/sir_graph/sir_graph'
command = SIR_GRAPH_EXE+'  -f -p -o tmp.ps '+basename(tmp_file_name)+'_1.ct'
system( command )

image_name = basename(tmp_file_name).replace('.fasta','.png')
command = 'convert tmp.ps '+image_name
system( command )

command = 'mv '+image_name+' '+CWD
system( command )

chdir( CWD )
command = 'open -a /Applications/Preview.app '+image_name
system( command )
