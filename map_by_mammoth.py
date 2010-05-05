#!/usr/bin/python

from sys import argv,stderr
from os.path import exists,expanduser
from os import system
import string

main_file = argv[1]
other_files = argv[2:]

HOMEDIR = expanduser( '~rhiju' )
MAMMOTH_EXE = HOMEDIR + '/mammoth2/mammoth_rna' # or whatever
assert( exists( MAMMOTH_EXE ) )

MAMMOTH2SEQUENCE = HOMEDIR + '/python/mammoth2sequence.py'
assert( exists( MAMMOTH2SEQUENCE ) )

for file in other_files:

    tmp_file = 'tmp.mammoth'
    command = '%s -p %s -e %s > %s' % (MAMMOTH_EXE , main_file, file, tmp_file )
    system( command )

    if file.find( '.pdb' ) > 0:
        seqfile = file.replace('.pdb','.mapping')
    else:
        seqfile = file + '.mapping'

    command = '%s %s > %s' % (MAMMOTH2SEQUENCE, tmp_file, seqfile )
    system( command )
    stderr.write( 'Generated ... %s\n' % seqfile )
