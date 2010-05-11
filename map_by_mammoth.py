#!/usr/bin/python

from sys import argv,stderr
from os.path import exists,expanduser
from os import system,popen
import string
from get_sequence import get_sequence

main_file = argv[1]
other_files = argv[2:]

HOMEDIR = expanduser( '~rhiju' )
MAMMOTH_EXE = HOMEDIR + '/mammoth2/mammoth_rna' # or whatever
assert( exists( MAMMOTH_EXE ) )

MAMMOTH2SEQUENCE = HOMEDIR + '/python/mammoth2sequence.py'
assert( exists( MAMMOTH2SEQUENCE ) )

def add_last_residue( line, last_residue ):
    for i in range( len( line ) - 2 ):
        pos = len(line) - 2 - i
        if line[ pos ] != '-':
            line = line[:(pos+1)] + last_residue + line[(pos+1):]
            break
    return line

def silly_fix_for_last_residue( lines, main_sequence, sequence ):
    lines[0] = add_last_residue( lines[0], main_sequence[-1] )
    lines[1] = add_last_residue( lines[1], sequence[-1] )


main_sequence = get_sequence( main_file )
for file in other_files:

    tmp_file = 'tmp.mammoth'
    command = '%s -p %s -e %s > %s' % (MAMMOTH_EXE , main_file, file, tmp_file )
    system( command )

    if file.find( '.pdb' ) > 0:
        seqfile = file.replace('.pdb','.mapping')
    else:
        seqfile = file + '.mapping'

    command = '%s %s ' % (MAMMOTH2SEQUENCE, tmp_file )
    lines = popen( command ).readlines()

    sequence = get_sequence( file )
    silly_fix_for_last_residue( lines, main_sequence, sequence )

    fid = open( seqfile, 'w')
    fid.write( lines[0] )
    fid.write( lines[1] )
    fid.close()

    stderr.write( 'Generated ... %s\n' % seqfile )
