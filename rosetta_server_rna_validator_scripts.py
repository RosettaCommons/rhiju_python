#!/usr/bin/python

from sys import stderr, stdout
from os.path import exists
import string
from rna_conversion import make_rna_rosetta_ready
from get_sequence import get_sequence

def convert_fasta_to_rosetta_format( fasta_file_name ):
    # read lines from fasta file
    if not exists( fasta_file_name ): return None
    lines = open( fasta_file_name ).readlines()


    line = lines[0]
    if not line[0] == '>':
        stderr.write( 'First line of fasta must start with \'>\'\n' )
        return None

    sequence = string.join( lines[1:] ).replace( ' ','').replace('\n','')
    sequence = sequence.lower()

    # quick validation
    goodchars = ['a','c','g','u']
    for char in sequence:
        if char not in goodchars:
            stderr.write( 'Character %s is not a,c,g, or u\n'% char )
            return None

    # here's the output
    outstring = ''
    outstring += line
    outstring += sequence
    outstring += '\n'

    return outstring


def convert_pdb_to_rosetta_format( pdb_file_name ):

    outstring = make_rna_rosetta_ready( pdb_file_name )

    return outstring


def does_PDB_match_fasta( PDB_file_name, fasta_file_name ):
    seq_PDB = get_sequence( PDB_file_name )
    seq_fasta = string.join( open( fasta_file_name ).readlines()[1:] ).replace( '\n', '' ).replace(' ','' )
    return ( seq_PDB == seq_fasta )
