#!/usr/bin/python

from os import system,getcwd,chdir
from glob import glob
from os.path import basename,exists,expanduser
import string
from sys import argv
from parse_options import parse_options
from make_tag import make_tag
from get_sequence import get_sequence

EXE = expanduser('~rhiju') + '/src/mini/bin/stepwise_protein_test.macosgccrelease'
DB = expanduser('~rhiju') + '/minirosetta_database'
EXTRACT = expanduser('~rhiju')+'/python/extract_lowscore_decoys.py'

pack_pdbs = argv[1:]

for pdb in pack_pdbs:

    sequence = get_sequence( pdb )
    NRES = len( sequence )
    res_range = range( 1, NRES+1 )

    fasta = pdb+'.fasta'
    fid = open( fasta, 'w' )
    fid.write( '> '+pdb+'\n')
    fid.write( sequence + '\n' )
    fid.close()

    outfile = pdb+'.PACK.out'

    command = '%s -database %s -s1 %s -fasta %s  -global_optimize -fixed_res %s  -input_res1 %s -pack_weights pack_no_hb_env_dep.wts -score:weights score12_no_hb_env_dep.wts -align_pdb %s -out:file:silent %s -use_packer_instead_of_rotamer_trials ' % (EXE, DB,pdb,fasta,make_tag(res_range), make_tag(res_range),pdb,outfile)

    print command
    system( command )

    assert( exists( outfile ) )

    command = '%s %s 1' % (EXTRACT,outfile)

    print command
    system( command )

    system( 'rm -rf '+fasta )
    system( 'rm -rf '+outfile )
