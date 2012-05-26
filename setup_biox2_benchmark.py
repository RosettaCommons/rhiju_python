#!/usr/bin/python

from sys import argv,exit
from os.path import exists
from os import system
from parse_options import parse_options

torsions_file = parse_options( argv, 'vall_torsions', '' )

infiles = argv[1:]

fid = open( 'bsubMINI', 'w' )

EXE = '/home/rhiju/src/rosetta_TRUNK/rosetta_source/bin/rna_denovo.linuxgccrelease'
DB = '/home/rhiju/src/rosetta_TRUNK/rosetta_database/'

if ( not infiles[-4:] == '.pdb' ):
    lines = open( infiles[0] ).readlines()
    infiles = []
    for line in lines:
        infiles.append( line[:-1] )

#fid.write( 'DELAY =$$([ $(Process)*5] )\n\n')

OUTFILE_DIR = 'OUTFILES'
QUEUE_NUM = 10;
NSTRUCT = 200;

process_num  = -1
for file in infiles:
    fasta_file = file

    assert( file[-6:] == '.fasta' )
    fourlettercode = fasta_file[ -11:-7 ]
    fullcode = fasta_file[:-7]
    native_file = fullcode+'_RNA.pdb'

    IN_PATH = '../bench_final/'


    print fasta_file
    assert( exists( IN_PATH + fasta_file ) )
    params_file = fullcode + '_.prm'

    if 1:

        for queue in range( QUEUE_NUM):
            process_num += 1
            system( 'mkdir -p %s/%d' % (OUTFILE_DIR, process_num ) )

            decoys_silent_file = fullcode + '_5kcycles_decoys_nonativefrags.out'
            command = 'bsub -W 24:0 -o /dev/null -e /dev/null %s  -database %s -native %s -fasta %s -params_file %s -nstruct %d -out::file::silent %s/%s/%s -minimize_rna -cycles 5000 -mute all -filter_lores_base_pairs -vary_geometry  -output_lores_silent_file -in:path %s  ' % \
                      (  EXE, DB, native_file, fasta_file, params_file, NSTRUCT, OUTFILE_DIR, process_num, decoys_silent_file , IN_PATH )
            if len(torsions_file) > 0: command += ' -vall_torsions %s' % torsions_file
            fid.write( command +'\n' )

    fid.write( '\n')

fid.close()
