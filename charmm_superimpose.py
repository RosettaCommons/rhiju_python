#!/usr/bin/python
# Uses tinker superimpose to get all-atom rmsd
# between two pdbs.


from sys import argv,exit
from os.path import exists,expanduser
from os import system, popen
import string

native_file = argv[1]
pdb_files_in = argv[2:]

HOMEDIR = expanduser('~')

pdb_files_convert = []
for pdb_file in pdb_files_in:
    pdb_file_convert = string.lower(pdb_file).replace( '.pdb', '_RNA.pdb' )

    if not exists( pdb_file_convert ):
        command = 'python "+HOMEDIR+"/python/make_rna_rosetta_ready.py '+pdb_file
        system( command )

    pdb_files_convert.append( pdb_file_convert )


MINI_ROSETTA_EXE = HOMEDIR+'/src/mini/bin/rna_test.macosgccrelease'
if not exists( MINI_ROSETTA_EXE ):
    MINI_ROSETTA_EXE = HOMEDIR+'/src/mini/bin/rna_test.linuxgccrelease'
assert( exists( MINI_ROSETTA_EXE ) )

command = MINI_ROSETTA_EXE+ ' -calc_rmsd -mute all -database ~/minirosetta_database -native '+native_file + ' -s '+string.join( pdb_files_convert )

lines = popen( command ).readlines()
rmsd = {}
for line in lines:
    cols = string.split( line )
    if len( cols ) > 2 and  cols[0] == 'RMSD' :
        print cols[1], cols[2]
        rmsd[ cols[1] ] = cols[2]


for pdb_file in pdb_files_in:
    pdb_file_convert = string.lower(pdb_file).replace( '.pdb', '_RNA.pdb' )
    if pdb_file_convert in rmsd.keys():
        fid = open( pdb_file+'.rms.txt','w')
        fid.write( ' RMSD : %10.6f\n' %  float(rmsd[ pdb_file_convert]) )
        fid.close()
    command = 'rm -rf '+pdb_file_convert
    system( command )



