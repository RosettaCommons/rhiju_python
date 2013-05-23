#!/usr/bin/python

from get_sequence import get_sequence
from make_tag import make_tag
from os import system
from sys import argv

pdbs = argv[1:]

def new_tag( file, change_pdb_to_this_tag ):
    if ( file.count( '.pdb' ) ):
        return file.replace( '.pdb', change_pdb_to_this_tag )
    else:
        return file + change_pdb_to_this_tag

for pdb in pdbs:
    cst_file = new_tag( pdb ,'.cst' )
    command = 'generate_CA_constraints.py  %s > %s ' % (pdb, cst_file )
    #command = 'generate_coordinate_CA_constraints.py  %s -anchor_res 1 > %s ' % (pdb, cst_file )
    print command
    system( command )

fid = open( 'README_MINIMIZE', 'w' )
for pdb in pdbs:

    fasta_file = new_tag( pdb, '.fasta')
    outfile = new_tag( pdb, '_minimize.out')
    cst_file = new_tag( pdb, '.cst' )

    sequence = get_sequence( pdb )
    nres = len( sequence )
    system( 'pdb2fasta.py %s > %s' % (pdb, fasta_file) )

    command = 'stepwise_protein_test.macosgccrelease  -database ~/minirosetta_database/ -s1 %s  -input_res1 %s  -use_packer_instead_of_rotamer_trials -global_optimize -cst_file  %s  -out:file:silent %s -fasta %s -score:weights score12_no_hb_env_dep.wts -pack_weights pack_no_hb_env_dep.wts  -align_pdb %s ' % ( pdb, make_tag( range(1,nres+1) ), cst_file, outfile, fasta_file, pdb )

    fid.write( command + '\n\n' )

fid.close()

###############################
# Early exit -- no relax
###############################
exit( 0 )


fid = open( 'README_RELAX', 'w' )
for pdb in pdbs:

    outfile = new_tag( pdb, '_relax.out')
    cst_file = new_tag( pdb, '.cst' )

    command = 'relax.macosgccrelease  -database ~/minirosetta_database/ -s %s  -constraints:cst_fa_file  %s  -out:file:silent %s -score:weights score12_no_hb_env_dep.wts -relax:fast -fastrelax_repeats 32 -nstruct 1' % ( pdb, cst_file, outfile )

    fid.write( command + '\n\n' )


fid.close()
