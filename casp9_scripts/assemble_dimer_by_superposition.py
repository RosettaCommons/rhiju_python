#!/usr/bin/python

from sys import argv
from os.path import basename,exists,expanduser
from os import system
from get_sequence import get_sequence

monomer_pdb = argv[1]
dimer_pdb = argv[2]
dimer_chain1 = argv[3]
dimer_chain2 = argv[4]

assert( exists( monomer_pdb ) )
assert( exists( dimer_pdb ) )
assert( monomer_pdb.count( '.pdb' ) == 1 )
assert( dimer_pdb.count( '.pdb' ) == 1 )
assert( len( dimer_chain1) == 1 )
assert( len( dimer_chain2) == 1 )
dimer_chains = [ dimer_chain1, dimer_chain2 ]

HOMEDIR = expanduser( '~' )
PYDIR = HOMEDIR + '/python/'
assert( exists( PYDIR ) )

count = 0
chain_extract_pdbs = []
for dimer_chain in dimer_chains:
    command = '%s/extract_chain.py  %s %s' % ( PYDIR, dimer_pdb, dimer_chain )
    print( command )
    system( command )

    chain_extract_pdb = basename( dimer_pdb ).replace('.pdb',dimer_chain+'.pdb')
    chain_extract_pdbs.append( chain_extract_pdb )
    assert( chain_extract_pdb )
    assert( len(get_sequence( chain_extract_pdb)) > 0 )

    count += 1

sup_pdbs = []
for count in range(2):

    sup_pdb = "tmp_align%s.pdb" % dimer_chains[count]

    command = "%s/superimpose.py %s %s -R 2.0 > %s" % (PYDIR,chain_extract_pdbs[count],monomer_pdb,sup_pdb)
    print( command )
    system( command )

    assert( exists( sup_pdb ) )
    command = "%s/parse_NMR_models.py %s" % (PYDIR, sup_pdb )
    print( command )
    system( command )

    sup_pdbs.append( sup_pdb )

monomer_dimerized = basename(monomer_pdb).replace('.pdb','_DIMER.pdb' )
command = "cat %s %s > %s" % \
    ( sup_pdbs[0].replace('.pdb','_002.pdb'),\
          sup_pdbs[1].replace('.pdb','_002.pdb'),\
          monomer_dimerized )
print command
system( command )

command = "renumber_pdb_in_place.py "+monomer_dimerized
print command
system( command )

for file in chain_extract_pdbs: system( 'rm -rf '+file )
for file in sup_pdbs:
    system( 'rm -rf '+file )
    system( 'rm -rf '+file.replace('.pdb','_001.pdb' ) )
    system( 'rm -rf '+file.replace('.pdb','_002.pdb' ) )

print
print "Created: ", monomer_dimerized
print


