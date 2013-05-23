#!/usr/bin/python

from os import system,getcwd,chdir
from glob import glob
from os.path import basename,exists
import string
from sys import argv
from parse_options import parse_options
from make_tag import make_tag
from get_sequence import get_sequence

loop_res = parse_options( argv, "loop_res", [-1] )
skip_res = parse_options( argv, "skip_res", [-1] )
input_res = parse_options( argv, "input_res", [-1] )
insert_res = parse_options( argv, "insert_res", [-1] )
fasta = parse_options( argv, "fasta", "" )
tag = parse_options( argv, "tag", "domain_insertion" )
input_pdb = parse_options( argv, "s", "" )
insert_pdb = parse_options( argv, "insert", "")
frag_files = parse_options( argv, "frag_files",[""] )
centroid = parse_options( argv, "centroid",0 )
endpoints = parse_options( argv, "endpoints", [-1] )

#print argv
assert( len (argv ) == 1 )

if input_res[0] > 1:
    print "Hey need to fix up this script a bit to allow for input_res that does not have Nterminus! E.g., frag files need to be trimmed."
    exit( 0 )

if not exists( tag ): system( "mkdir -p "+tag )

#########################################################
if len ( frag_files ) > 0:
    for file in frag_files:
        if not exists( tag+'/'+basename( file ) ):
            command = 'rsync -avzL '+file+' '+tag
            print command
            system( command )
    frag_files = map( lambda x:basename(x), frag_files )

for file in [fasta, input_pdb, insert_pdb]:
    assert( exists( file ) )
    if not exists( tag + '/' + basename( file ) ):
        system( 'rsync -avzL %s %s ' % ( file, tag ) )
chdir( tag )



input_pdb = basename( input_pdb )
insert_pdb = basename( insert_pdb )
fasta = basename( fasta )

########################################################
# Do sequence checks
########################################################
sequence_lines = open( fasta  ).readlines()[1:]
sequence = string.join(  map( lambda x : x[:-1], sequence_lines) ,  '' )
NRES = len( sequence )

sequence_base = get_sequence( input_pdb )
sequence_insert = get_sequence( insert_pdb )

sequence_sub = ""
for m in input_res: sequence_sub += sequence[ m-1 ]
assert( sequence_base == sequence_sub )

sequence_sub = ""
for m in insert_res: sequence_sub += sequence[ m-1 ]
assert( sequence_insert == sequence_sub )


#########################################################
if len( skip_res ) == 0:
    insert_res.sort()
    for i in range( 1, len( insert_res )-1 ): skip_res.append( insert_res[i] )
    print "Will use skip_res of ", make_tag( skip_res )

#########################################################
all_model_res = []
for m in loop_res:    all_model_res.append( m )
for m in input_res:   all_model_res.append( m )
all_model_res.sort()

model_res = [] # All residues to be modeled ... without any repeats
excise_res = [] # Any redundancy? Should cut out of starting conformation
for m in all_model_res:
    if m in model_res:
        if m not in excise_res: excise_res.append( m )
    else:
        model_res.append( m )

print "EXCISE_RES:", excise_res

for m in model_res: assert( model_res.count( m ) == 1 )
for i in range( len( model_res) - 1 ): assert( model_res[i]+1 == model_res[i+1] )

#########################################################
noloop_input_pdb = "noloop_" + input_pdb
system( "renumber_pdb_in_place.py "+input_pdb )
working_excise_res = []
input_res_sliced = []
for i in range( len( input_res ) ):
    if input_res[i] in excise_res:
        working_excise_res.append( i+1 )
    else:
        input_res_sliced.append( input_res[i] )

if len( working_excise_res ) == 0:
    system( 'cp %s %s' % (input_pdb, noloop_input_pdb) )
else:
    command =  'pdbslice.py %s -excise %s noloop_' % ( input_pdb, make_tag( working_excise_res )  )
    print command
    system( command )

working_input_sequence_sliced = get_sequence( noloop_input_pdb )
sequence_sub = ""
for m in input_res_sliced: sequence_sub += sequence[ m-1 ]
assert( working_input_sequence_sliced == sequence_sub )



#########################################################
if ( not len( model_res ) == NRES):
    assert len( model_res ) < NRES
    full_fasta = fasta
    fasta = "subsequence_" + fasta
    fid = open( fasta, 'w' )
    fid.write( '>Domain insertion run for '+basename( input_pdb)+' and '+basename( insert_pdb)+'\n' )
    for m in model_res:   fid.write( '%s' % sequence[m-1] )
    fid.write( '\n' )
    fid.close()


#########################################################
# Constraint file setup
#########################################################
for i in range( len( input_res_sliced ) ):
    if ( input_res_sliced[i] > loop_res[0] ):
        input_res_Nterm_to_loop = input_res_sliced[ i - 1]
        input_res_Cterm_to_loop = input_res_sliced[ i ]
        break

fixed_insert_res = []
for m in insert_res:
    if ( m not in skip_res ): fixed_insert_res.append( m )
fixed_insert_res.sort()
insert_res_near_Nterm_of_loop = fixed_insert_res[ 0 ]
insert_res_near_Cterm_of_loop = fixed_insert_res[ 1 ]

cst_file = 'close_chain_harmonic.cst'
fid = open( cst_file, 'w' )
fid.write( '[ atompairs ]\n' )
fid.write( ' CA %3d CA %3d  HARMONIC 0.0 4.0\n' % ( input_res_Nterm_to_loop, insert_res_near_Nterm_of_loop) )
fid.write( ' CA %3d CA %3d  HARMONIC 0.0 4.0\n' % ( input_res_Cterm_to_loop, insert_res_near_Cterm_of_loop) )
fid.close()
#########################################################

fid = open( 'README_SETUP','w')
fid.write('rm -rf CONDOR* SLAVE*  core*  REGION*  # region*\n')

run_type = " -denovo 1 "
if len( frag_files ) > 0:
    run_type = " -frag_files "+string.join( frag_files )
if centroid: run_type += " -centroid -score_diff_cut 200 -cluster_radius 1.0"

command = 'grinder_dagman.py  -loop_start_pdb %s  -loop_res %s -fasta %s  -nstruct 50  -skip_res %s   -s %s  -input_res %s  %s '  % \
                  ( noloop_input_pdb,
                    make_tag( loop_res ),
                    fasta,
                    make_tag( skip_res ),
                    insert_pdb,
                    make_tag( insert_res ),
                    run_type )
if len( endpoints ) > 0 : command += ' -endpoints '+make_tag( endpoints )
command += ' -cst_file ' + cst_file

fid.write( command + '\n' )


fid.close()

fid = open( 'README_SUB','w')
fid.write( 'rm blah.* \n')
fid.write('bsub -W 96:0 -o blah.out -e blah.err SWA_pseudo_dagman_continuous.py  -j 100 protein_build.dag\n')
fid.close()

#chdir( cwd )


