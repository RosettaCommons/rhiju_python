#!/usr/bin/python


from sys import argv
from os import chdir,getcwd,popen,system
from os.path import exists
from parse_options import parse_options

############################
CUTOFF = parse_options( argv, 'coord_cutoff', -1.0 )
segment_residues = parse_options( argv, 'segments', [ -1 ] )
start_file  = parse_options( argv, 's', '' )
native_file = parse_options( argv, 'native', '' )
loop = parse_options( argv, 'loop', [ -1 ] )
NJOBS =  parse_options( argv, 'j', 400 )

assert( len( start_file ) > 0 )
assert( exists( start_file ) )
assert( len( native_file ) > 0 )
assert( exists( native_file ) )
assert( len ( loop ) == 2 )
loop_start = loop[ 0 ]
loop_end = loop[ 1 ]

subset_residues = []
if len( segment_residues ) > 0:
    for i in range( len(segment_residues)/2):
        for j in range( segment_residues[2*i], segment_residues[2*i+1]+1 ):
            subset_residues.append( j )
    print subset_residues

lines = popen( 'pdb2fasta.py '+start_file ).readlines()
numres = int( lines[0].split()[-1] )
pdb_sequence = lines[1][:-1]


lines = popen( 'pdb2fasta.py '+native_file ).readlines()
native_sequence = lines[1][:-1]
# Important for now!
assert( pdb_sequence == native_sequence )

if len( subset_residues ) == 0: subset_residues = range( 1, numres+1 )

assert( loop_start in subset_residues )
assert( loop_end   in subset_residues )

# Make directory
if ( native_file == start_file ):
    dir_name = 'loop%d-%d/from_native' % ( loop_start, loop_end )
else:
    dir_name = 'loop%d-%d/from_template' % ( loop_start, loop_end )

system( 'mkdir -p '+dir_name )
CWD = getcwd()

# copy in PDB, go into directory.
for pdb_file in [start_file, native_file]:
    system( 'rsync -avz %s %s ' % (pdb_file, dir_name ) )
chdir( dir_name )

# Reduce, and strip out sidechains.
for pdb_file in [start_file, native_file]:

    system( 'renumber_pdb_in_place.py '+pdb_file )

    if pdb_file.count('_H.pdb') > 0:
        pdb_H = pdb_file
    else:
        pdb_H = pdb_file.replace('.pdb','_H.pdb')
        system( 'reduce %s > %s ' % (pdb_file, pdb_H) )

    # Slice out pdb file.
    command = 'pdbslice.py %s -subset ' % pdb_H
    for k in subset_residues: command += ' %d' % k
    command += ' mini_ '
    print( command )
    system( command )

    if ( pdb_file == start_file ):  mini_start_file = 'mini_'+pdb_H
    if ( pdb_file == native_file ): mini_native_file = 'mini_'+pdb_H


# In the future, could input fasta instead of start pdb.
fasta_file = mini_start_file.replace('.pdb','.fasta')
system( 'pdb2fasta.py %s > %s ' % ( mini_start_file, fasta_file ) )


pdb_stripsidechain = mini_start_file.replace('.pdb','_stripsidechain.pdb')
system( '~/python/strip_sidechain.py ' + mini_start_file )

# Slice out "noloop" file.
command = 'pdbslice.py %s -excise ' % pdb_stripsidechain
loop_res = range( loop_start, loop_end+1)
for k in loop_res: command += ' %d' % k
command += ' noloop_ '
print( command )
system( command )
noloop_file = 'noloop_'+pdb_stripsidechain

# Prepare "README_SETUP"
input_res = []
working_loop_res = []
for k in range( len(subset_residues) ):
    if ( subset_residues[k] not in loop_res ):
        input_res.append( k+1 )
    else:
        working_loop_res.append( k+1 )

input_res_tag = ''
for k in input_res: input_res_tag += ' %d' % k
working_loop_res_tag = ''
for k in working_loop_res: working_loop_res_tag += ' %d' % k

fid = open( 'README_SETUP', 'w' )
fid.write( 'rm -rf CONDOR* SLAVE*  core*  REGION*  # region*\n' )

cst_file_tag = ''
if CUTOFF > 0:
    cst_file =  mini_start_file.replace( '.pdb','_coordinate%3.1f.cst' % CUTOFF )
    system( 'generate_coordinate_constraints.py %s -fade -stdev %3.1f -cst_res %s -anchor_res 1  > %s\n' % ( mini_start_file, CUTOFF, working_loop_res_tag, cst_file ) )
    cst_file_tag =  '-cst_file ' + cst_file

fid.write( 'protein_loop_build_dagman.py -s %s  -input_res %s  -native %s -fasta %s  -score_diff_cut 10 -nstruct 200  -weights score12_no_hb_env_dep.wts -pack_weights pack_no_hb_env_dep.wts  -one_loop  -cluster_radius_sample 0.1 -cluster_radius 0.25  %s  -align_pdb %s \n' % ( noloop_file, input_res_tag, mini_native_file, fasta_file, cst_file_tag, mini_start_file ) )
fid.close()

# Prepare "README_SUB"
fid = open( 'README_SUB', 'w' )
fid.write(' rm blah.* \n' )
fid.write( 'bsub -W 24:0 -o blah.out -e blah.err SWA_pseudo_dagman_continuous.py  -j %d rna_build.dag  \n' % NJOBS )
fid.close()

chdir( CWD )

print
print "Setup directory: ",dir_name

