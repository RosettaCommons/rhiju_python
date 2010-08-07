#!/usr/bin/python


import string
from sys import argv
from os import chdir,getcwd,popen,system
from os.path import exists, basename
from parse_options import parse_options
from get_sequence import get_sequence
from make_tag import make_tag
############################

start_file = parse_options( argv, 's', '' )
segment_residues = parse_options( argv, 'segments', [ -1 ] )
template_files = parse_options( argv, 'templates', [''] )
loop = parse_options( argv, 'loop', [ -1 ] )
CUTOFF = parse_options( argv, 'coord_cutoff', -1.0 )
tag = parse_options( argv, 'tag', '' )
input_cst_file = parse_options( argv, 'cst_file','' )
NJOBS =  parse_options( argv, 'j', 0 )
final_number =  parse_options( argv, 'final_number', 50 )
loop_force_Nsquared =  parse_options( argv, 'loop_force_Nsquared', 0 )
native_file = parse_options( argv, 'native', '' )
virtual_res = parse_options( argv, 'virtual_res', [-1] )
superimpose_res = parse_options( argv, 'superimpose_res', [-1] )
endpoints = parse_options( argv, 'endpoints', [-1] )
tight_fade = parse_options( argv, 'tight_fade', 0 )
frag_files = parse_options( argv, 'frag_files', [""] )
no_fixed_res = parse_options( argv, 'no_fixed_res', 0 )
filter_native_big_bins = parse_options( argv, 'filter_native_big_bins', 0 )

###############################
# Basic setup and checks.
###############################
assert( len( start_file ) > 0 )
assert( exists( start_file ) )
assert( len( template_files ) > 0 )
for template_file in template_files:
    print template_file
    assert( exists( template_file ) )
assert( len ( loop ) == 2 )
loop_start = loop[ 0 ]
loop_end = loop[ 1 ]

###############################
subset_residues = []
working_cutpoints_open = []
if len( segment_residues ) > 0:
    num_segments = len(segment_residues)/2
    for i in range( num_segments ):
        for j in range( segment_residues[2*i], segment_residues[2*i+1]+1 ):
            subset_residues.append( j )
        if ( i < num_segments-1): working_cutpoints_open.append( len( subset_residues ) )
    #print subset_residues


################################################################
# FragFile -- need to parse out segments?
################################################################
if len( segment_residues ) > 0 and len( frag_files ) > 0:
    frag_files_input = []
    for file in frag_files: frag_files_input.append( file )

    prefix = 'segment'
    for i in segment_residues: prefix += '_%d' % i
    new_frag_dir = prefix + '_frags'

    if not exists( new_frag_dir ): system( 'mkdir '+new_frag_dir )

    frag_files = []
    for file in frag_files_input:
        new_frag_file = '%s/%s_%s' % (new_frag_dir, prefix, basename(file) )
        frag_files.append( new_frag_file )
        if exists( new_frag_file ): continue
        command = 'fragfile_subset.py %s %s %s' % ( file,  new_frag_file, make_tag( subset_residues ) )
        print command
        system( command )


###############################

pdb_sequence = get_sequence( start_file )

# Important for now!
for template_file in template_files:
    template_sequence = get_sequence( template_file )
    assert( pdb_sequence == template_sequence )

numres = len( pdb_sequence )
if len( subset_residues ) == 0: subset_residues = range( 1, numres+1 )

assert( loop_start in subset_residues )
assert( loop_end   in subset_residues )

#######################
# Make directory
#######################
if len( tag ) == 0:
    if ( template_file == start_file ):
        tag = 'from_native'
    else:
        tag = 'from_template'

dir_name = 'loop%d-%d/%s' % ( loop_start, loop_end, tag )
if len( frag_files ) > 0:
    dir_name = 'loop%d-%d_frags/%s' % ( loop_start, loop_end, tag )


system( 'mkdir -p '+dir_name )
CWD = getcwd()

#######################
# PDB setup
#######################
# copy in PDB, go into directory.
all_files = [ start_file ]
for pdb_file in template_files: all_files.append( pdb_file )
if len( native_file ) > 0: all_files.append( native_file )

for pdb_file in all_files:    system( 'rsync -avz %s %s ' % (pdb_file, dir_name ) )
chdir( dir_name )

# Reduce, and strip out sidechains.
mini_files = []
rm_files = []
for pdb_file in all_files:

    system( 'renumber_pdb_in_place.py '+pdb_file )

    if ( pdb_file.count('.pdb') ):
        pdb_H = pdb_file.replace( '.pdb','_H.pdb' )
    else:
        pdb_H = pdb_file + '_H'

    system( 'reduce %s > %s 2> /dev/null ' % (pdb_file, pdb_H) )

    # Slice out pdb file.
    command = 'pdbslice.py %s -subset ' % pdb_H
    for k in subset_residues: command += ' %d' % k
    command += ' mini_ '
    print( command )
    system( command )

    mini_files.append( 'mini_' + pdb_H )

    # Remove all these extra files?
    rm_files.append( pdb_file )
    rm_files.append( pdb_H )


for file in rm_files: system( 'rm -f '+file )


mini_start_file = mini_files[0]
mini_template_files = mini_files[1:]

if len( native_file ) > 0:
    mini_native_file = mini_template_files[-1]
    del( mini_template_files[-1] )

print mini_template_files


# In the future, could input fasta instead of start pdb.
if ( mini_start_file.count('.pdb' ) ):
    fasta_file = mini_start_file.replace('.pdb','.fasta')
else:
    fasta_file = mini_start_file + '.fasta'

system( 'pdb2fasta.py %s > %s ' % ( mini_start_file, fasta_file ) )


######################################################################
# Good hygiene?  Perhaps not necessary
#pdb_stripsidechain = mini_start_file.replace('.pdb','_stripsidechain.pdb')
#system( '~/python/strip_sidechain.py ' + mini_start_file )
# actual_start_pdb = pdb_stripsidechain

actual_start_pdb = mini_start_file

# Slice out "noloop" file.
command = 'pdbslice.py %s -excise ' % actual_start_pdb
loop_res = range( loop_start, loop_end+1)
for k in loop_res: command += ' %d' % k
command += ' noloop_ '
print( command )
system( command )

noloop_file = 'noloop_' + actual_start_pdb


##############################################################
# Residue renumbering
##############################################################
input_res = []
working_loop_res = []
working_superimpose_res = []
working_endpoints = []
for k in range( len(subset_residues) ):
    if ( subset_residues[k] not in loop_res ):
        input_res.append( k+1 )
    else:
        working_loop_res.append( k+1 )

    if subset_residues[k] in superimpose_res:
        working_superimpose_res.append( k+1 )

    if subset_residues[k] in endpoints:
        working_endpoints.append( k+1 )

input_res_tag = make_tag( input_res )
working_loop_res_tag = make_tag( working_loop_res )



###############################################################
# Constraints.
###############################################################
align_file = mini_template_files[ 0 ]
print "ASSUMING alignment to: ", align_file

cst_file_tag = ''

if CUTOFF > 0:

    subset_res_tag = ''
    if len(working_superimpose_res) > 0:
        subset_res_tag = ' -subset'
        # Not working_superimpose_res -- superimpose.py will call pdbslice.py which uses original numbering.
        for m in superimpose_res: subset_res_tag += ' %d' % m

    command =  "superimpose.py %s %s  %s  > sup.pdb" % ( align_file, string.join( mini_template_files), subset_res_tag)
    print command
    system( command )
    system( "parse_NMR_models.py sup.pdb" )

    cst_depth = -10.0 / len( mini_template_files )
    count = 0
    cst_files = []
    for file in mini_template_files:
        count += 1
        if ( file.count('.pdb') ):
            cst_file =  file.replace( '.pdb', '_coordinate%3.1f.cst' % CUTOFF )
        else:
            cst_file =  file + '_coordinate%3.1f.cst' % CUTOFF

        tight_fade_tag = ''
        if tight_fade: tight_fade_tag = ' -tight_fade'

        command =  'generate_coordinate_CA_constraints.py  sup_%03d.pdb  -fade -stdev %3.1f -cst_res %s -anchor_res 1 -cst_depth  %8.3f %s > %s\n' % ( count+1, CUTOFF, working_loop_res_tag, cst_depth, tight_fade_tag, cst_file )
        print command
        system( command )
        assert( exists( cst_file ) )
        cst_files.append( cst_file )

    if len( input_cst_file ) > 0:
        command = 'rsync ../../%s .' % input_cst_file
        print command
        system( command )
        cst_files.append( basename( input_cst_file ) )

    cat_cst_file = 'cat_coordinate%3.1f.cst' % CUTOFF
    command = 'cat_cst.py '+string.join( cst_files )+' > '+cat_cst_file
    print command
    system( command )

    system( 'rm '+string.join( cst_files ) )
    cst_file_tag =  '-cst_file ' + cat_cst_file

    count = 0
    for file in mini_template_files:
        count += 1
        sup_file = file
        if (sup_file.count( '.pdb' ) == 0):  sup_file += '.pdb'
        sup_file = sup_file.replace('.pdb','.sup.pdb')
        system( 'renumber_pdb.py sup_%03d.pdb > %s' % (count+1, sup_file ) )


    system( 'rm sup_*.pdb')
    #system( 'rm sup.pdb')

###############################################################
# Prepare "README_SETUP"
###############################################################
fid = open( 'README_SETUP', 'w' )
fid.write( 'rm -rf CONDOR* SLAVE*  core*  REGION*  # region*\n' )

loop_force_Nsquared_tag = ''
if loop_force_Nsquared:
    loop_force_Nsquared_tag = ' -loop_force_Nsquared'

native_tag = ''
if len( native_file ) > 0:
    native_tag = ' -native %s ' % mini_native_file

virtual_res_tag = ''
if len( virtual_res ) > 0 :
    virtual_res_tag = ' -virtual_res '+make_rag( virtual_res )

superimpose_res_tag = ''
if len(working_superimpose_res) > 0:
    superimpose_res_tag = ' -superimpose_res'
    for m in working_superimpose_res: superimpose_res_tag += ' %d' % m

if ( working_loop_res[0] == 1 ) or ( working_loop_res[-1] == len( subset_residues ) ): # N-terminus or C-terminus
    start_pdb_tag = ' -start_pdb %s' % noloop_file
else: #internal loop
    start_pdb_tag = ' -loop_start_pdb %s' % noloop_file

cutpoints_open_tag = ''
if len(working_cutpoints_open) > 0 :
    cutpoints_open_tag = " -cutpoint_open %s" % make_tag( working_cutpoints_open )

command =  'grinder_dagman.py %s -loop_res %s -align_pdb %s -fasta %s  -nstruct 200 %s  %s %s %s %s -final_number %d' % ( start_pdb_tag, working_loop_res_tag, align_file, fasta_file, cst_file_tag, loop_force_Nsquared_tag, native_tag, superimpose_res_tag, cutpoints_open_tag, final_number )

if len( frag_files ) > 0:
    command += ' -frag_files'
    for file in frag_files:
        system( 'rsync ../../%s .' % file )
        command += ' %s' % basename( file )
else:
    command += ' -denovo 1'

if no_fixed_res: command += ' -no_fixed_res -loop_force_Nsquared'
if filter_native_big_bins: command += ' -filter_native_big_bins'

if len( working_endpoints ) > 0: command += ' -endpoints '+make_tag( working_endpoints )

fid.write( command  + '\n' )

fid.close()

###############################################################
# Prepare "README_SUB"
###############################################################
if ( NJOBS == 0 ):
    NJOBS = 20
    if len( working_loop_res ) > 20:
        NJOBS = 35
    if len( working_loop_res ) > 30:
        NJOBS = 50

#print NJOBS, len( working_loop_res), working_loop_res

fid = open( 'README_SUB', 'w' )
fid.write(' rm -rf blah.* \n' )
fid.write( 'bsub -W 96:0 -o blah.out -e blah.err SWA_pseudo_dagman_continuous.py  -j %d protein_build.dag  \n' % NJOBS )
fid.close()

chdir( CWD )

print
print "Setup directory: ",dir_name

