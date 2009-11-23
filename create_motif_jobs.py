#!/usr/bin/python

import string
from sys import argv,exit
from os import system

fasta_file = argv[1]
params_file = argv[2]

native_exists = 0
data_exists = 0
cst_exists = 0
torsions_exists = 0
for i in range( 3, len( argv ) ):
    if argv[i][-4:] == '.pdb':
        native_pdb_file = argv[i]
        native_exists = 1
    if argv[i][-5:] == '.data':
        data_file = argv[i]
        data_exists = 1
    if argv[i][-4:] == '.cst':
        cst_file = argv[i]
        cst_exists = 1
    if argv[i][-9:] == '.torsions':
        torsions_file = argv[i]
        torsions_exists = 1

# Read in files
lines = open( fasta_file ).readlines()
sequence = lines[1][:-1]
print sequence
numres = len( sequence )

44# Read in data information
data_info = []
if data_exists:
    backbone_burial_info = []
    lines = open( data_file ).readlines()
    for line in lines:
        if len( line ) > 6 and line[:6]=='EXPOSE':
            cols = string.split( line )
            for i in range( len(cols)/3 ):
                data_info.append( [int( cols[ 3*i+1] ), cols[3*i+2],cols[3*i+3]] )
        if len( line ) > 15 and line[:15]=='BACKBONE_BURIAL':
            cols = string.split( line )
            for i in range( 1,len(cols) ):
                backbone_burial_info.append( int( cols[i] ) )

pyrimidines = ['c','u']
purines = ['a','g']
cst_info = []
if cst_exists:
    lines = open( cst_file ).readlines()
    for line in lines:
        if len( line ) > 6 and line[0]!='[':
            cols = string.split( line )
            atom_name1 = cols[0]
            res1 = int( cols[1] )
            atom_name2 = cols[2]
            res2 = int( cols[3] )
            if ( sequence[ res1 - 1 ] in pyrimidines and atom_name1=='N1'):
                atom_name1 = 'N3'
                print 'correcting atom name for ', res1
            if ( sequence[ res2 - 1 ] in pyrimidines and atom_name2=='N1'):
                atom_name2 = 'N3'
                print 'correcting atom name for ', res2
            if ( sequence[ res1 - 1 ] in purines and atom_name1=='N3'):
                atom_name1 = 'N1'
                print 'correcting atom name for ', res1
            if ( sequence[ res2 - 1 ] in purines and atom_name2=='N3'):
                atom_name2 = 'N1'
                print 'correcting atom name for ', res2

            cst_info.append( [ atom_name1, res1, atom_name2, res2, string.join( cols[4:] )] )

pair_map = {}
all_pairs = []

complement = {'a':['u'], 'u':['a','g'], 'c':['g'], 'g':['c','u']};

stem_out_files = []
motif_out_files = []

# Parse out stems
lines = open( params_file ).readlines()
for line in lines:
    if line[:4] == 'STEM':
        cols = string.split( line )
        for i in range( len( cols )):
            if cols[i] == 'PAIR':
                #Offset to get to python numbering (starts with zero)
                res1 = int(cols[i+1])-1
                res2 = int(cols[i+2])-1
                pair_map[ res1 ] = res2
                pair_map[ res2 ] = res1
                all_pairs.append( [res1,res2] )
                assert ( sequence[res1] in complement[ sequence[res2] ] )
    else:
        try:
            cols = string.split( line[:-1] )
            res1 = int( cols[ 0 ] ) - 1
            res2 = int( cols[ 1 ] ) - 1
            pair_map[ res1 ] = res2
            pair_map[ res2 ] = res1
            all_pairs.append( [res1,res2] )
            assert ( sequence[res1] in complement[ sequence[res2] ] )
        except:
            continue

print pair_map

# Parse out stems
already_in_stem = {}
for i in range( numres ): already_in_stem[ i ] = 0

stems = []
stem_count = 0
for i in range( numres ):
    if pair_map.has_key( i ) and not already_in_stem[ i ]:  # In a base pair
        k = i
        stem_count += 1
        stem_res = []

        stem_res.append( [k, pair_map[k]] )
        already_in_stem[ k ] = 1
        already_in_stem[ pair_map[k] ] = 1

        # Can we extend in one direction?
        while( pair_map.has_key( k + 1 ) and  pair_map[ k+1 ] == pair_map[ k ] - 1 ):
            k += 1
            stem_res.append( [k, pair_map[k]] )
            already_in_stem[ k ] = 1
            already_in_stem[ pair_map[k] ] = 1

        # Do not allow single WC base pairs.
        if ( len( stem_res ) <2 ):
            print 'All stems must have length > 1 bp '
            print stem_res
            exit()
        stems.append( stem_res )

# Parse out motifs
already_in_motif = {}
for i in range( numres ): already_in_motif[ i ] = 0

motif_count = 0
motif_cutpoints = []
motif_res_maps = []
motif_stem_sets = []
motifs = []
for i in range( numres ):

    if ( not already_in_stem[ i ] and not already_in_motif[ i ]):

        motif_count += 1
        motif_res = []
        motif_stem_set = []
        cutpoints = []

        if ( i > 1 ):
            k = i-2

            motif_stem = []
            motif_stem.append( [ k, pair_map[k] ] )
            motif_stem.append( [ k+1, pair_map[k+1] ] )
            motif_stem_set.append( motif_stem )

            #Couple base pairs at end.
            motif_res.append( k )
            k += 1
            assert( already_in_stem[ k ] )
            motif_res.append( k )

            k = pair_map[ k ]
            assert( already_in_stem[ k ] )
            motif_res.append( k )
            k += 1
            assert( already_in_stem[ k ] )
            motif_res.append( k )

            cutpoints.append( k )

        k = i
        while ( k not in motif_res and k < numres):
            # Move forward to next stem:
            while ( k < numres and not already_in_stem[ k ] ):
                if (already_in_motif[ k ] ):
                    print 'Hey cant deal with pseudoknots!'
                    exit()
                motif_res.append( k )
                already_in_motif[ k ] = 1
                k += 1


            if k >= numres:
                cutpoints.append( k-1 )
                break

            if k in motif_res : break

            cutpoints.append( k+1 )

            motif_stem = []
            motif_stem.append( [ k, pair_map[k] ] )
            motif_stem.append( [ k+1, pair_map[k+1] ] )
            motif_stem_set.append( motif_stem )


            #Couple base pairs at end.
            motif_res.append( k )
            k += 1
            assert( already_in_stem[ k ] )
            motif_res.append( k )

            k = pair_map[ k ]
            assert( already_in_stem[ k ] )
            motif_res.append( k )
            k += 1
            assert( already_in_stem[ k ] )
            motif_res.append( k )

            # Next non-helical part..
            k += 1

        motif_res.sort()

        motif_res_map = {}
        for k in range( len( motif_res ) ):
            motif_res_map[ motif_res[k] ] = k

        motifs.append( motif_res )
        motif_stem_sets.append( motif_stem_set )
        motif_cutpoints.append( cutpoints )
        motif_res_maps.append( motif_res_map )
        print 'CUTPOINTS ', cutpoints

#print motifs
#print motif_stem_sets
#print motif_cutpoints

# Output stem definition jobs
fid_README_STEMS = open( 'README_STEMS','w')
for i in range( stem_count ):

    # Fasta
    tag = 'stem%d_%s' % (i+1, fasta_file)
    fid = open( tag , 'w' )
    fid.write( '>'+tag+'\n')

    stem_res = stems[i]
    stem_length = len( stem_res )

    for k in range( stem_length ):
        fid.write( sequence[stem_res[k][0]] )
    for k in range( stem_length ):
        fid.write( sequence[stem_res[stem_length-k-1][1]] )
    fid.write('\n')
    fid.close()

    # pdb_file
    if native_exists:
        command = 'pdbslice.py  %s -segments %d %d %d %d stem%d_' %(
            native_pdb_file,
            stem_res[0][0]+1,
            stem_res[-1][0]+1,
            stem_res[-1][-1]+1,
            stem_res[0][-1]+1,
            i+1 )
        print command
        system( command )

    outfile = 'stem%d_%s.out' % (i+1, fasta_file.replace('.fasta',''))
    stem_out_files.append( outfile )
    command = 'rna_assemble_test.macosgccrelease -database  /Users/rhiju/minirosetta_database -nstruct 1 -build_helix_test -fasta %s -out:file:silent %s' % (tag, outfile)
    fid_README_STEMS.write(command+'\n')

fid_README_STEMS.close()


# Output motif jobs
fid_README_MOTIFS = open( 'README_MOTIFS','w')
for i in range( motif_count ):

    # Fasta
    motif_fasta_file = 'motif%d_%s' % (i+1, fasta_file)
    fid = open( motif_fasta_file , 'w' )
    fid.write( '>'+motif_fasta_file+'\n')

    motif_res = motifs[i]
    motif_length = len( motif_res )

    for k in range( motif_length ):
        fid.write( sequence[motif_res[k]] )
    fid.write('\n')
    fid.close()

    # params file
    motif_stem_set = motif_stem_sets[ i ]
    motif_res_map = motif_res_maps[ i ]
    motif_cutpoint = motif_cutpoints[ i ]

    motif_params_file = 'motif%d_%s.params' % (i+1, fasta_file.replace('.fasta',''))
    fid = open( motif_params_file , 'w' )

    for k in range( len(motif_stem_set) ):
        motif_stem = motif_stem_set[ k ]
        fid.write( 'STEM   PAIR %d %d W W A   PAIR %d %d W W A \n' % \
                       ( motif_res_map[ motif_stem[0][0] ]+1,
                         motif_res_map[ motif_stem[0][1] ]+1,
                         motif_res_map[ motif_stem[1][0] ]+1,
                         motif_res_map[ motif_stem[1][1] ]+1 ) )

    motif_cutpoint.sort()
    if ( len( motif_cutpoint ) > 1 ):
        fid.write( 'CUTPOINT_OPEN ' )
        for k in range( len( motif_cutpoint ) ):
            if motif_res_map[ motif_cutpoint[k] ] < (len( motif_res )-1):
                fid.write( ' %d' % (motif_res_map[ motif_cutpoint[k] ]+1) )
    fid.write('\n')
    fid.close()

    # pdb_file
    native_tag = ''
    if native_exists:
        command = 'pdbslice.py  %s -subset ' %  native_pdb_file
        for k in range( motif_length ): command += ' %d' % (motif_res[k]+1)
        command += ' motif%d_' % (i+1)
        print command
        system( command )
        native_tag = '-native motif%d_%s' % (i+1, native_pdb_file )

    if data_exists:
        motif_data_file = 'motif%d_%s' % ( i+1, data_file )
        fid_data = open( motif_data_file, 'w' )
        fid_data.write( 'EXPOSE' )
        for data in data_info:
            if data[0]-1 in motif_res_map.keys():
                fid_data.write( '   %d %s %s ' % (motif_res_map[data[0]-1]+1,data[1],data[2]) )
        fid_data.write('\n')

        if len( backbone_burial_info ) > 0:
            fid_data.write( 'BACKBONE_BURIAL ' )
            for k in backbone_burial_info:
                if k-1 in motif_res_map.keys():
                    fid_data.write( ' %d' % (motif_res_map[ k-1 ] + 1) )
            fid_data.write( '\n' )
        fid_data.close()

    cst_found = 0;

    if cst_exists:
        motif_cst_file = 'motif%d_%s' % ( i+1, cst_file )
        fid_cst = open( motif_cst_file, 'w' )
        fid_cst.write( '[ atompairs ]\n' )
        for cst in cst_info:
            if cst[1]-1 in motif_res_map.keys() and cst[3]-1 in motif_res_map.keys():
                fid_cst.write( '%s %d %s %d %s\n' % (cst[0], motif_res_map[cst[1]-1]+1,cst[2],motif_res_map[cst[3]-1]+1,cst[4]) )
                cst_found = 1
        fid_cst.close()

    motif_out_file = motif_params_file.replace( '.params','.out')
    motif_out_files.append( motif_out_file )
    NSTRUCT = 40
    command = 'rna_denovo.macosgccrelease -database  /Users/rhiju/minirosetta_database %s -fasta %s -params_file %s -nstruct %d -out::file::silent %s -cycles 5000 -mute all -vary_geometry' % \
        ( native_tag, motif_fasta_file, motif_params_file, NSTRUCT, motif_out_file )

    if data_exists: command += ' -data_file %s ' % motif_data_file
    if cst_exists and cst_found: command += ' -cst_file %s ' % motif_cst_file
    if torsions_exists: command += ' -vall_torsions %s ' % torsions_file
    fid_README_MOTIFS.write( command+'\n' )


fid_README_MOTIFS.close()


# Output assembly job
#Where are the jumps and chainbreaks?
jumps = []
cutpoints  = []
for i in range( motif_count ):

    motif_stem_set = motif_stem_sets[ i ]

    for k in range( len(motif_stem_set) ):
        motif_stem = motif_stem_set[ k ]
        possible_cutpoints =  [ motif_stem[ 0 ][ 0 ], motif_stem[ 1 ][ 1 ] ]
        possible_cutpoints.sort()
        if ( possible_cutpoints[0] not in cutpoints):
            cutpoints.append( possible_cutpoints[ 0 ] )


params_file = fasta_file.replace('.fasta','_assemble.params' )
fid = open( params_file, 'w')
fid.write( 'CUTPOINT_OPEN ' )
cutpoints.sort()
for cutpoint in cutpoints:
    fid.write( ' %d'  % (cutpoint+1) )
fid.write( '\n' )

for cutpoint in cutpoints:
    fid.write( 'OBLIGATE   PAIR %d %d W W A\n' % (cutpoint+1, pair_map[cutpoint]+1) )

for i in range( stem_count ):
    stem_res = stems[i]
    fid.write( 'STEM ')
    for k in range( len( stem_res )):
        fid.write( ' PAIR %d %d W W A ' % \
                       ( stem_res[k][ 0 ]+1, stem_res[k][ 1 ]+1 ) )
    fid.write('\n')

fid.close()
########
assemble_cst_file = params_file.replace('.params','.cst')
if cst_exists:
    assemble_cst_file = cst_file.replace('.cst','_assemble.cst')
fid = open( assemble_cst_file,'w')
fid.write('[ atompairs ]\n')
for cutpoint in cutpoints:
    fid.write( 'O3*  %d  P     %d  HARMONIC  1.619  2.0\n' % \
                   ( cutpoint+1, cutpoint+2 ) )
if cst_exists:
    for cst in cst_info:
            fid.write( '%s %d %s %d %s \n' % (cst[0], cst[1], cst[2], cst[3], cst[4]) )


fid.close()


#########
fid = open( 'README_ASSEMBLE', 'w' )

native_tag = ''
if native_exists: native_tag = '-native '+native_pdb_file

outfile = params_file.replace( '.params','.out' )
command = 'rna_denovo.macosgccrelease -database  /Users/rhiju/minirosetta_database %s -fasta %s -in:file:silent_struct_type binary_rna  -cycles 10000 -nstruct 200 -out:file:silent %s -params_file %s -cst_file %s -close_loops -in:file:silent ' % \
( native_tag, fasta_file, outfile, params_file, assemble_cst_file )

for stem_out_file in stem_out_files:
    command += ' '+stem_out_file
for motif_out_file in motif_out_files:
    command += ' '+motif_out_file
if torsions_exists: command += ' -vall_torsions %s ' % torsions_file
if data_exists:
    command += ' -data_file '+data_file

fid.write( command+'\n')
fid.close()


