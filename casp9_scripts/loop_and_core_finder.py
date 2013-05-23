#!/usr/bin/python

from sys import argv,stdout,stderr
from os import system,getcwd
from os.path import abspath
import string
from math import sqrt
from parse_options import parse_options
from make_rhiju_color import make_rhiju_color
from make_tag import make_tag
from get_sequence import get_sequence

loop_res = parse_options( argv, 'loop_res', [-1] )
segment_res = parse_options( argv, 'segments', [-1] )
coord_cutoff = parse_options( argv, 'coord_cutoff', 2.0 )
final_number = parse_options( argv, 'final_number', 50 )
tight_fade = parse_options( argv, 'tight_fade', 0 )
cst_file = parse_options( argv, 'cst_file', "" )
frag_files = parse_options( argv, 'frag_files', [""] )
endpoints = parse_options( argv, 'endpoints', [-1] )
CUTOFF = parse_options( argv, 'core_cutoff', 0.5 )
njobs_input = parse_options( argv, 'j', 0 )
native_file = parse_options( argv, 'native', '' )
no_fixed_res = parse_options( argv, 'no_fixed_res', 0 )

#ONLY_USE_CLOSE_SUPERIMPOSE_RES = parse_options( argv, 'only_use_close_superimpose_res', 0 )
ONLY_USE_CLOSE_SUPERIMPOSE_RES = 1

#stderr.write( 'CUTOFF!!! %f\n'% CUTOFF )

#What are the files?
files = argv[1:]
if len( files ) == 0:
    stderr.write( 'Need to supply some starting templates\n')
    exit( 0 )


#######################
subset_res = []
subset_res_tag = ''
if len( segment_res ) > 0:
    for i in range( len(segment_res)/2):
        for j in range( segment_res[2*i], segment_res[2*i+1]+1 ):
            subset_res.append( j )
        segments_tag  = ' -segments'
    for i in segment_res: segments_tag += ' %d' % i
    subset_res_tag = ' -subset ' + make_tag( subset_res )
else:
    #sequence = get_sequence( files[ 0 ] )
    #nres = len( sequence )
    #subset_res = range(1, nres+1 )
    segments_tag = ''


# superimpose them.
command = "superimpose.py "+string.join( files ) + subset_res_tag + " -R 2.0 > sup.pdb "
#stderr.write( command + '\n' )
system( command )
# parse_models
system( "parse_NMR_models.py sup.pdb" )

def load_CA( file ):
    lines = open( file ).readlines()
    model_xyzs = []
    for line in lines:
        if len( line ) > 16 and line[12:16] == ' CA ':
            model_xyzs.append( [float(line[30:38]), float(line[38:46]), float(line[46:54])] )
    return model_xyzs

# Figure out potential core residues ... read in Calphas
CA_coords = []
for count in range( len( files ) ):
    sup_file = 'sup_%03d.pdb' % (count+1)
    CA_coords.append( load_CA( sup_file ) )
    #print "LENGTH ", len( CA_coords[count] )

# Calc rmsds (3D)
rmsds = []
nres =  len(CA_coords[0])


def get_dev2( coord1, coord2 ):
    dev2 = 0.0
    for k in range(3):
        dev = ( coord1[k] - coord2[k] )
        dev2 += dev*dev
    return dev2


for i in range( nres ):

    mean_coord = CA_coords[0][i]
    for j in range( 1,len(files) ):
        for k in range(3):  mean_coord[k] += CA_coords[j][i][k]
    for k in range(3): mean_coord[k] /= len( files )

    totdev2 = 0.0
    for j in range( len(files) ):
        totdev2 += get_dev2(  CA_coords[j][i], mean_coord )

    totdev2 /= len( files )
    rmsd = sqrt( totdev2 )

    rmsds.append( rmsd )

# output res vs. rmsds
outfile = 'rmsds.txt'
fid = open( outfile, 'w' )
for i in range( len( rmsds) ):
    fid.write( '%d %8.3f\n' % (i+1, rmsds[i]) )
fid.close()
stderr.write( 'Made '+ outfile+'\n')

##################################
outfile = 'per_res.gplot'
gplot_file = "per_res.gplot"
fid = open( gplot_file,'w' )
ps_file = "per_res.ps"

fid.write( "plot 'rmsds.txt' u 1:2 ps 1\n" )
fid.write( "replot 'rmsds.txt' u 1:2 w lines lt 1 \n" )
fid.write( "set yrange [0:8]\n" )
fid.write( "set xtics 10\n" )
fid.write( "set ytics 0.5\n" )
fid.write( "set grid \n" )
fid.write( "set title '%s'\n" % abspath( getcwd() ) )
fid.write( "set term post color\n" )
fid.write( "set out '%s'\n" % ps_file )
fid.write( "replot\n"  )
fid.write( "set term x11\n" )
fid.write( "set out\n"  )
fid.write( "replot\n"  )
fid.close()

system( "gnuplot "+gplot_file )


# apply cutoff to define core.
core_res = []
for i in range( len( rmsds) ):
    if ( rmsds[i] < CUTOFF ):
        core_res.append( i+1 )
        #stderr.write( '%d \n' % i )

# Do a screen for chainbreaks. Those are totally bad.
chainbreak_res = []
CA_CA_CUTOFF = 4.2
CA_CA_CUTOFF2 = CA_CA_CUTOFF * CA_CA_CUTOFF
for j in range( len(files ) ):
    CA = CA_coords[j]
    for i in range( nres-1 ):
        if ( i+1 not in chainbreak_res and get_dev2( CA[i], CA[i+1] )  > CA_CA_CUTOFF2):
            chainbreak_res.append( i+1 )

chainbreak_res.sort()
for i in chainbreak_res:
    stderr.write('Chainbreak residue! %d\n' % (i) )



hide_res = []
for i in range( 1, nres+1 ):
    if len( subset_res ) > 0 and i not in subset_res: hide_res.append( i )

superimpose_res = []
for i in core_res:
    if len( subset_res ) == 0  or i in subset_res:  superimpose_res.append( i )


#######################
#Setup pymol script.

make_rhiju_color( files, core_res, hide_res, subset_res )

files_including_native = []
for file in files: files_including_native.append( file )
if len( native_file ) > 0 : files_including_native.append( native_file )
if len( native_file ) == 0:  native_file = files[ 0 ] # pseudo_native

#######################
njobs = njobs_input
if len( loop_res ) > 0 :

    CUTOFF = 12.0
    CUTOFF2 = CUTOFF*CUTOFF

    tight_fade_tag = ''
    if tight_fade: tight_fade_tag = ' -tight_fade'

    fid = stdout

    for i in range( len( loop_res )/2 ):
        loop_start = loop_res[ 2*i ]
        loop_end = loop_res[ 2*i+1 ]

        stderr.write( 'Loop residues: %d %d [spacing: %d]\n' % (loop_start, loop_end, (loop_end-loop_start+1) ) )

        if ONLY_USE_CLOSE_SUPERIMPOSE_RES:
            superimpose_res = []
            for m in range( loop_start, loop_end+1 ):
                loop_coord = CA_coords[0][m-1]
                # assert( m not in core_res )
                for i in range( 1, nres+1 ):
                    if i not in subset_res: continue
                    if i not in core_res: continue
                    if i in superimpose_res: continue
                    if i in range( loop_start, loop_end+1): continue
                    check_coord = CA_coords[0][i-1]
                    dev2 = get_dev2( loop_coord, check_coord )
                    if ( dev2 < CUTOFF2 ):
                        superimpose_res.append( i )

        superimpose_res.sort()
        superimpose_res_tag = " -superimpose_res"
        for m in superimpose_res: superimpose_res_tag += " %d" % m

        if njobs_input == 0:
            njobs = 20
            if ( loop_end - loop_start + 1 ) > 20:                njobs = 35
            if ( loop_end - loop_start + 1 ) > 30:                njobs = 50

        for start_file in files_including_native:
            tag = start_file.replace('.pdb','')
            command = 'setup_consensus_loop_build_job.py  -tag %s   -s %s   -templates %s  %s  -native %s  -loop %d %d %s -coord_cutoff %6.2f %s -final_number %d  -j %d ' % ( tag, start_file, string.join( files ), segments_tag, native_file, loop_start, loop_end, superimpose_res_tag, coord_cutoff, tight_fade_tag, final_number, njobs )

            if len( cst_file ) > 0:  command += ' -cst_file ' + cst_file
            if len( frag_files ) > 0:  command += ' -frag_files ' + string.join( frag_files )
            if len( endpoints ) > 0: command += ' -endpoints ' + make_tag( endpoints )
            if no_fixed_res: command += ' -no_fixed_res'

            command += '\n'

            fid.write( command )
        fid.write( '\n' )

    fid.close()
    #stderr.write( '\nType:  source '+outfile+'\n')
