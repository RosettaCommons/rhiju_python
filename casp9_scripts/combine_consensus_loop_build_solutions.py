#!/usr/bin/python

import string
from sys import argv
from glob import glob
from os import system
from os.path import exists, basename
from parse_options import parse_options
from get_sequence import get_sequence
from make_tag import make_tag

newdir = parse_options( argv, "tag", 'consensus_combine')
if not exists( newdir ): system( 'mkdir -p '+newdir )

dirs = argv[1:]

loopdirs = []
for dir in dirs:
    if not exists( dir):
        print 'does not exist: ', dir
        exit( 0 )
    globdirs = glob( dir+'/loop*/' )
    for dir in globdirs:    loopdirs.append( dir )
loopdirs.sort()

def strip_tag( file, loop_tag ):
    return basename(file)[ (len(loop_tag)+1) : ].replace( '.pdb','')

###################################################################
# Figure out names of working directories,
#  loop residue numbering (in full pose and in each working pose)
###################################################################

template_names = []
all_loop_files = {}
all_loop_res = {}
all_working_loop_res = {}
cutpoint_open = []
for dir in loopdirs:

    loop_tag = string.split( dir, '/' )[1]
    globfiles = glob( dir+"/pdb/"+loop_tag+"*" )
    globfiles.sort()

    if ( len(template_names) == 0 ):
        for file in globfiles:
            template_names.append(  strip_tag( file, loop_tag)  )
        #print template_names

    for file in globfiles:
        #print "HEY!! ",  file, loop_tag
        template_name = strip_tag( file, loop_tag )
        if not template_name in template_names:
            print template_name, 'not in', template_names
            exit( 0 )

        if template_name not in all_loop_files.keys():
            all_loop_files[ template_name ] = []
            all_loop_res[ template_name ] = []
            all_working_loop_res[ template_name ] = []

        all_loop_files[ template_name ].append( file )

        check_loop_tag = loop_tag[4:]
        check_loop_tag = string.split( check_loop_tag, '_' )[0] # in case of _frags tag
        loop_bounds = string.split( check_loop_tag, '-' )
        loop_start = int( loop_bounds[0] )
        loop_end   = int( loop_bounds[1] )
        loop_res =  range( loop_start, loop_end+1)
        all_loop_res[ template_name ].append( loop_res )

        lines = string.join( open( dir + '/' + template_name.replace('.pdb','') + '/README_SETUP' ).readlines() )
        cols = string.split( lines )
        working_loop_res = parse_options( cols, "loop_res", [-1] )
        assert( len( working_loop_res ) > 0 )
        all_working_loop_res[ template_name ].append( working_loop_res )
        assert( len( working_loop_res) == len( loop_res ) )
        working_cutpoint_open = parse_options( cols, "cutpoint_open", [-1] )

        # Super-special case -- modeled an internal domain as if it was a terminus, with a free end.
        # Need to record this information and tell Rosetta exectuable so that it can create an appropriate fold_tree!
        if len( working_cutpoint_open ) > 0:
            # Map to full_pose.
            working_loop_res.sort()
            min_working_loop_res = working_loop_res[  0 ]
            max_working_loop_res = working_loop_res[ -1 ]
            for cut in working_cutpoint_open:
                if ( cut == min_working_loop_res-1 ):
                    cutpoint_open.append( loop_res[0] - 1 ) # cutpoint in "full pose" numbering
                elif ( cut == max_working_loop_res ):
                    cutpoint_open.append( loop_res[-1] ) # cutpoint in "full pose" numbering

cutpoint_open = set( cutpoint_open ) #uniquify cutpoint


##############################################################
# Check original starting template files are in each directory
#  and copy over working versions.
##############################################################
all_working_files = {}
for dir in dirs:
    for template_name in template_names:
        orig_file = dir+'/'+template_name
        if not exists( orig_file ):
            orig_file += '.pdb'
            if not exists( orig_file ):
                print 'Missing: ', orig_file
                exit( 0 )

        # Just need to do copying once
        if not( dir == dirs[0]): continue

        working_file = newdir + '/' + template_name + '_COMBINE.pdb'
        command = 'cp %s %s' % (orig_file, working_file )
        print command
        system( command )

        all_working_files[ template_name ] = working_file


first_sequence = ''
for template_name in template_names:

    working_file = all_working_files[ template_name ]

    sequence = get_sequence( working_file )
    if first_sequence == '': first_sequence = sequence
    assert( sequence == first_sequence )
    nres = len( sequence )

    loop_files = all_loop_files[ template_name ]
    loop_res = all_loop_res[ template_name ]
    working_loop_res = all_working_loop_res[ template_name ]

    temp_file = newdir+'/temp.pdb'

    #command =  'stepwise_protein_test.macosgccrelease  -database ~/minirosetta_database/  -combine_loops -s1 %s  -input_res1 %s    -s2 %s  -input_res2 %s  -slice_res2 %s  -o %s' % ( working_file, make_tag( range(1,nres+1) ), loop_files[i], make_tag(loop_res[i]), make_tag(working_loop_res[i]), temp_file)

    loop_bounds = []
    working_loop_bounds = []
    for i in range( len( loop_files ) ):
        loop_bounds.append( loop_res[i][0] )
        loop_bounds.append( loop_res[i][-1] )
        working_loop_bounds.append( working_loop_res[i][0] )
        working_loop_bounds.append( working_loop_res[i][-1] )


    cutpoint_tag = ''
    if ( len( cutpoint_open ) > 0 ):
        cutpoint_tag = ' -cutpoint_open'+make_tag( cutpoint_open )

    command =  'stepwise_protein_test.macosgccrelease  -database ~/minirosetta_database/  -combine_loops -start_pdb %s  -s %s  -loop_bounds %s  -working_loop_bounds %s %s -mute all  -o %s' % ( working_file, string.join( loop_files ), make_tag(loop_bounds), make_tag(working_loop_bounds), cutpoint_tag, temp_file)

    command += ' -pack_weights pack_no_hb_env_dep.wts'

    print command
    system( command )

    assert( exists( temp_file ) )

    system( 'mv %s %s' % (temp_file, working_file ) )

    print 'Made: ', working_file

