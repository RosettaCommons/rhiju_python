#!/usr/bin/python

# Run this from inside casp7/t0523  etc...

from sys import argv
from os.path import basename,dirname,exists,expanduser
from os import system,chdir,getcwd, popen
import string
from glob import glob
from parse_options import parse_options
from get_sequence import get_sequence

target_num = int( argv[ 1 ] )

HOMEDIR = expanduser( '~' )
PYDIR = HOMEDIR + "/python"
assert( exists( PYDIR) )
IDEALIZE_EXE = HOMEDIR+"/src/mini/bin/idealize.macosgccrelease"
assert( exists( IDEALIZE_EXE ) )
DB = HOMEDIR+"/minirosetta_database"
assert( exists( DB ) )

workdir = HOMEDIR + '/projects/casp9/t%04d/' % target_num
assert( exists( workdir ) )
chdir( workdir )

pathway = parse_options( argv, "pathway", [-1] )
pick_endpoints = parse_options( argv, "pick_endpoints", 0 )
optimal_length = parse_options( argv, "optimal_length", 10 )
if ( len( pathway ) == 0 and not pick_endpoints ):
    print "Must specify -pathway or -pick_endpoints"
    exit( 0 )

endpoint_spacing = parse_options( argv, "endpoint_spacing", 0 )

dirtag = parse_options( argv, "tag", "round1" )

cst_coords = parse_options( argv, "cst_coords",0)
cst_HB = parse_options( argv, "cst_HB",0)
cst_CA_CA = parse_options( argv, "cst_CA_CA",0)

if ( not cst_coords and not cst_HB and not cst_CA_CA ):
    print " Need to specify -cst_coords, -cst_HB, or -cst_CA_CA. "
    print " In an early version of this script, just used coords by default, "
    print "  but that is not recommended..."
    exit( 0 )


assert( exists( "rosetta_frags/" ) )
assert( exists( "server_predictions/" ) )

################
# Following aren't easy to get from my "parse_options" class, because they could
# be a bunch of "-loose_res " options after the "-cst_args" tag.
important_args = ["-cst_args","-grind_args"]
cst_args = ""
grind_args = ""

native_pdb = parse_options( argv, "native", "HHpredA_TS1" )

if argv.count( "-cst_args" ) > 0 :
    pos = argv.index( "-cst_args" )
    del( argv[ pos ] )
    while pos < len(argv) and argv[ pos ] not in important_args:
        cst_args += ' '+argv[ pos ]
        del( argv[ pos ] )

if argv.count( "-grind_args" ) > 0 :
    pos = argv.index( "-grind_args" )
    del( argv[ pos ] )
    while pos < len(argv) and argv[ pos ] not in important_args:
        grind_args += ' '+argv[ pos ]
        del( argv[ pos ] )

#print "GRIND_ARGS ", grind_args
#print "CST_ARGS", cst_args
#exit( 0 )


################
FINAL_NUMBER = 50
NSTRUCT = 50

def setup_job( jobdir, template_pdbs, FRAGS = 0 ):

    cwd = getcwd()

    newdir = jobdir

    system( "mkdir -p " +newdir )

    chdir( "server_predictions")
    globfiles = []
    for pdb in string.split(template_pdbs):
        someglobfiles = glob( pdb )
        someglobfiles.sort()
        for globfile in someglobfiles: globfiles.append( globfile )

    assert( len( globfiles ) > 0 )

    chdir( "../"+newdir )

    ###################################
    # (pseudo)native and fasta
    ###################################

    file = "../../server_predictions/"+native_pdb
    assert( exists( file ) )
    system( "rsync "+file+" . " )

    fasta_file = "t%3d_.fasta" % ( target_num )
    file = "../../rosetta_frags/"+fasta_file
    if not exists( file ):
        fasta_file = "T%04d.fasta" % ( target_num )
        file = "../../rosetta_frags/"+fasta_file


    assert( exists( file ) )
    system( "rsync "+file+" . " )

    print "rsynced native and fasta"

    sequence_lines = open( fasta_file  ).readlines()[1:]
    sequence = string.join(  map( lambda x : x[:-1], sequence_lines) ,  '' )

    idealized_template_pdbs = []
    globfiles_good = []
    for pdb in globfiles:
        if not exists( pdb ): system( "rsync ../../server_predictions/"+pdb+" . " )
        pdb_sequence = get_sequence( pdb )
        if not ( pdb_sequence  == sequence ): continue
        globfiles_good.append( pdb )
        idealized_file = pdb+"_ideal.pdb"
        if not exists( idealized_file ):
            outfile = pdb+"_0001.pdb"
            system( "rm -rf "+outfile )

            command = "%s -database %s -s %s -fast " % (IDEALIZE_EXE, DB, pdb )
            system( command )
            system( "mv "+outfile+" "+idealized_file )
            assert( exists( idealized_file ) )

        idealized_template_pdbs.append( idealized_file )


    endpoint_tag = ''
    if pick_endpoints:
        #assert( len( pathway) == 0 )
        assert( not cst_coords )
        psipred_files = glob( '../../rosetta_frags/*psipred_ss2' )
        assert( len( psipred_files ) == 1 )
        psipred_ss2_file = psipred_files[ 0 ]
        endpoint_tag = popen( 'endpoint_picker.py -optimal_length %s %s ' % (optimal_length, psipred_ss2_file)  ).readlines()[-1][:-1]

    if len( pathway ) > 0:
        pathway_file = "t%04d_guess.pathway" % target_num
        command = PYDIR+"/generate_pathway.py "
        for m in pathway: command += ' %d' % m
        command += " > "+pathway_file
        system( command )
        assert( exists( pathway_file ) )

        anchor_res = int(  string.split( open( pathway_file ).readlines()[0] )[1] )

    ###################################
    # constraint setup
    ###################################
    first_pdb = globfiles_good[0]
    print "Basing constraints off: ", first_pdb
    cst_pdb = first_pdb
    cst_file = first_pdb+".cst"

    separate_cst_files = []

    if cst_HB:
        cst_HB_file = first_pdb+".HB.cst"
        separate_cst_files.append( cst_HB_file )
        print "#################################################################"
        print "WARNING: by default HB csts ignore cst_res, loose_res, no_cst_res"
        print "#################################################################"
        command = PYDIR + "/generate_backbone_HB_constraints.py %s -fade  -stdev 1.0  > %s" % (first_pdb, cst_HB_file )
        print command
        system( command )

    if cst_coords:
        cst_coord_file = first_pdb+".coord.cst"
        separate_cst_files.append( cst_coord_file )
        command = PYDIR + "/generate_coordinate_CA_constraints.py %s %s -anchor_res %d -fade  > %s" % (first_pdb, cst_args, anchor_res, cst_coord_file )
        print command
        system( command )

    if cst_CA_CA:
        cst_CA_CA_file = first_pdb+".CA_CA.cst"
        separate_cst_files.append( cst_CA_CA_file )
        command = PYDIR + "/generate_CA_constraints.py %s %s -fade > %s" % (first_pdb, cst_args, cst_CA_CA_file )
        print command
        system( command )

    if cst_CA_CA and cst_coords:
        print "cannot specify both -cst_CA_CA and -cst_HB. If you want both, fix this script to properly concatenate cst files to combine [atompairs] blocks"
        exit( 0 )

    system( "cat %s > %s " % (string.join(separate_cst_files), cst_file ) )


    fid = open( "README_SETUP", 'w' )
    fid.write( "rm -rf REGION* *~ CONDOR core.* SLAVE*  \n")

    ###################################
    # any frags
    ###################################
    frag_tag = ""
    if FRAGS:
        globfiles = glob( '../../rosetta_frags/*_05.200_v1_3*' )
        assert( len( globfiles ) > 1 )
        system( "rsync " + string.join(globfiles) + " . " )
        globfiles.sort()

        frag_tag = " -frag_files "+string.join(  map( lambda x:basename(x), globfiles ) )

    print "copied in frags"
    print "CWD: ", getcwd()

    ##################################################
    # endpoint setup -- sparse sampling for speed.
    ##################################################
    if ( len( pathway ) > 0 and endpoint_spacing > 0 ):

        endpoints = [ pathway[0] ]
        start_path = pathway[ 0]
        end_path = pathway[ 0 ]

        for m in range( len(pathway) - 1 ):
            current_res = pathway[ m+1 ]
            if ( current_res > pathway[ m ] ):
                previous_res = end_path
                end_path = current_res
            else:
                previous_res = start_path
                start_path = current_res

            spacing = abs( current_res - previous_res )
            sign = 1
            if (current_res < previous_res): sign = -1
            num_stepping_stones = spacing / endpoint_spacing

            #print "SPACING",spacing,"  NUM_STEPPING_STONE",num_stepping_stones

            for k in range( num_stepping_stones ):
                stepping_stone = previous_res + sign * endpoint_spacing * (k+1)
                #print "stepping stone: ", stepping_stone
                endpoints.append(  stepping_stone )
            if current_res not in endpoints: endpoints.append( current_res )


        endpoints.sort()
        endpoint_tag = ' -endpoints'

        print 'ENDPOINTS: ', endpoints
        print 'number of endpoints: ', len( endpoints )

        for m in endpoints: endpoint_tag += ' %d' % m

    pathway_file_tag = ''
    if len( pathway ) > 0: pathway_file_tag = ' -pathway_file %s' % pathway_file

    fid.write( "grinder_dagman.py  -native %s -fasta %s -cluster_radius 0.5 -final_number %d -nstruct %d  %s  -cst_file %s -align_pdb %s -s  %s  %s %s %s\n" \
                   % (native_pdb, fasta_file, FINAL_NUMBER, NSTRUCT, pathway_file_tag, cst_file, cst_pdb, string.join( idealized_template_pdbs ), frag_tag, endpoint_tag, grind_args ) )

    fid.close()

    fid = open( "README_SUB", "w" )
    fid.write( "rm -rf blah.*\n" )

    N_JOBS = 10 * len( idealized_template_pdbs )
    if FRAGS: N_JOBS += 50

    fid.write( "bsub -W 96:0 -o blah.out -e blah.err SWA_pseudo_dagman_continuous.py  -j %d  protein_build.dag\n" % N_JOBS)
    fid.close()

    chdir(cwd)

    print "****************************************************"
    print "Made directory ", newdir
    print "****************************************************"
    print


###################################################
###################################################
###################################################


#######################
# Original 3-way split
#######################
#setup_job( dirtag+"/hhpred_and_frags", "HHpredA*",  1 )
#setup_job( dirtag+"/baker_multi", "BAKER*",  0 )
#setup_job( dirtag+"/zhang_midway_multi", "Zhan* Midw*",  0 )

#######################
# Grind everything
#######################
#setup_job( target_num, dirtag+"/all_and_frags", "HHpredA* BAKER* Zhan* Midw*",  cst_args, grind_args, 1 )

# Single server model refinement.
setup_job( dirtag+"/hhpred_and_frags", "HHpredA_TS1",  1 )
setup_job( dirtag+"/multi1" , "Zhan* Midway*"  ,  0 )
setup_job( dirtag+"/multi2" , "PconsR* MULTICOM-REF*"  ,  0 )
setup_job( dirtag+"/baker" , "BAKER*"  ,  0 )
