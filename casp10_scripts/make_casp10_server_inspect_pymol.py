#!/usr/bin/python

from sys import argv
from os.path import basename,dirname,exists,expanduser,abspath
from os import system,chdir,getcwd
import string
from glob import glob
from parse_options import parse_options
from get_sequence import get_sequence

target_nums = map( lambda x:int(x), argv[ 1: ] )

servers_all = ["BAKER","Zhang-Server","HHpredA","HHpred2","ROBETTA","Midway","RaptorX"]

for target_num in target_nums:

    HOMEDIR = expanduser( '~' )
    workdir = HOMEDIR + '/projects/casp10/t%04d/' % target_num
    if not exists( workdir ): system( 'mkdir -p '+workdir)
    chdir( workdir )

    PYDIR = HOMEDIR + "/python"
    assert( exists( PYDIR) )

    #####################################
    # Download the server predictions
    #####################################
    if not ( exists( "server_predictions/" ) ):
        print "... Could not find: "+workdir+"/server_predictions"
        print "... Trying to download it from CASP website."
        system( dirname( abspath(argv[0]) )+'/get_casp10_server_predictions.py %d ' % target_num )
        if not ( exists( "server_predictions/" ) ):
            print "Problem!"
            exit( 0 )

    chdir( "server_predictions" )


    #####################################
    # Superimpose the "select" servers
    #####################################
    select_server_pdbs = []
    pdbs = {}
    servers = []
    for server in servers_all:
        globfiles = glob( server+"*" )
        globfiles.sort()
        if len( globfiles ) > 0: servers.append( server )
        pdbs[ server ] = globfiles

        for file in globfiles: select_server_pdbs.append( file )

    print "Superimposing... "
    print select_server_pdbs

    command = PYDIR+"/superimpose.py "+string.join( select_server_pdbs )+"  -per_res > sup.pdb 2> per_res.txt"
    print command
    system( command )

    command = PYDIR+"/parse_NMR_models.py sup.pdb"
    system( command )


    #################################################
    # pymol script to visualize superimposition
    #################################################
    print
    print "LOOK IN: %s" % ( getcwd() )
    print

    outfile = "TEST.pml"
    fid = open( outfile, "w" )
    print "  making: ",outfile

    fid.write( "reinitialize\n")
    count = 0
    for pdb in select_server_pdbs:
        count += 1
        fid.write("load sup_%03d.pdb, %s\n" % ( count, pdb ) )

    fid.write( "hide everything, all\n" )
    fid.write( "show cartoon\n" )

    fid.write( "set cartoon_rect_length, 0.7\n" )
    fid.write( "set cartoon_oval_length, 0.5\n" )
    fid.write( "bg_color white\n" )

    #################################################
    # Gnuplot to visualize Calpha deviations
    #################################################
    count = 0
    colors = [ "red","green","blue","magenta"]
    for server in servers:
        for pdb in pdbs[ server ]:
            fid.write( "color %s, %s\n" % (colors[count],pdb ) )
        count += 1

    fid.close()

    outfile = "make_CA_dev.gplot"
    fid = open( outfile, "w" )
    print "  making: ",outfile

    ps_file = "CA_dev.ps"
    count = 0
    server_count = 1
    for server in servers[0:]:
        server_count += 1

        for pdb in pdbs[ server ]:
            count += 1

            if count == 1: continue

            if count == 2: plot_tag = "plot"
            else: plot_tag = "replot"

            if pdb == pdbs[ server ][0]: lw = 3
            else: lw = 1

            fid.write("%s 'per_res.txt' u 1:%d w lines lt %d lw %d t '%s' \n" % (plot_tag,count, server_count, lw, pdb))

    fid.write( "set yrange [0:8]\n" )
    fid.write( "set xlabel 'resnum'\n " )
    fid.write( "set ylabel 'CA deviation'\n " )
    fid.write( "set title 'target number t%04d'\n" % target_num )

    fid.write( "set term post color\n" )
    fid.write( "set out '%s'\n" % ps_file )
    fid.write( "replot\n"  )
    fid.write( "set term x11\n" )
    fid.write( "set out\n"  )
    fid.write( "replot\n"  )

    fid.close()

    ps_file = "CA_dev.ps"
    print "  making: ", ps_file
    system( "gnuplot "+outfile )


    ###############################################################
    # Create beautiful secondary structure/contact-map file....
    ###############################################################

    def create_contact_map( servers_info_file, server_pdbs, native_file ):

        if not exists( servers_info_file ):

            server_pdbs_reorder = server_pdbs

            assert( server_pdbs_reorder.count( ref_file ) > 0 )

            pos = server_pdbs_reorder.index( ref_file )
            del( server_pdbs_reorder[ pos ] )
            server_pdbs_reorder.insert(0,ref_file)

            native_tag = ""
            if len( native_file ) > 0:
                native_tag = " -native %s" % native_file

            command = "%s/src/mini/bin/get_info.macosgccrelease -database ~/minirosetta_database -s %s -o %s %s -ignore_unrecognized_res -mute all " % \
                ( HOMEDIR, string.join(server_pdbs_reorder), servers_info_file, native_tag )
            print command
            system( command )

            print "  making: ", servers_info_file
        else:
            print " already exists: ", servers_info_file

        servers_info_ps = servers_info_file + ".ps"
        if not exists( servers_info_ps ):

            ss_tag = ''
            frag_tag = ''

            print "  making: ", servers_info_ps

            command = "%s/src/casp7/python/make_new_plot.py % s%s %s" % \
                ( HOMEDIR, ss_tag, frag_tag, servers_info_file )
            system( command )

        else:
            print " already exists: ", servers_info_ps

    if False:
        #ref_file = "HHpredA_TS1"
        ref_file = "Zhang-Server_TS1"
        if ref_file.replace( '_TS1','') not in servers: ref_file = "HHpred1_TS1"
        if ref_file.replace( '_TS1','')  not in servers: ref_file = "HHpred2_TS1"
        assert( exists( ref_file ) )

        all_server_pdbs = glob( "*_TS?" )
        ref_sequence = get_sequence( ref_file )
        filtered_server_pdbs = []
        for pdb in all_server_pdbs:
            if ( get_sequence( pdb ) == ref_sequence ):
                lines = open( pdb ).readlines()
                count = 0
                for line in lines:
                    if len( line) > 40  and line[12:16] == " N  ": count += 1
                if count == len( ref_sequence ): filtered_server_pdbs.append( pdb )
        print "After filtering, number of server models dropped from ", len( all_server_pdbs), " to ", len( filtered_server_pdbs )

        create_contact_map( "select_servers.info", select_server_pdbs, "" )
        create_contact_map( "all_servers.info", filtered_server_pdbs, "" )

    print
    print '--------------------------------------------'
    print "LOOK IN: %s" % ( getcwd() )
    print '--------------------------------------------'
