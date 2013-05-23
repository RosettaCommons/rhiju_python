#!/usr/bin/python

from sys import argv
from os import system,chdir,getcwd
from os.path import exists,abspath,basename,dirname
from get_sequence import get_annotated_sequence
import string
from glob import glob

OVERWRITE = 0
if argv.count('-overwrite'): OVERWRITE = 1

#dirnames = [ "baker_multi" ]
globdirs = glob( '*/' )
dirnames = []
for dir in globdirs:
    if dir not in ["GOOD/"]: dirnames.append( dirname(dir) )
dirnames.sort()

startpos = int( argv[1] )
endpos = int( argv[2] )


cwd = getcwd()
if OVERWRITE: system( "rm -rf GOOD/*TS? GOOD/*out*pdb" )

original_files = {}
refine_files = {}
for dirname in dirnames:
    copy_files = []

    if not exists( dirname ): continue
    chdir( dirname )

    orig_outfile = "region_FINAL.out"
    if not exists( orig_outfile ):
        orig_outfile = "region_%d_%d_sample.cluster.out" % (startpos,endpos)
    assert( exists(orig_outfile) )
    outfile = dirname+".out"

    if not exists( outfile ) or OVERWRITE:
        system( "ln -fs %s %s" % (orig_outfile,outfile) )
    if not exists( outfile+".1.pdb" ) or OVERWRITE:
        system( "rm -rf "+outfile+".*.pdb" )
        system( "extract_lowscore_decoys.py %s 10 -no_virtual" % outfile )

    ###################
    # Check for virtual
    ###################
    virtual_res = []

    if ( endpos > startpos ):
        print orig_outfile
        ( sequence, all_variant_types ) = get_annotated_sequence( orig_outfile )

        assert( len(sequence) == endpos - startpos + 1 )

        for i in range( len(sequence) ):
            for v in all_variant_types[ i ]:
                if v.count( 'Virtual_Protein_Residue' ):
                    virtual_res.append( startpos + i )

        if len( virtual_res ) > 0:
            outfile_with_virtual_res = outfile.replace('.out','.WITH_VIRTUAL_RES.out')
            if not exists( outfile_with_virtual_res ) or OVERWRITE:
                system( "ln -fs %s %s" % (orig_outfile,outfile_with_virtual_res) )
            if not exists( outfile_with_virtual_res+".1.pdb" ) or OVERWRITE:
                system( "rm -rf "+outfile_with_virtual_res+".*.pdb" )
                system( "extract_lowscore_decoys.py %s 9 " % outfile_with_virtual_res )


    ##############################
    # pdb slice any server models
    ##############################
    subset = []
    subset_tag = ''
    for m in range( startpos, endpos+1):
        if m not in virtual_res:
            subset.append( m )
            subset_tag += ' %d' % m

    globfiles = glob( "*TS?" )
    server_files = []
    for file in globfiles:
        if file[:4] != "mini":  server_files.append( file )

    for file in server_files:
        if ( exists( "mini_"+file ) and not OVERWRITE ): continue
        command = "pdbslice.py %s -subset %s mini_" % (file, subset_tag)
        print( command )
        system( command )

    globfiles = glob("mini_*TS?")
    globfiles.sort()
    for file in globfiles:
        if file.count("HHpred") and not dirname.count( "hhpred" ): continue
        copy_files.append( file )

    globfiles = glob("mini_mini_*TS?_H.pdb")
    globfiles.sort()
    for file in globfiles:
        if file.count("BAKER") and not dirname.count( "baker" ): continue
        copy_files.append( file )

    local_copy_files = []
    for file in copy_files: local_copy_files.append( file )
    original_files[ dirname ] = local_copy_files

    globfiles = glob( outfile+".*.pdb" )
    globfiles.sort()
    refine_files[ dirname ] = globfiles
    for file in globfiles: copy_files.append( file )

    #################################################
    # Copy into one directory for easy comparison
    #################################################
    copy_dir = "../GOOD"
    if not exists( copy_dir ): system( "mkdir -p "+copy_dir )

    for file in copy_files:
        if exists( copy_dir+"/"+file ) and not OVERWRITE: continue
        system( "rsync -avz "+file+" "+copy_dir )

    chdir( cwd )

##################################
# Carry out all_vs_all comparison
##################################
all_copy_files = []
for dirname in dirnames:
    for file in original_files[ dirname ]: all_copy_files.append( file )
    for file in refine_files[ dirname ]: all_copy_files.append( file )

chdir( "GOOD" )
fid = open( "list", 'w' )
for file in all_copy_files:
    fid.write( file+'\n' )
fid.close()

if not exists( "all_vs_all.txt" ) or OVERWRITE:
    command = "all_vs_all.py list > all_vs_all.txt"
    print command
    system( command )

gplot_file = "score_rms.gplot"
fid = open( gplot_file,'w' )
ps_file = "score_rms.ps"
fid.write( "plot '../%s/%s' u 26:2\n" % (dirnames[0], orig_outfile) )
for dirname in dirnames[1:]:
    fid.write( "replot '../%s/%s' u 26:2\n" % (dirname, orig_outfile) )
fid.write( "set term post color\n" )
fid.write( "set out '%s'\n" % ps_file )
fid.write( "replot\n"  )
fid.write( "set term x11\n" )
fid.write( "set out\n"  )
fid.write( "replot\n"  )
fid.close()

system( "gnuplot "+gplot_file )


################################
# Deviation plot
################################
gplot_file = "per_res.gplot"
fid = open( gplot_file,'w' )
ps_file = "per_res.ps"

pseudo_native = "mini_HHpredA_TS1"
if not exists( pseudo_native ):
    pseudo_native = "mini_mini_HHpredA_TS1_H.pdb"
    assert( exists( pseudo_native ) )

if not exists( 'per_res.txt' ) or OVERWRITE:
    command =  "~rhiju/python/superimpose.py "+ pseudo_native +" "+string.join( all_copy_files )+"  -R 2.0 -per_res > sup.pdb 2> per_res.txt"
    print command
    system( command )


all_copy_files = []
count = 0
ltcount = 0
for dirname in dirnames:
    ltcount += 1
    for file in original_files[ dirname ]:
        count += 1
        plot_tag = 'replot'
        if (count == 1): plot_tag = 'plot'
        fid.write("%s 'per_res.txt' u 1:%d w lines lt %d t '%s'\n" % (plot_tag, (1+count), ltcount, file) )
        print 1, file
    ltcount += 1
    for file in refine_files[ dirname ]:
        count += 1
        fid.write("%s 'per_res.txt' u 1:%d w lines lt %d t '%s'\n" % (plot_tag, (1+count), ltcount, file) )
        print 2, file

fid.write( "set yrange [0:16]\n" )
fid.write( "set title '%s'\n" % abspath( getcwd() ) )
fid.write( "set term post color\n" )
fid.write( "set out '%s'\n" % ps_file )
fid.write( "replot\n"  )
fid.write( "set term x11\n" )
fid.write( "set out\n"  )
fid.write( "replot\n"  )
fid.close()

system( "gnuplot "+gplot_file )


#########################################
# Need 3D viewer... pymol file.
#########################################
system( "parse_NMR_models.py sup.pdb" )

outfile = "TEST.pml"
fid = open( outfile, "w" )
print "  making: ",outfile

fid.write( "reinitialize\n")
fid.write( "cd %s\n" % getcwd() )
ltcount = 0
count = 0
colors_original = [ "red","blue","forest","orange"]
colors_refine = [ "magenta","cyan","green","wheat"]
for dirname in dirnames:
    for file in original_files[ dirname ]:
        count += 1
        fid.write("load sup_%03d.pdb, %s\n" % ( count+1, file ) )
        fid.write( "color %s, %s\n" % (colors_original[ltcount], file ) )
    for file in refine_files[ dirname ]:
        count += 1
        fid.write("load sup_%03d.pdb, %s\n" % ( count+1, file ) )
        fid.write( "color %s, %s\n" % (colors_refine[ltcount], file ) )
    ltcount += 1

fid.write( "hide everything, all\n" )
fid.write( "show cartoon\n" )

fid.write( "set cartoon_rect_length, 0.7\n" )
fid.write( "set cartoon_oval_length, 0.5\n" )
fid.write( "bg_color white\n" )
fid.close()
