#!/usr/bin/python

from sys import stdout,stderr
from os.path import abspath
from os import system,getcwd
import string

def make_rhiju_color( files, core_res = [], hide_res = [], subset_res = []):

    subset_res_tag = ""
    if len( subset_res ) > 0:
        subset_res_tag = " -subset"
        for m in subset_res: subset_res_tag += " %d" % m

    # superimpose them.
    system( "superimpose.py -R 2.0 "+string.join( files )+" "+subset_res_tag+" > sup.pdb " )
    system( "parse_NMR_models.py sup.pdb " )

    #Setup pymol script.
    outfile = 'TEST.pml'
    fid = open( outfile, 'w' )

    fid.write( 'reinitialize\n')
    fid.write( 'cd ' + abspath( getcwd() ) + '\n' )

    for count in range( len( files ) ):
        sup_file = 'sup_%03d.pdb' % (count+1)
        fid.write( 'load %s, %s\n' % (sup_file, files[count]) )

    # Color in my favorite coloring
    fid.write( 'hide everything\n' )
    fid.write( 'show ribbon\n' )
    fid.write( 'set ribbon_sampling,1\n' )

    fid.write( 'bg_color white\n')
    fid.write( 'color gray\n')
    fid.write( 'color blue, resn his+lys+arg\n')
    fid.write( 'color red, resn asp+glu\n')
    fid.write( 'color green, resn gln+asn+thr+ser\n')
    fid.write( 'color orange, resn gly\n')
    fid.write( 'color purple, resn cys\n')

    fid.write( 'show spheres, resn pro and name CA\n')

    fid.write( 'set cartoon_rect_length,1.0\n')
    fid.write( 'set cartoon_oval_length,1.0\n')

    # cartoon for core residues.
    #  or if loops are defined, cartoon for all non-loop residues
    #  and white for core.

    fid.close()
    stderr.write( 'Made ' + outfile + '\n' )


    outfile = 'TEST.script'
    fid = open( outfile, 'w' )
    fid.write( 'load sup.pdb\n')

    fid.write('wireframe off\n')
    fid.write('select all\n')
    fid.write('backbone 70\n')
    fid.write('select hydrophobic\n')
    fid.write('color gray\n')
    fid.write('select polar\n')
    fid.write('color green\n')
    fid.write('select positive\n')
    fid.write('color blue\n')
    fid.write('select negative\n')
    fid.write('color red\n')
    fid.write('select gly\n')
    fid.write('color gold\n')
    fid.write('select cys\n')
    fid.write('color purple\n')

    fid.write('select *.CA and pro\n')
    fid.write('spacefill 200\n')

    fid.write('select *.CA and 1\n')
    fid.write('spacefill 200\n')
    fid.write('color blue\n')


    for m in hide_res:
        fid.write( 'select %d\n' % m)
        fid.write( 'color cyan\n')

    fid.write('set vectps on\n')
    fid.write('set specular on\n')
    fid.write('restrict not hydrogen\n')

    for m in core_res:
        fid.write( 'select %d\n' % m)
        fid.write( 'color white\n')
    fid.write('select all\n')

    fid.write('\n')
    fid.close()
    stderr.write( 'Made ' + outfile + '\n' )

    ####################################################################
