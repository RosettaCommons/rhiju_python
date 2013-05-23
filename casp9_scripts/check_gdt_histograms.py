#!/usr/bin/python

from sys  import argv
from os import system
from os.path import abspath


gdt_file = argv[1]

gdt_servers_file = gdt_file.replace( '.txt','_servers.txt')
system( 'grep _TS %s > %s' % ( gdt_file, gdt_servers_file ) )

gdt_submissions_file = gdt_file.replace( '.txt','_submissions.txt')
system( 'grep casp %s > %s' % ( gdt_file, gdt_submissions_file ) )

gnuplot_file =  'make_gdt_histograms.gplot'
fid = open( gnuplot_file, 'w' )

fid.write( 'plot "< histo.py %s 4 0.005 0 1.0" u 1:2 w lines\n' % gdt_servers_file )
fid.write( 'replot "< histo.py %s 4 0.005 0 1.0" u 1:2 w lines lt 3\n' % gdt_submissions_file )
fid.write( 'set title "%s"\n' % abspath( gdt_file) )

fid.write( 'set term post color\n' )
fid.write( 'set out "gdt_histograms.ps" \n')

fid.write( 'replot\n' )
fid.write( 'set term x11\n' )
fid.write( 'set out\n' )

fid.close()

system( 'gnuplot %s ' % gnuplot_file )


