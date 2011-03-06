#!/usr/bin/python

from os import system
from sys import argv
import string
from os.path import exists

silent_files = argv[1:]

for silent_file in silent_files:

    assert( silent_file.find( ".out" ) )

    silent_file_out = silent_file.replace( ".out", ".cluster.out" )

    if exists( silent_file_out): system( "rm -rf " + silent_file_out )

    command = '~/src/mini_swa_rna/bin/rna_swa_test.macosgccrelease -in:file:silent %s -out:file:silent %s -algorithm cluster_old -database ~/minirosetta_database/ -silent_read_through_errors -cluster:radius 0.5 ' % (silent_file, silent_file_out )

    print command
    system( command )


