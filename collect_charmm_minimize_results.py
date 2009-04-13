#!/usr/bin/python

from sys import argv,exit


from glob import glob
from os.path import basename,exists,dirname,abspath
from os import popen,system
import string

outdirs = argv[1:]

native_supplied = 0
if (outdirs[0][-4:] == '.pdb' ): # Native specified
    native_pdb = outdirs[ 0 ]
    del( outdirs[ 0 ] )
    native_supplied = 1


# Do not use outfile!!

for outdir in outdirs:

    assert( outdir.count('_OUT' ) )


    native_exists = native_supplied
    if not native_exists:
        # Still look for native in rhiju's directory.
        pos = outdir.index( 'chunk')
        rna_name = outdir[pos:(pos+13)]
        native_pdb = '/home/rhiju/projects/rna_new_benchmark/bench_final/%s_RNA.pdb' % rna_name
        if exists( native_pdb ):
            native_exists = 1


    #############################################
    # Go through scorefile from tinker minimization

    globfiles = glob( outdir+'/S*/*.scores' )
    globfiles.sort()

    all_scores = {}
    score_terms = []

    obligate_score_terms =  [ 'total' ]
    for score_term in obligate_score_terms:
        score_terms.append( score_term )
        all_scores[ score_term ] = {}

    files_to_calc_rms = []
    for file in globfiles:

        for score_term in obligate_score_terms: all_scores[score_term][file] = 0.0 #Always have one of these terms.
        lines = open( file ).readlines()
        for line in lines:
            cols = string.split( line )
            try:
                score = float( cols[1] )
                score_term = cols[0]
                if score_term not in score_terms:
                    score_terms.append( score_term )
                    all_scores[ score_term ] = {}
                all_scores[ score_term ][ file ] = score
            except:
                continue


        if native_exists:
            pdbfile = file.replace('.scores','')
            rmsfile = pdbfile + '.rms.txt'
            if  exists( pdbfile ) and not exists( rmsfile ):
                files_to_calc_rms.append( pdbfile )

            pdbfile = file.replace('.scores','.min_pdb')
            rmsfile = pdbfile + '.rms.txt'
            if  exists( pdbfile ) and not exists( rmsfile ):
                files_to_calc_rms.append( pdbfile )


    ##########################################################################
    # superimpose any files that haven't been supermposed
    num_files = len( files_to_calc_rms)
    N_DECOYS = 500
    if native_exists and num_files > 0 :
        #This is stupid, there's a limit to how many can go into a command line
        for x in range(  int( num_files/N_DECOYS ) + 1 ):
            command = '~rhiju/python/charmm_superimpose.py '+native_pdb+' '+string.join( \
                files_to_calc_rms[ (N_DECOYS*x): (N_DECOYS*x+N_DECOYS)] )
            print( command )
            system( command )


    ##########################################################################
    # read out rms's from saved files.
    rms_vals = {}
    rms_vals_init = {}


    for file in globfiles:

        rmsfile = file.replace('.scores','') +'.min_pdb.rms.txt'

        if not exists( rmsfile): continue
        rmslines = popen( 'tail -n 1 '+rmsfile ).readlines()
        if len( rmslines ) < 1: continue
        rmsline = rmslines[-1]
        if len( string.split(rmsline) ) < 1: continue
        rms_vals[ file ] = float( string.split( rmsline )[-1] )

        rmsfile = file.replace('.scores','') +'.rms.txt'
        if not exists( rmsfile): continue
        rmslines = popen( 'tail -n 1 '+rmsfile ).readlines()
        if len( rmslines ) < 1: continue
        rmsline = rmslines[-1]
        if len( string.split(rmsline) ) < 1: continue
        rms_vals_init[ file ] = float( string.split( rmsline )[-1] )
        #print rmsfile



    ######################################################
    # Output unified scorefile.
    new_scorefile = basename( abspath(outdir) ).replace( '_OUT', '_CHARMM_minimize.sc' )
    print ' Creating '+new_scorefile
    fid = open( new_scorefile ,'w')

    fid.write( 'SCORE: ' )
    for score_term in score_terms:
        fid.write( ' %10s' % score_term[:10] )

    if native_exists:
        fid.write( ' %8s' % 'rms_init' )
        fid.write( ' %8s' % 'rms' )

    fid.write( ' description\n' )

    for file in globfiles:

        if not file in rms_vals.keys(): continue
        #if not file in rms_vals_init.keys(): continue


        score_term_present = 1
        for score_term in score_terms:
            if not file in all_scores[score_term].keys():
                score_term_present = 0
                break
        if not score_term_present: continue

        fid.write( 'SCORE: ' )

        for score_term in score_terms:
            fid.write( ' %10.4f' % all_scores[ score_term ][file ] )

        if native_exists:
            #tag = basename( file ).replace('min_','').replace('.sc','').replace('minimize_','')
            #fid.write( ' %8.4f' % rms_vals_init[ file ] )
            fid.write( ' %8.4f' % rms_vals[ file ] )

        fid.write( ' '+basename(file).replace('.scores','') + '\n' )

    fid.close()
