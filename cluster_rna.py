#!/usr/bin/python

from sys import argv,exit
from os import system,popen
from os.path import exists,basename
import string

verbose = 0
if argv.count( '-verbose'):
    pos = argv.index( '-verbose' )
    del( argv[ pos ] )
    verbose = 1

outfiles = argv[1:]

num_models = '0.01'

RMS_THRESHOLD = 2.0

numfiles = len( outfiles )


for outfile in outfiles:
    if not exists( outfile ):
        print "Cannot find", outfile
        continue

    ##################################################################################
    # Lowscore decoys.
    ##################################################################################

    TINKER_CHARMM_SCOREFILE = 0
    CHARM_SCOREFILE = 0
    if outfile[-3:] == '.sc': # TINKER output?
        TINKER_CHARMM_SCOREFILE = 1
    if outfile.count( 'CHARMM' ):
        CHARMM_SCOREFILE = 1

    cluster_logfile = outfile.replace('.out','.cluster.log' )
    if TINKER_CHARMM_SCOREFILE:
        cluster_logfile = outfile.replace('.sc','.cluster.log' )

    if TINKER_CHARMM_SCOREFILE:
        listfile_scorecut = outfile.replace('.sc','.list' )

        if not exists ( cluster_logfile ):
            lines = open( outfile ).readlines()
            score_and_tag = []
            for line in lines[1:]:
                cols = string.split( line[:-1] )
                score_and_tag.append( ( float(cols[1]), cols[-1] ) )
            score_and_tag.sort()
            NUM_TAGS = int( len( score_and_tag ) * float( num_models ) + 0.5 )
            print NUM_TAGS
            fid = open( listfile_scorecut,'w')
            for i in range( NUM_TAGS ):
                tag = score_and_tag[i][-1]
                pdbname =  outfile.replace('_minimize.sc','_OUT') + '/' + \
                    tag.replace('minimize_','')+'_OUT/'+tag + '.pdb'
                print pdbname
                if not exists( pdbname ):
                    pdbname =  outfile.replace('_minimize.sc','_OUT') + '/' + \
                        tag.replace('minimize_','')+'.min_pdb'
                    pdbname = pdbname.replace( '_CHARMM','')
                    if not exists( pdbname ):
                        pdbname =  outfile.replace('_minimize.sc','_OUT') + '/' + \
                            tag.replace('minimize_','').replace('.pdb','')+'_OUT/'+tag + '.min_pdb'
                        pdbname = pdbname.replace( '_CHARMM','')
                assert( exists( pdbname ) )

                pdbname_RNA = pdbname.replace('.pdb','_RNA.pdb').replace('S_','s_')

                if not exists( pdbname_RNA ):
                    system( '~rhiju/python/make_rna_rosetta_ready.py '+pdbname )
                #print pdbname_RNA
                assert( exists( pdbname_RNA ) )
                fid.write( pdbname_RNA+'\n' )
            fid.close()
    else:
        outfile_scorecut = outfile.replace('.out','.low%s.out' % num_models )

    if not( exists( cluster_logfile ) ):
        if TINKER_CHARMM_SCOREFILE:
            lines = open( outfile ).readlines()
            score_and_tag = []
            for line in lines[1:]:
                cols = string.split( line[:-1] )
                score_and_tag.append( ( float(cols[1]), cols[-1] ) )
                score_and_tag.sort()
            NUM_TAGS = int( len( score_and_tag ) * float( num_models ) )
            listfile_scorecut = outfile.replace('.sc','.list' )
            fid = open( listfile_scorecut,'w')
            for i in range( NUM_TAGS ):
                tag = score_and_tag[i][-1]
                pdbname =  outfile.replace('_minimize.sc','_OUT').replace('_CHARMM','') + '/' + \
                    tag.replace('minimize_','').replace('.pdb','')+'_OUT/'+tag.replace('.pdb','') + '.pdb'
                assert( exists( pdbname ) )
                pdbname_RNA = pdbname.replace('.pdb','_RNA.pdb')
                if not exists( pdbname_RNA ):
                    system( '~rhiju/python/make_rna_rosetta_ready.py '+pdbname )
                assert( exists( pdbname_RNA ) )
                fid.write( pdbname_RNA+'\n' )
            fid.close()
        else:
            outfile_scorecut = outfile.replace('.out','.low%s.out' % num_models )

            if not exists( outfile_scorecut ):
                command = '~rhiju/python/extract_lowscore_decoys_outfile.py %s %s > %s '% (outfile, num_models, outfile_scorecut)
                print( command )
                system( command )

            if verbose: print 'Extracting low energy decoys into ', outfile_scorecut,
            numdecoys = int(popen( 'grep SCORE '+outfile_scorecut+' | wc ' ).readlines()[-1].split()[0]) - 1
            if verbose: print " ==> ", numdecoys, " decoys"

        ##################################################################################
        # Cluster
        ##################################################################################
        CLUSTER_EXE = '/users/rhiju/src/mini/bin/cluster.macosgccrelease'
        if not( exists( CLUSTER_EXE ) ):
            CLUSTER_EXE = '/work/rhiju/src/mini/bin/cluster.linuxgccrelease'
        if not( exists( CLUSTER_EXE ) ):
            CLUSTER_EXE = '/home/rhiju/src/mini/bin/cluster.linuxgccrelease'
        assert( exists( CLUSTER_EXE) )

        if TINKER_CHARMM_SCOREFILE:
            native_tag = ''

            pos = outfile.index( 'chunk' )
            rna_name = outfile[pos:(pos+13)]
            native_pdb = '~rhiju/projects/rna_new_benchmark/bench_final/%s_RNA.pdb' % rna_name
            print native_pdb
            native_tag = '-native '+native_pdb

            command = '%s -database ~/minirosetta_database  -l %s -in:file:fullatom -score:weights rna_hires.wts   -radius %f %s > %s' % ( CLUSTER_EXE, listfile_scorecut, RMS_THRESHOLD, native_tag, cluster_logfile )
            print( command )
            system( command )
        else:
            remark_line = popen( 'head -n 3  '+outfile_scorecut ).readlines()[2]
            if ( len(remark_line.split())  > 1 and remark_line.split()[1] == 'BINARY_SILENTFILE' ) :
                binary_tag = ' -in:file:silent_struct_type binary_rna '
            else:
                binary_tag = ''
            command = '%s -database ~/minirosetta_database  -in:file:silent %s -in:file:fullatom -score:weights rna_hires.wts  %s -radius %f > %s' % ( CLUSTER_EXE, outfile_scorecut, binary_tag,  RMS_THRESHOLD, cluster_logfile )
            print( command )
            system( command )

    lines = open( cluster_logfile ).readlines()
    rmsds = []
    NUM_CLUSTERS = 5
    num_members = []
    for i in range( len( lines) ):
        line = lines[i]
        cols =  string.split( line )
        if len( cols ) > 2:
            if cols[0] == 'Cluster:':
                cluster_num = int( cols[1] )
                num_members.append( int( cols[3] ) )
                if cluster_num < NUM_CLUSTERS:
                    rmsds.append( [string.split( lines[i+1] )[2], cluster_num+1] )

    for i in range( NUM_CLUSTERS ):
        for j in range(3):
            clusterfile = 'c.%s.%d.pdb' % (i,j)
            if exists( clusterfile ):
                command = 'mv %s %s.cluster%s.%s.pdb' % (clusterfile, outfile.replace('.out',''),i+1,j )
                system( command )
            #else:
            #    print clusterfile, 'missing!'

    command = 'rm -rf c.*pdb'
    system( command )

    if verbose:
        for i in range( NUM_CLUSTERS ):
            if i < len( rmsds ) :
                best_cluster_file = '%s.cluster%s.pdb' % ( outfile.replace('.out',''),i+1)
                print best_cluster_file, '==>', rmsds[i][0], "    [ N =",num_members[i],"]"
        print

    rmsds.sort()
    best_cluster_file = '%s.cluster%s.pdb' % ( outfile.replace('.out',''),rmsds[0][1])
    print best_cluster_file, '==>', rmsds[0][0], "    [ N =",num_members[0],"]"
