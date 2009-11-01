#!/usr/bin/python

from sys import argv,exit
from os import system,popen
from os.path import exists,basename,dirname
from glob import glob
import string

verbose = 0
if argv.count( '-verbose'):
    pos = argv.index( '-verbose' )
    del( argv[ pos ] )
    verbose = 1

outfiles = argv[1:]

num_models = '0.01'

RMS_THRESHOLD = 2.5

numfiles = len( outfiles )


CLUSTER_EXE = '/users/rhiju/src/mini/bin/cluster.macosgccrelease'
if not( exists( CLUSTER_EXE ) ):
    CLUSTER_EXE = '/work/rhiju/src/mini/bin/cluster.linuxgccrelease'
if not( exists( CLUSTER_EXE ) ):
    CLUSTER_EXE = '/home/rhiju/src/mini/bin/cluster.linuxgccrelease'
assert( exists( CLUSTER_EXE) )

RNA_TEST_EXE = CLUSTER_EXE.replace( 'cluster','rna_test' )
assert( exists( RNA_TEST_EXE) )

all_cluster_info_lines = []

for outfile in outfiles:
    if not exists( outfile ):
        print "Cannot find", outfile
        continue

    pos = outfile.index( 'chunk' )
    rna_name = outfile[pos:(pos+13)]
    native_pdb = '/Users/rhiju/projects/rna_new_benchmark/bench_final/%s_RNA.pdb' % rna_name
    #print native_pdb

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

                pdbname_RNA = dirname(pdbname)+'/'+string.lower(basename(pdbname)).replace('.pdb','_RNA.pdb').replace('minimize_S_','minimize_s_')

                if not exists( pdbname_RNA ):
                    system( '~rhiju/python/make_rna_rosetta_ready.py '+pdbname )
                print pdbname_RNA
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
                pdbname_RNA = dirname(pdbname)+'/'+string.lower(basename(pdbname)).replace('.pdb','_RNA.pdb').replace( 'minimize_S','minimize_s' )
                if not exists( pdbname_RNA ):
                    system( '~rhiju/python/make_rna_rosetta_ready.py '+pdbname )
                print pdbname_RNA
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
        if TINKER_CHARMM_SCOREFILE:
            native_tag = ''

            native_tag = '-native '+native_pdb

            command = '%s -database ~/minirosetta_database  -l %s -in:file:fullatom -score:weights rna_hires.wts   -radius %f %s > %s' % ( CLUSTER_EXE, listfile_scorecut, RMS_THRESHOLD, native_tag, cluster_logfile )
            print( command )
            system( command )
        else:
            remark_line = popen( 'head -n 3  '+outfile_scorecut ).readlines()[2]
            if ( len(remark_line.split())  > 1 and remark_line.split()[1] == 'BINARY_SILENTFILE' ) :
                binary_tag = ' -in:file:silent_struct_type binary_rna '
            else:
                binary_tag = ' -in:file:silent_struct_type rna'
            command = '%s -database ~/minirosetta_database  -in:file:silent %s -in:file:fullatom -score:weights rna_hires.wts  -mute all %s -radius %f -sort_groups_by_energy -remove_singletons > %s' % ( CLUSTER_EXE, outfile_scorecut, binary_tag,  RMS_THRESHOLD, cluster_logfile )
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
            new_clusterfile = '%s.cluster%s.%s.pdb' % (outfile.replace('.out',''),i+1,j )
            clusterfile = 'c.%s.%d.pdb' % (i,j)
            if exists( clusterfile ):
                command = 'mv %s %s' % (clusterfile, new_clusterfile)
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

    rmsds_to_sort = rmsds
    rmsds_to_sort.sort()
    best_cluster_file = '%s.cluster%s.pdb' % ( outfile.replace('.out',''),rmsds_to_sort[0][1])
    #print best_cluster_file, '==>', rmsds[0][0], "    [ N =",num_members[0],"]"

    ##################################
    # New ... actually do a rescore.
    globfiles = glob( outfile.replace('.out','')+'.cluster*pdb' )
    globfiles.sort()

    cluster_scorefile = outfile+'.cluster_rms.out'
    if not exists( cluster_scorefile ):
        command = '%s -rna_stats -database ~/minirosetta_database  -s %s -in:file:fullatom -native %s -out:file:silent %s' % ( RNA_TEST_EXE, string.join( globfiles ), native_pdb, cluster_scorefile )
        print( command )
        system( command )

    assert( exists( cluster_scorefile ) )
    lines = open(cluster_scorefile).readlines()
    cols = lines[1].split()
    rms_index = cols.index( 'rms' )
    f_natNWC_index = cols.index( 'f_natNWC' )


    find_tags = []
    for line in lines[2:]:
        cols = line.split()
        tag = cols[-1]
        #if ( tag.count( '.1.pdb' ) or tag.count( '.0.pdb' ) ) :
        if ( tag.count( '.0.pdb' ) ) :
            find_tags.append( tag )

    filter_tags = []
    for tag in find_tags:
        if tag.count( '.0.pdb' ) and (tag.replace( '.0.pdb' , '.1.pdb' ) in find_tags): continue
        filter_tags.append( tag )


    cluster_info = []
    filter_cluster_info = []
    for line in lines[2:]:
        cols = line.split()
        tag = cols[-1]
        f_natNWC = float( cols[ f_natNWC_index ]  )
        f_rms    = float( cols[ rms_index ] )

        if tag in filter_tags:
            filter_cluster_info.append( ( f_rms, -1 * f_natNWC, tag ) )
        cluster_info.append( ( f_rms, -1 * f_natNWC, tag ) )

    filter_cluster_info.sort()

    #####################################
    # How many NWC base pairs were there anyway?
    native_stats_file = native_pdb+'.stats.out'
    if not exists( native_stats_file ):
        command = '%s -rna_stats -database ~/minirosetta_database  -s %s -in:file:fullatom -native %s -out:file:silent %s' % ( RNA_TEST_EXE, native_pdb, native_pdb, native_stats_file )
        print( command )
        system( command )

    lines = popen( 'cat '+native_stats_file ).readlines()
    cols = string.split( lines[1] )
    n_nwc_col = cols.index( 'N_NWC' )
    N_NWC = float( string.split( lines[2] )[ n_nwc_col ] )

    best_cluster_info = filter_cluster_info[ 0 ]

    #Need cluster size info
    best_cluster_name = best_cluster_info[2]
    pos = best_cluster_name.index('.cluster')
    cluster_num = int( best_cluster_name[ pos+8 ] ) - 1
    cluster_size = num_members[ cluster_num ]

    cluster_center_rms = -999
    cluster_center_fNWC = -0.99
    cluster_lowE_rms = -999
    cluster_lowE_fNWC = -0.99

    for info in cluster_info:
        if info[-1] == best_cluster_name.replace('.1.','.0.'):
            cluster_center_rms = info[0]
            cluster_center_fNWC = -1 * info[1]
        if info[-1] == best_cluster_name.replace('.0.','.1.'):
            cluster_lowE_rms = info[0]
            cluster_lowE_fNWC = -1 * info[1]

    all_cluster_info_lines.append( [ max( cluster_center_fNWC, cluster_lowE_fNWC ), N_NWC, cluster_center_rms, best_cluster_name] )

    print "Cluster %1d [Size %3d]. CENTER: fNWC %6.3f  rms %5.2f ; LOWEST_ENERGY: fNWC %6.3f  rms %5.2f. [max fNWC %2d ] %s " % \
        ( cluster_num+1,
          cluster_size,
          cluster_center_fNWC,
          cluster_center_rms,
          cluster_lowE_fNWC,
          cluster_lowE_rms,
          N_NWC, best_cluster_name )


N_NWC_TOT = 0
recovered_NWC_TOT = 0
N_motifs = 0

N_rmsd2 = 0
#N_rmsd2_5 = 0
N_good = 0
for info_line in all_cluster_info_lines:
    N_motifs += 1
    N_NWC_TOT += info_line[ 1 ]
    recovered_NWC_TOT += info_line[ 0 ] * info_line[ 1 ]
    if info_line[2] < 2.00:
        N_rmsd2 += 1
        N_good += 1
        print info_line[-1]
    elif info_line[ 0 ] > 0.99:
        N_good += 1
        print info_line[-1]
print

frac_NWC_TOT = 0
if N_NWC_TOT > 0: frac_NWC_TOT = float( recovered_NWC_TOT)/N_NWC_TOT

print 'Frac. NWC recovered: %6.3f  (total = %3d)' % ( frac_NWC_TOT, N_NWC_TOT )
print 'N (rmsd < 2.0 )             :  %2d out of %2d' % ( N_rmsd2, N_motifs )
print 'N (rmsd < 2.0 or f_NWC = 1) :  %2d out of %2d' % ( N_good, N_motifs )

