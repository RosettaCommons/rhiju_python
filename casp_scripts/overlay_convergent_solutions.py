#!/usr/bin/python

import string
from glob import glob
from sys import argv
from os import system,popen,chdir,getcwd
from os.path import exists,basename,dirname
from parse_options import parse_options
from math import sqrt,acos,pi

OVERRIDE = parse_options( argv, "override", 0 )
ignore_parent = parse_options( argv, "ignore_parent", 0 )

globdirs = glob( '*/')
dirs = []
for dir in globdirs:
    if dir == 'pdb/': continue
    dirs.append( dir )
dirs.sort()

outfile = 'region_FINAL.out'

def load_CA( file ):
    lines = open( file ).readlines()
    model_xyzs = []
    for line in lines:
        if len( line ) > 16 and line[12:16] == ' CA ':
            model_xyzs.append( [float(line[30:38]), float(line[38:46]), float(line[46:54])] )
    return model_xyzs

all_pdbs = []
all_loop_CAs = []
which_dir = []
count = 0

cwd = getcwd()

def get_dev2( coord1, coord2 ):
    dev2 = 0.0
    for k in range(3):
        dev = ( coord1[k] - coord2[k] )
        dev2 += dev*dev
    return dev2

def get_angle( coord1, coord2, coord3 ):
    dotproduct = 0.0
    norm1 = 0.0
    norm2 = 0.0
    for k in range(3):
        dotproduct += ( coord1[k] - coord2[k] )*(coord3[k] - coord2[k] )
        norm1 += ( coord1[k] - coord2[k] )*(coord1[k] - coord2[k] )
        norm2 += ( coord3[k] - coord2[k] )*(coord3[k] - coord2[k] )
    angle = (180.0/pi) * acos(  dotproduct/ sqrt( norm1*norm2) )
    return angle

CA_CA_CUTOFF = 4.2
CA_CA_CUTOFF2 = CA_CA_CUTOFF * CA_CA_CUTOFF
def is_there_a_chainbreak( CA ):

    for i in range( len(CA) - 1 ):
        if ( get_dev2( CA[i], CA[i+1] )  > CA_CA_CUTOFF2): return True

    #for i in range( len(CA) - 2 ):
    #    if ( get_angle( CA[i], CA[i+1], CA[i+2] ) < 45.0 ): return True

    return False

maxmodels = 50
for dir in dirs:
    chdir( dir )
    if not exists( outfile ):
        print 'Need to download: region_FINAL.out'
        if not OVERRIDE: exit( 0 )

    if exists( outfile ) and not exists( outfile+'.1.pdb' ):
        system( 'extract_lowscore_decoys.py %s %d' % (outfile, maxmodels) )

    ################################
    #What are the loop residues?
    ################################
    assert( exists( 'README_SETUP' ) )
    lines = string.join( open( 'README_SETUP' ).readlines() )
    cols = string.split( lines )
    loop_res = parse_options( cols, "loop_res", [-1] )
    assert( len( loop_res ) > 0 )

    pdbs = []

    # Starting pdb might be OK -- keep it as a candidate.
    if ( not ignore_parent ):
        pdb = 'mini_' + dir.replace('/','') + '_H.sup.pdb'
        if not( exists( pdb ) ):
            print "Hey! where is: %s/%s" % (dir,pdb)

            if not OVERRIDE: exit( 0 )

            pdb = 'mini_' + dir.replace('/','') + '_H.pdb'
            assert( exists( pdb ) )

        pdbs.append( pdb )

    for i in range(1, maxmodels + 1):
        pdb = outfile+'.%d.pdb' %  i
        if exists( pdb  ): pdbs.append( pdb )

    assert( len( pdbs ) > 0 )

    loop_CAs = []
    not_ok = []
    for m in range( len( pdbs ) ):
        pdb = pdbs[ m ]
        ########################################
        # Read out Calphas in loop
        #  -- assume they are already aligned!!
        ########################################
        CA = load_CA( pdb )
        loop_CA = []
        for i in range( len( CA ) ):
            if (i+1) in loop_res:
                loop_CA.append( CA[i] )
        loop_CAs.append( loop_CA )

        if ( is_there_a_chainbreak( loop_CA ) ):
            print dir,pdb, "CHAINBREAK", is_there_a_chainbreak( loop_CA )
            not_ok.append( m )

    if len( not_ok ) == len( pdbs ):
        print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        print "  WARNING!  WARNING! WARNING! WARNING! "
        print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        print "All pdbs have a chainbreak -- "
        print " using as placeholder: ",  pdbs[0]
        del( not_ok[0] )

    for i in range( len( pdbs ) ):
        pdb = pdbs[ i ]
        if i in not_ok: continue
        all_pdbs.append( dir+'/' + pdb )
        which_dir.append( count )
        all_loop_CAs.append( loop_CAs[ i ] )

    chdir( cwd )
    count += 1

# Calc pairwise rmsd between models.
def get_rmsd( CA1, CA2 ):
    assert( len( CA1 ) == len( CA2 ) )
    dev2 = 0.0
    for i in range( len( CA1 ) ):
        dev2 += get_dev2( CA1[i], CA2[i] )
    dev2 /= len( CA1 )
    return sqrt( dev2 )

rmsds = []
for i in range( len( all_loop_CAs ) ):
    loop_CA1 = all_loop_CAs[i]

    rmsds.append( [] )
    for j in range( len( all_loop_CAs ) ):
        loop_CA2 = all_loop_CAs[j]

        rmsds[i].append( get_rmsd( loop_CA1, loop_CA2 ) )

# output...
PRINT_ALL = 0
if PRINT_ALL:
    for i in range( len( all_pdbs ) ):
        print '%50s' % all_pdbs[i],

        for j in range( len(all_pdbs) ):
            print '%4.2f' % rmsds[i][j] ,
        print


# Need to find representative for each server
mean_min_rmsds_to_other_servers = []
for i in range( len( dirs ) ):

    for m in range( len( all_pdbs ) ):
        if (which_dir[m] != i ):continue

        rmsds_to_other_server = {}
        for n in range( len( all_pdbs ) ):

            dirnum = which_dir[n]
            if ( dirnum == i ):continue

            if dirnum not in rmsds_to_other_server.keys():
                rmsds_to_other_server[ dirnum ] = []

            rmsds_to_other_server[ dirnum ].append( rmsds[m][n] )

        min_rmsd_to_other_server = []
        for other_dir in rmsds_to_other_server.keys():
            rmsds_to_other_server[ other_dir ].sort()
            min_rmsd_to_other_server.append( rmsds_to_other_server[ other_dir ][0] )

        product = 1.0
        for x in min_rmsd_to_other_server: product *= x

        #print product, all_pdbs[m]
        mean_min_rmsds_to_other_servers.append( [product, all_pdbs[m] ] )

mean_min_rmsds_to_other_servers.sort()
central_pdb = mean_min_rmsds_to_other_servers[0][1]
central_index = all_pdbs.index( central_pdb )

representative_pdbs = []
print
for i in range( len( dirs ) ):

    rmsds_to_central = []
    for m in range( len( all_pdbs ) ):
        if (which_dir[m] != i ):continue
        rmsds_to_central.append( [ rmsds[m][central_index], all_pdbs[m] ] )

    rmsds_to_central.sort()
    representative_pdb = rmsds_to_central[0][1]
    representative_pdbs.append(  representative_pdb )
    print representative_pdb
print


if not exists( 'pdb' ): system( 'mkdir -p pdb/' )

save_pdbs = []
for i in range( len( dirs) ):
    save_pdb = '%s_%s.pdb' % ( basename( cwd ), dirs[i].replace('/','') )
    # print save_pdb
    system( 'cp -rf %s pdb/%s' % (representative_pdbs[i], save_pdb) )
    save_pdbs.append( save_pdb )

chdir( 'pdb' )
command = 'superimpose.py '
for pdb in save_pdbs: command += ' '+pdb
command += ' > q'
print command
system( command )
print
print 'rs pdb/TEST.script'

fid = open( 'TEST.script','w' )

fid.write( 'load '+getcwd()+'/q\n')

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

fid.write('set vectps on\n')
fid.write('set specular on\n')
fid.write('restrict not hydrogen\n')

for m in range(1,len( CA )+1 ):
    if m in loop_res: continue
    fid.write( 'select %d\n' % m)
    fid.write( 'color white\n')
fid.write('select all\n')

fid.close()
