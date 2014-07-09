#!/usr/bin/python
## make mammoth structure alignments

import string
from glob import glob
from sys import argv,stderr,exit
from os import popen,system
from os.path import exists,basename
from operator import add
from math import sqrt
from parse_options import parse_options

#############################
def Help():
    print '\n'
    print '-'*75
    print 'USAGE: %s <pdb1> <pdb2> {... <pdbN>} > <superposition-pdb>'%argv[0]
    print '\n will superimpose pdbs 2-N onto pdb1 using maxsub, so seqs should agree'
    print '-'*75
    print '\n\n'
    exit()

if len(argv) <=2:
    Help()

MAMMOTH = '~rhiju/src/mammoth2/mammoth_rna'

RENUMBER_ATOMS = 0
SHOW_MODEL_0 = 1

model_count = 0
atom_count = 0

args = argv[1:]
R_DEFINED = 0
if args.count('-R'):
    pos = args.index('-R')
    rmsd_threshold = float(args[pos+1])
    del args[pos]
    del args[pos]
    distance_threshold = rmsd_threshold  #+ 3
    R_DEFINED = 1
else:
#    rmsd_threshold = 0.0
    rmsd_threshold = 4.0

if args.count('-D'):
    pos = args.index('-D')
    distance_threshold = float(args[pos+1])
    assert(R_DEFINED)

CALC_PER_RESIDUE_DEVIATIONS = parse_options( args, "per_res", 0 )
DUMP = parse_options( args, "dump", 0 )
COPY_RESNUM = parse_options( args, "copy_resnum", 0 )
rmsd_threshold = parse_options( args, "R", 0.0 );
RENUMBER_ATOMS = parse_options( args, "renumber_atoms", 0.0 );
COPY_HETATM = parse_options( args, "copy_hetatm", 0.0 );
slicenum = parse_options( args, "N", -1 )

if ( rmsd_threshold > 0.0 ):
    R_DEFINED = True
    distance_threshold = rmsd_threshold
    stderr.write( 'using distance threshold %5.1f A; reporting number of residues within %5.1f A \n' % (distance_threshold, rmsd_threshold))

if slicenum > 0:
    slice  = 1
else:
    rmsd_threshold = 4.0
    slice = 0

subset_residues = parse_options( args, "subset", [-1] );
subset_residues1 = parse_options( args, "subset_residue1", [-1] );
subset_residues2 = parse_options( args, "subset_residue2", [-1] );
if len( subset_residues ) > 0:
    subset_residues1 = subset_residues
    subset_residues2 = subset_residues
use_subset = ( len( subset_residues1 ) > 0 and len( subset_residues2 ) > 0 )


# this is kind of dangerous...
if args.count('-1'):
    del args[args.index('-1')]
    SHOW_MODEL_0 = 0

pdb_list = args
pdb1 = pdb_list[0]
all_per_res_dev = []

for pdb in pdb_list:
    if not exists( pdb ):
        stderr.write( 'Could not find '+pdb+'\n' )
        Help()



for pdb in pdb_list:

    pdb1_to_superimpose = pdb1
    pdb_to_superimpose = pdb

    if slice:
        command = 'pdbslice.py %s 1 %d blah_' % (pdb1, slicenum)
        system(command)
        pdb_to_superimpose = 'blah_'+pdb1

    if use_subset:

        command = 'pdbslice.py '+pdb1
        command += ' -subset '
        for i in subset_residues1:
            command += ' %d ' % i
        command += ' blah_ '
        system(command)
        pdb1_to_superimpose = 'blah_'+basename(pdb1)

        command = 'pdbslice.py '+pdb
        command += ' -subset '
        for i in subset_residues2:
            command += ' %d ' % i
        command += ' blah_ '
        system(command)
        pdb_to_superimpose = 'blah_'+basename(pdb)


    if R_DEFINED:
        command = '%s -R %f -D %f -p %s -e %s 2> /dev/null | grep PSI.end'\
                          %(MAMMOTH, rmsd_threshold,distance_threshold,pdb1_to_superimpose,pdb_to_superimpose)

    else:
        command = '%s -p %s -e %s 2> /dev/null | grep PSI.end'\
             %(MAMMOTH, pdb1_to_superimpose,pdb_to_superimpose)
#        command = '/work/pbradley/maxsub/maxsub -p %s -e %s 2> /dev/null | grep PSI.end'\
#                  %(pdb_to_superimpose,pdb)

#    if rmsd_threshold:
#        command = '/work/pbradley/mammoth2/test.out -R %f -p %s -e %s 2> /dev/null | grep PSI.end'\
#              %(rmsd_threshold,pdb1,pdb)
#    else:
#        command = '/work/pbradley/mammoth/mastodon -p %s -e %s 2> /dev/null | grep PSI.end'\
#                  %(pdb1,pdb)

    #stderr.write(command+'\n')
    print command
    lines = popen(command).readlines()

    if slice:
        command = 'rm blah_'+basename(pdb1)
        system(command)

    if use_subset:
        command = 'rm blah_'+basename(pdb1)
        system(command)
        command = 'rm blah_'+basename(pdb)
        system(command)

    if not lines:
        stderr.write('empty file? %s\n'%pdb)
        continue

    l = string.split(lines[0])

    stderr.write('%s -vs- %s: %s over %s residues\n'%(pdb_list[0],pdb,l[7],l[3]))


    file = 'maxsub_sup.pdb'
    matrix = map(lambda x:map(float,string.split(x)[1:]),
                 popen('grep -A3 "Transformation Matrix" %s'%file).readlines()[1:4])


    P_translation = map(float,
                        string.split(popen('grep -A1 "Translation vector (Pred" %s'\
                                           %file).readlines()[1])[1:])

    E_translation = map(float,
             string.split(popen('grep -A1 "Translation vector (Exp" %s'\
                                %file).readlines()[1])[1:])

    def E_transform(v,matrix,tP,tE):
        ans = [0.0]*3
        for i in range(3):
            for j in range(3):
                ans[i] = ans[i] + matrix[i][j]*(v[j]-tE[j])
            ans[i] = ans[i] + tP[i]

        return ans

    if model_count == 0:
        model0_resnum = []
        hetatm_lines = []

        model0_xyzs = []

        if SHOW_MODEL_0:
            print 'MODEL     %4d'%model_count
            model_count = model_count+1

            data = open(pdb1,'r')
            line = data.readline()

            prev_resnum = ''
            while line:
                if line[:4] in ['ATOM','HETA']:
                    atom_count = atom_count + 1
                    print '%s%5d%s'%(line[:6],atom_count,line[11:-1])

                    if (line[12:16]==' CA ' or line[12:16]==' C4*' or line[12:16]==' C4''') \
                            and CALC_PER_RESIDUE_DEVIATIONS:
                        model0_xyzs.append( [float(line[30:38]), float(line[38:46]), float(line[46:54])] )

                    if COPY_RESNUM:
                        resnum = line[22:26]
                        if resnum != prev_resnum:
                            model0_resnum.append( resnum )
                        prev_resnum = resnum

                    if line[:4] == 'HETA':
                        hetatm_lines.append(line[:66])
                elif line[:6] == 'CONECT': print line[:-1]
                elif line[:6] == 'ENDMDL':break
                line = data.readline()
            data.close()
            print 'ENDMDL'
        else:
            model_count = model_count + 1

    print 'MODEL     %4d'%model_count
    model_count = model_count+1
    data = open(pdb,'r')
    line = data.readline()

    if DUMP:
        sup_model_filename = pdb.replace( '.pdb','' ) + '.sup.pdb'
        fid_model = open( sup_model_filename, 'w' )

    prev_resnum = ''
    rescount = -1

    atom_count = 0
    per_res_dev = []
    while line:
        if line[:4] in ['ATOM','HETA']:
                atom_count = atom_count + 1
                pos = E_transform(map(float,[line[30:38],line[38:46],line[46:54]]),
                                  matrix,
                                  P_translation,
                                  E_translation)

                current_resnum = line[22:26]
                if prev_resnum != current_resnum:
                    rescount += 1
                prev_resnum = current_resnum
                if COPY_RESNUM:
                    new_resnum = model0_resnum[rescount]
                    line = line[0:22]+new_resnum+line[26:]

                if (line[12:16]==' CA ' or line[12:16]==' C4*')  and CALC_PER_RESIDUE_DEVIATIONS:
                    per_res_dev.append( sqrt( ( model0_xyzs[rescount][0] - pos[0] )*( model0_xyzs[rescount][0] - pos[0] ) + \
                                              ( model0_xyzs[rescount][1] - pos[1] )*( model0_xyzs[rescount][1] - pos[1] ) + \
                                              ( model0_xyzs[rescount][2] - pos[2] )*( model0_xyzs[rescount][2] - pos[2] ) ) )


                if RENUMBER_ATOMS:
                    new_line = '%s%5d%s%8.3f%8.3f%8.3f%s'\
                          %(line[:6],atom_count,line[11:30],pos[0],pos[1],pos[2],line[54:-1])
                else:
                    new_line = '%s%s%8.3f%8.3f%8.3f%s'\
                          %(line[:6],line[6:30],pos[0],pos[1],pos[2],line[54:-1])
                print new_line
                if DUMP: fid_model.write( new_line+'\n' )

        elif line[:6] == 'CONECT': print line[:-1]
        elif line[:6] == 'ENDMDL':break
        line = data.readline()

    if COPY_HETATM:
        for line in hetatm_lines:
            atom_count += 1
            if RENUMBER_ATOMS:
                print '%s%5d%s' % (line[:6],atom_count,line[11:])
            else:
                print line[:-1]

    data.close()

    print 'ENDMDL'
    if DUMP:
        fid_model.close()
        stderr.write(  'Created: %s\n' % sup_model_filename )

    all_per_res_dev.append( per_res_dev )

if CALC_PER_RESIDUE_DEVIATIONS:
    for n in range(  len( all_per_res_dev[ 0 ]  ) ):
        stderr.write( '%d ' % (n+1) )
        for i in range( len( all_per_res_dev ) ) :
            stderr.write( ' %8.5f' % all_per_res_dev[ i ][ n ] )
        stderr.write('\n')


system('rm -rf maxsub_sup.pdb maxsub_sup2.pdb rasmol.tcl')
