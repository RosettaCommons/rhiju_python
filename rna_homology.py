#!/usr/bin/python
"""
rna_homology
9/16/2011
Matthew Seetin
Rhiju Das Lab

rna_homology performs all-atom homology modeling of RNAs.  It takes as the first input a fasta-format file containing the alignment of
two homologous RNAs, template first, target second, and a dot-bracket notation of the secondary structure common to both RNAs. Simple pseudoknots
are supported in the secondary structure using [] instead of ().  The second input is a PDB format file of the coordinates
of the template structure.  Third, it takes an integer listing the desired number of output structures.  Finally, it takes a
file name prefix for the temporary files generated during the run.  This prefix may contain a file path if it is desired that
these files be written to somewhere other than the present working directory.

The program first goes through and quickly remodels the Watson-Crick pairs that have changed nucleotide identity one at a time
in a series of short calls to the rna_denovo of the rosetta package.  Next, it identifies which loop regions that have changed
nucleotide identity are adjascent to one another and remodels all those nucleotides together with a series of longer calls to
rna_denovo.  This program assumes that all the relevant Rosetta executables and python scripts are in your PATH and PYTHONPATH.

Usage:
rna_homology.py input.fasta template.pdb #_models prefix

"""

def Usage():
    print "Input Error."
    print "Usage: "
    print "rna_homology.py input.fasta template.pdb #_decoys prefix [conserved tertiary restraints]"

rosetta_build = "macosgccrelease"
rosetta_root = "~/src/rosetta_TRUNK/"
use_tertiary = False
tertiary = ""

import sys, os, math
import fasta_for_rna_homology, rna_for_rna_homology
from util_for_rna_homology import make_tag_with_dashes


if(len(sys.argv)  < 5 or len(sys.argv) > 6):
    Usage()
    sys.exit()

infasta = fasta_for_rna_homology.Fasta_RNA(sys.argv[1])
infasta.init_rnas()
infasta.rnas[0].readcoms(sys.argv[2])
numdecoys = int(sys.argv[3])
prefix = sys.argv[4]
if len(sys.argv) == 6:
    tertiary = sys.argv[5]
    use_tertiary = True

infasta.align_by_number(0, 1)
infasta.finddifferences(0, 1)
current_pdb = sys.argv[2]

tertiary_bases = []
tertiary_restraints = []
insert_pair_bases = []
insert_pair_restraints = []
do_not_remodel = []
force_remodel = []

if(use_tertiary):
    print "Reading tertiary restraints from:", tertiary
    for line in open(tertiary,"r"):
	column = line.split()
	if line[:14] == "CONSERVED PAIR":
	    tertiary_bases.append(int(column[2]))
	    tertiary_bases.append(int(column[3]))
	    tertiary_restraints.append(column[2:])
	elif line[:11] == "INSERT PAIR":
	    insert_pair_bases.append(int(column[2]))
	    insert_pair_bases.append(int(column[3]))
	    insert_pair_restraints.append(column[2:])
	elif line[:7] == "REMODEL":
	    for j in range(1,len(column)):
		force_remodel.append(int(column[j]))
	elif line[:3] == "FIX":
	    for j in range(1,len(column)):
		do_not_remodel.append(int(column[j]))

if(use_tertiary):
    for i in range(len(tertiary_restraints)):
	tertiary_restraints[i][0] = int(tertiary_restraints[i][0])
	tertiary_restraints[i][1] = int(tertiary_restraints[i][1])
    for i in range(len(insert_pair_restraints)):
	insert_pair_restraints[i][0] = int(insert_pair_restraints[i][0])
	insert_pair_restraints[i][1] = int(insert_pair_restraints[i][1])


current_seq = []
for i in range(1,len(infasta.rnas[0].seq)):
    current_seq.append(infasta.rnas[0].seq[i])

wc_pair_mutation_pairs = []

EXE = "rna_denovo."+rosetta_build+" -database "+rosetta_root+"rosetta_database/"
EXE += " -score:weights rna/rna_hires_07232011_with_intra_base_phosphate.wts  -rna_torsion_potential RNA11_based_new  -geom_sol_correct_acceptor_base "

for base in infasta.wc_pair_mutations:
    wc_pair_mutation_pairs.append(infasta.rnas[0].basepr[base])

##########################################################################################################
# 1. Remodel WC base pairs first.  This is done one base pair at a time to preserve homology with input pdb.
##########################################################################################################
if len( infasta.wc_pair_mutations )  > 0:
    print "####################################################################################"
    print "Remodeling these pairs:\n", infasta.wc_pair_mutations, "\n", wc_pair_mutation_pairs
    print "####################################################################################"

os.system("rm -f "+prefix+".*.cut.out")

for i in range(0,len(infasta.wc_pair_mutations)):

    print "###############################"
    print "Remodel pair:",infasta.wc_pair_mutations[i],"--",infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]]
    print "###############################"
    mutstring = ".mut"+str(infasta.wc_pair_mutations[i])+"-"+str(infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]])
    outpdb = open(prefix+mutstring+".pdb","w")
    currnuc = 0
    prevres = ""
    for line in open(current_pdb,"r"):
	column = line.split()
	if len(column) > 6:
	    if column[5] != prevres:
		    prevres = column[5]
		    currnuc += 1
	    if currnuc != infasta.wc_pair_mutations[i] and currnuc != infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]]:
		outpdb.write(line)
	else:
	    outpdb.write(line)
    outpdb.close()

#    os.system("convert_pdb_to_silent_file.linuxgccrelease  -output_silent_file "+ \
    os.system("convert_pdb_to_silent_file."+rosetta_build+"  -output_silent_file "+ \
    prefix + mutstring + ".cut.out -database " + rosetta_root + "rosetta_database/ -s " + prefix + mutstring + \
    ".pdb")

    #write main prm file
    outprm = open(prefix+mutstring+".prm","w")
    outprm.write("STEM PAIR "+str(infasta.wc_pair_mutations[i])+" "+str(infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]])+" W W A\n")

    outprm.write("ALLOW_INSERT_RES "+str(infasta.wc_pair_mutations[i])+" "+str(infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]])+"\n")
    if infasta.wc_pair_mutations[i] > 1:
        # The X X X was dangerous as it could lead to jump-atom connections that were not at the bases.
        # This occasionally led to atom trees where fragment insertions went outside the loop! Horrible!
	outprm.write("OBLIGATE PAIR "+str(infasta.wc_pair_mutations[i]-1)+" "+str(infasta.wc_pair_mutations[i]+1)+" W W A\n")
    if infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]] < infasta.rnas[0].length:
	outprm.write("OBLIGATE PAIR "+str(infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]]-1)+" "+str(infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]]+1)+" W W A\n")

    if infasta.wc_pair_mutations[i] > 1:
	outprm.write("CUTPOINT_CLOSED "+str(infasta.wc_pair_mutations[i])+"\n")
    if infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]] < infasta.rnas[0].length:
	outprm.write("CUTPOINT_CLOSED "+str(infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]])+"\n")

    if infasta.wc_pair_mutations[i] > 1:
	outprm.write("VIRTUAL_ANCHOR 1\n")
    else:
	outprm.write("VIRTUAL_ANCHOR 2\n")

    #outprm.write("CUTPOINT_OPEN "+str(infasta.rnas[0].length)+"\n")

    #outprm.write("CUTPOINT_CLOSED "+str(infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]])+"\n")
    outprm.close()

    #write new fasta file with mutated sequence
    outfasta = open(prefix+mutstring+".mut.fasta","w")
    outseq = ""
    for j in range(1,len(infasta.rnas[0].seq)):
	if j != infasta.wc_pair_mutations[i] and j != infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]]:
	    outseq += current_seq[j-1]
	else:
	    nuc = infasta.rnas[1].seq[infasta.alignnum1.index(j)]
	    outseq += nuc
	    current_seq[j-1] = nuc
    outfasta.write(infasta.labels[0])
    outfasta.write(outseq+"\n")
    outfasta.close()

    chunk_res = []
    for j in range(1, len(infasta.rnas[0].seq)):
	if j != infasta.wc_pair_mutations[i] and j != infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]]:
            chunk_res.append( j )
    chunkstr = make_tag_with_dashes( chunk_res )

    #Now actually mutate the base pair:

    ndecoys = 5
    command = EXE+" -fasta "+prefix+mutstring+".mut.fasta  -params_file " \
        + prefix+mutstring+".prm -nstruct "+str(ndecoys)+" -out::file::silent "+ \
        prefix+mutstring+".out -cycles 20000  -in:file:silent " + prefix + mutstring+ \
        ".cut.out -close_loops -dump -output_virtual -output_lores_silent_file  -chunk_res"+chunkstr
    command  += " -close_loops_after_each_move"
    command  += " -minimize_rna"
    print command
    os.system( command )


    #Extract the structure with the lowest RMS to the native structure (or the previous mutant)
    os.system("extract_lowscore_decoys.py "+prefix+mutstring+".out 1 -no_virtual")

    cleanpdb = open(prefix+mutstring+".1.pdb","w")
    for line in open(prefix+mutstring+".out.1.pdb","r"):
	column = line.split()
	if column[0] == "TER":
	    cleanpdb.write(line)
	elif len(column) > 4:
	    if "OV" not in column[2]:
		cleanpdb.write(line)
    cleanpdb.close()

    #now perform this procedure on the next base pair using the PDB from the previous mutant
    current_pdb = prefix+mutstring+".1.pdb"


##########################################################################################################
# 2. Second, we remodel the user-defined conserved tertiary interactions, if they need remodeling.
#     The user specifying these when the sequence is unchanged will prevent them from being remodeled in a loop calculation.
##########################################################################################################
if(use_tertiary):
    for i in range(0,len(tertiary_restraints)):
	col1 = 0
	col2 = 0
	tempseq = ""
	for k in range(0,len(infasta.alignments[0])):
	    if infasta.alignments[0][k] in ["a","c", "g", "u"]:
		tempseq += infasta.alignments[0][k]
	    if len(tempseq) == tertiary_restraints[i][0]:
		col1 = k
	    if len(tempseq) == tertiary_restraints[i][1]:
		col2 = k

	if infasta.alignments[0][col1] not in ["a","c", "g", "u"] or infasta.alignments[0][col2] not in ["a","c", "g", "u"] \
	or infasta.alignments[1][col1] not in ["a","c", "g", "u"] or infasta.alignments[1][col2] not in ["a","c", "g", "u"]:
	    print "Error: The user-specified tertiary interaction between bases", tertiary_restraints[i][0], "and", tertiary_restraints[i][1],\
	    "does not appear to be conserved but rather is an insertion.  Please check the input alignment and the numbering of the interaction."
	    sys.exit()

	if infasta.alignments[0][col1] == infasta.alignments[1][col1] and  infasta.alignments[0][col2] == infasta.alignments[1][col2]:
	    continue

	mutstring = ".mut"+str(infasta.alignnum0.index(tertiary_restraints[i][0]))+"-"+str(infasta.alignnum0.index(tertiary_restraints[i][1]))
	outpdb = open(prefix+mutstring+".pdb","w")
	currnuc = 0
	prevres = ""
	for line in open(current_pdb,"r"):
	    column = line.split()
	    if len(column) > 6:
		if column[5] != prevres:
			prevres = column[5]
			currnuc += 1
		if currnuc != infasta.alignnum0.index(tertiary_restraints[i][0]) and currnuc != infasta.alignnum0.index(tertiary_restraints[i][1]):
		    outpdb.write(line)
	    else:
		outpdb.write(line)
	outpdb.close()

	#os.system("convert_pdb_to_silent_file.linuxgccrelease  -output_silent_file "+ \
	os.system("convert_pdb_to_silent_file."+rosetta_build+"  -output_silent_file "+ \
	prefix + mutstring + ".cut.out -database " + rosetta_root + "rosetta_database/ -s " + prefix + mutstring + \
	".pdb")

	#write main prm file
	outprm = open(prefix+mutstring+".prm","w")
	outprm.write("OBLIGATE PAIR "+str(infasta.alignnum0.index(tertiary_restraints[i][0]))+" "+str(infasta.alignnum0.index(tertiary_restraints[i][1]))+" "+tertiary_restraints[i][2]+" "+tertiary_restraints[i][3]+" "+tertiary_restraints[i][4]+"\n")

	outprm.write("ALLOW_INSERT_RES "+str(infasta.alignnum0.index(tertiary_restraints[i][0]))+" "+str(infasta.alignnum0.index(tertiary_restraints[i][1]))+"\n")

	if tertiary_restraints[i][0] > 1 and tertiary_restraints[i][0] < infasta.rnas[0].length:
	    #outprm.write("OBLIGATE PAIR "+str(infasta.alignnum0.index(tertiary_restraints[i][0])-1)+" "+str(infasta.alignnum0.index(tertiary_restraints[i][1])+1)+" X X X\n")
	    outprm.write("OBLIGATE PAIR "+str(infasta.alignnum0.index(tertiary_restraints[i][0])-1)+" "+str(infasta.alignnum0.index(tertiary_restraints[i][1])+1)+" W W A\n")
	if tertiary_restraints[i][1] > 1 and tertiary_restraints[i][1] < infasta.rnas[0].length:
	    #outprm.write("OBLIGATE PAIR "+str(infasta.alignnum0.index(tertiary_restraints[i][1])-1)+" "+str(infasta.alignnum0.index(tertiary_restraints[i][1])+1)+" X X X\n")
	    outprm.write("OBLIGATE PAIR "+str(infasta.alignnum0.index(tertiary_restraints[i][1])-1)+" "+str(infasta.alignnum0.index(tertiary_restraints[i][1])+1)+" W W A\n")

	if tertiary_restraints[i][0] > 1 and tertiary_restraints[i][0] < infasta.rnas[0].length:
	    outprm.write("CUTPOINT_CLOSED "+str(infasta.alignnum0.index(tertiary_restraints[i][0]))+"\n")
	if tertiary_restraints[i][1] > 1 and tertiary_restraints[i][1] < infasta.rnas[0].length:
	    outprm.write("CUTPOINT_CLOSED "+str(infasta.alignnum0.index(tertiary_restraints[i][1]))+"\n")
	"""if infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]] < infasta.rnas[0].length:
	    outprm.write("VIRTUAL_ANCHOR "+str(infasta.rnas[0].length)+"\n")
	else:
	    outprm.write("VIRTUAL_ANCHOR 1\n")"""
	#outprm.write("CUTPOINT_OPEN "+str(infasta.rnas[0].length)+"\n")

	#outprm.write("CUTPOINT_CLOSED "+str(infasta.rnas[0].basepr[infasta.wc_pair_mutations[i]])+"\n")
	outprm.close()

	#write new fasta file with mutated sequence
	outfasta = open(prefix+mutstring+".mut.fasta","w")
	outseq = ""
	for j in range(1,len(infasta.rnas[0].seq)):
	    if j != infasta.alignnum0.index(tertiary_restraints[i][1]) and j != tertiary_restraints[i][1]:
		outseq += current_seq[j-1]
	    else:
		nuc = infasta.rnas[1].seq[infasta.alignnum1.index(j)]
		outseq += nuc
		current_seq[j-1] = nuc
	outfasta.write(infasta.labels[0])
	outfasta.write(outseq+"\n")
	outfasta.close()

        chunk_res = []
	for j in range(1, len(infasta.rnas[0].seq)):
	    if j != tertiary_restraints[i][0] and j != tertiary_restraints[i][1]:
                chunk_res.append( j )
        chunkstr = make_tag_with_dashes( chunk_res )

	#Now actually mutate the base pair:
        command = EXE+" -fasta "+prefix+mutstring+".mut.fasta  -params_file " \
            + prefix+mutstring+".prm -nstruct 10 -out::file::silent "+ \
            prefix+mutstring+".out -cycles 5000 -minimize_rna "\
            " -in:file:silent " + prefix + mutstring+ \
            ".cut.out -close_loops -chunk_res"+chunkstr
        command += " -close_loops_after_each_move"
        os.system( command )

	#Extract the structure with the lowest RMS to the native structure (or the previous mutant)
	os.system("extract_lowscore_decoys.py "+prefix+mutstring+".out 1")

	cleanpdb = open(prefix+mutstring+".1.pdb","w")
	for line in open(prefix+mutstring+".out.1.pdb","r"):
	    column = line.split()
	    if column[0] == "TER":
		cleanpdb.write(line)
	    elif len(column) > 4:
		if "OV" not in column[2]:
		    cleanpdb.write(line)
	cleanpdb.close()

	#now perform this procedure on the next base pair using the PDB from the previous mutant
	current_pdb = prefix+mutstring+".1.pdb"
    for i in range(0,len(do_not_remodel)):
	if do_not_remodel[i] in infasta.insertions:
	    print "Error.  Inserted bases must be remodeled.  Please restate your restraints so that base",do_not_remodel[i],"is not fixed."
	    sys.exit()
	if infasta.rnas[1].basepr[do_not_remodel[i]] > 0:
	    continue #base pairs will have already been handled.
	elif infasta.rnas[1].seq[do_not_remodel[i]] == infasta.rnas[0].seq[infasta.alignnum0.index(infasta.rnas[1].seq[do_not_remodel[i]])]:
	    continue #bases are idential, no remodeling necessary
	else:
	    print "Warning.  Base",do_not_remodel[i],"in the target sequence was set as fixed, but the nucleotide identity has changed. \
	    Remodeling this base on its own and then leaving it fixed."

	    mutstring = ".mut"+str(infasta.alignnum0.index(do_not_remodel[i]))
	    outpdb = open(prefix+mutstring+".pdb","w")
	    currnuc = 0
	    prevres = ""
	    for line in open(current_pdb,"r"):
		column = line.split()
		if len(column) > 6:
		    if column[5] != prevres:
			    prevres = column[5]
			    currnuc += 1
		    if currnuc != infasta.alignnum0.index(do_not_remodel[i]):
			outpdb.write(line)
		else:
		    outpdb.write(line)
	    outpdb.close()

            #os.system("convert_pdb_to_silent_file.linuxgccrelease  -output_silent_file "+ \
	    os.system("convert_pdb_to_silent_file."+rosetta_build+"  -output_silent_file "+ \
	    prefix + mutstring + ".cut.out -database " + rosetta_root + "rosetta_database/ -s " + prefix + mutstring + \
	    ".pdb")

	    #write main prm file
	    outprm = open(prefix+mutstring+".prm","w")
	    outprm.write("ALLOW_INSERT_RES "+str(infasta.alignnum0.index(do_not_remodel[i]))+"\n")
	    if infasta.alignnum0.index(do_not_remodel[i]) > 1 and infasta.alignnum0.index(do_not_remodel[i]) < infasta.rnas[0].length:
		outprm.write("OBLIGATE PAIR "+str(infasta.alignnum0.index(do_not_remodel[i])-1)+" "+str(infasta.alignnum0.index(do_not_remodel[i])+1)+" X X X\n")
	    if infasta.alignnum0.index(do_not_remodel[i]) > 1 and infasta.alignnum0.index(do_not_remodel[i]) < infasta.rnas[0].length:
		outprm.write("CUTPOINT_CLOSED "+str(infasta.alignnum0.index(do_not_remodel[i]))+"\n")
	    outprm.close()

	    #write new fasta file with mutated sequence
	    outfasta = open(prefix+mutstring+".mut.fasta","w")
	    outseq = ""
	    for j in range(1,len(infasta.rnas[0].seq)):
		if j != infasta.alignnum0.index(do_not_remodel[i]):
		    outseq += current_seq[j-1]
		else:
		    nuc = infasta.rnas[1].seq[do_not_remodel[i]]
		    outseq += nuc
		    current_seq[j-1] = nuc
	    outfasta.write(infasta.labels[0])
	    outfasta.write(outseq+"\n")
	    outfasta.close()

	    chunk_res = []
	    for j in range(1, len(infasta.rnas[0].seq)):
		if j != infasta.alignnum0.index(do_not_remodel[i]):
		    chunkstr += " " + str(j)
                    chunk_res.append( j )
            chunkstr = make_tag_with_dashes( chunk_res )

	    #Now actually mutate the base pair:
            # remove native for now
            #	    os.system("rna_denovo."+rosetta_build+" -fasta "+prefix+mutstring+".mut.fasta -native "+current_pdb+" -params_file " \
	    os.system(EXE+" -fasta "+prefix+mutstring+".mut.fasta  -params_file " \
	    + prefix+mutstring+".prm -nstruct 20 -out::file::silent "+ \
	    prefix+mutstring+".out -cycles 10000 -minimize_rna  -in:file:silent " + prefix + mutstring+ \
	    ".cut.out -close_loops -chunk_res"+chunkstr)

	    #Extract the structure with the lowest RMS to the native structure (or the previous mutant)
	    os.system("extract_lowscore_decoys.py "+prefix+mutstring+".out 1")

	    cleanpdb = open(prefix+mutstring+".1.pdb","w")
	    for line in open(prefix+mutstring+".out.1.pdb","r"):
		column = line.split()
		if column[0] == "TER":
		    cleanpdb.write(line)
		elif len(column) > 4:
		    if "OV" not in column[2]:
			cleanpdb.write(line)
	    cleanpdb.close()

	    #now perform this procedure on the next base pair using the PDB from the previous mutant
	    current_pdb = prefix+mutstring+".1.pdb"

###########################################################################################################################
# 3. Deletions at or continuous from either the 5' or the 3' ends require no remodeling.  They will just be deleted.
###########################################################################################################################
five_prime_dele=0
three_prime_dele=infasta.rnas[0].length+1
for i in range(1,infasta.rnas[0].length+1):
    if five_prime_dele+1 in infasta.deletions_in_template:
	five_prime_dele += 1
	idx = infasta.deletions_in_template.index(five_prime_dele)
	infasta.deletions.pop(idx)
	infasta.deletions_in_template.pop(idx)
    else:
	break

for i in range(infasta.rnas[0].length,1,-1):
    if three_prime_dele-1 in infasta.deletions_in_template:
	three_prime_dele -= 1
	idx = infasta.deletions_in_template.index(three_prime_dele)
	infasta.deletions.pop(idx)
	infasta.deletions_in_template.pop(idx)
    else:
	break

if five_prime_dele > 0 or three_prime_dele < infasta.rnas[0].length+1:
    outpdb = open(prefix+".delete5p3p.pdb","w")
    currnuc = 0
    prevres = ""
    for line in open(current_pdb,"r"):
	column = line.split()
	if len(column) > 6:
	    if column[5] != prevres:
		    prevres = column[5]
		    currnuc += 1
	    if currnuc > five_prime_dele and currnuc < three_prime_dele:
		outpdb.write(line)
	else:
	    outpdb.write(line)
    outpdb.close()
    current_pdb = prefix+".delete5p3p.pdb"


###########################################################################################################################
# 4. Next is remodeling loops. This takes place in two stages.  First is identifying which nucleotides need to be remodeled.
#    Second is identifying which nucleotides need to be remodeled simultaneously.
#    If two "loop regions" have at least two nucleotides with their centroids within 5.0 A of each other,
#    they will be remodeled simultaneously.
###########################################################################################################################
scratch_loop_regions = []
scratch_loop_regions.append([])
curr_region = 0
obligate_pairs = []
obligate_pairs.append([])
scratch_delete = []
scratch_delete.append([])
template_tertiary_bases = []
template_tertiary_bases.append(-1)

if(use_tertiary):
    for base in tertiary_bases:
	if base not in template_tertiary_bases:
	    template_tertiary_bases.append(infasta.alignnum0.index(base))
    template_tertiary_bases.pop(0)

template_tertiary_bases.sort()

print "######################################################"
print "Starting analysis of which bases need to be remodeled."
print "######################################################"

###############################################################################################################################
#Insertions.
if len( infasta.insertions ) > 0: print "Scanning through potential insert positions: ", infasta.insertions

# If following is True, insertions trigger the two neighboring bases (as long as those bases have not yet been
# set to be remodeled by another mutation).
#remodel_neighbors = True
remodel_neighbors = False

for insert in infasta.insertions:
    redundant_base = False

    for region in scratch_loop_regions:
	if insert in region:
            # Already slated for remodeling
	    redundant_base = True
	    break

    if redundant_base: continue

    scratch_loop_regions[curr_region].append(insert)
    #Don't append this base to scratch_delete as it doesn't correspond to any template base

    end_insert = insert+1
    while end_insert in infasta.insertions:
        scratch_loop_regions[curr_region].append(end_insert)
        end_insert += 1

    moving_base = False
    still_5p_base = insert - 1

    if remodel_neighbors and (insert-1 not in do_not_remodel) and (still_5p_base not in tertiary_bases) and still_5p_base > 0:
        scratch_loop_regions[curr_region].append(insert-1)
        scratch_delete[curr_region].append(infasta.alignnum0.index(insert-1))
        still_5p_base = insert-2

    if (still_5p_base in infasta.insertions or still_5p_base in infasta.loop_mutations or still_5p_base in infasta.deletions \
            or still_5p_base in infasta.insertions_plus or still_5p_base in infasta.insertions_minus or \
            still_5p_base in infasta.deletions_plus or still_5p_base in force_remodel) \
            and still_5p_base not in do_not_remodel and still_5p_base not in tertiary_bases:
        moving_base = True
        scratch_loop_regions[curr_region].append(still_5p_base)
        if still_5p_base in infasta.alignnum0:
            scratch_delete[curr_region].append(infasta.alignnum0.index(still_5p_base))
        still_5p_base -= 1

    while moving_base and still_5p_base >= 0:
        if (still_5p_base in infasta.insertions or still_5p_base in infasta.loop_mutations or still_5p_base in infasta.deletions \
                or still_5p_base in infasta.insertions_plus or still_5p_base in infasta.insertions_minus or \
                still_5p_base in infasta.deletions_plus or still_5p_base in force_remodel) \
                and still_5p_base not in do_not_remodel and still_5p_base not in tertiary_bases:
            scratch_loop_regions[curr_region].append(still_5p_base)
            if still_5p_base in infasta.alignnum0:
                scratch_delete[curr_region].append(infasta.alignnum0.index(still_5p_base))
            still_5p_base -= 1
        else:
            in_region = False
            for region in scratch_loop_regions:
                if still_5p_base in region:
                    still_5p_base -= 1
                    in_region = True
                    break
            moving_base = in_region

    still_3p_base = end_insert

    if remodel_neighbors and still_3p_base not in do_not_remodel and still_3p_base not in tertiary_bases and still_3p_base <= infasta.rnas[1].length:
        scratch_loop_regions[curr_region].append(end_insert)
        scratch_delete[curr_region].append(infasta.alignnum0.index(end_insert))
        still_3p_base = end_insert+1


    if (still_3p_base in infasta.insertions or still_3p_base in infasta.loop_mutations or still_3p_base in infasta.deletions \
            or still_3p_base in infasta.insertions_plus or still_3p_base in infasta.insertions_minus or \
            still_3p_base in infasta.deletions_plus or still_3p_base in force_remodel) \
            and (still_3p_base not in do_not_remodel and still_3p_base not in tertiary_bases):
        moving_base = True
        scratch_loop_regions[curr_region].append(still_3p_base)
        if still_3p_base in infasta.alignnum0:
            scratch_delete[curr_region].append(infasta.alignnum0.index(still_3p_base))
        still_3p_base += 1

    while moving_base and still_3p_base <= infasta.rnas[1].length:
        if (still_3p_base in infasta.insertions or still_3p_base in infasta.loop_mutations or still_3p_base in infasta.deletions\
                or still_3p_base in infasta.insertions_plus or still_3p_base in infasta.insertions_minus or \
                still_3p_base in infasta.deletions_plus or still_3p_base in force_remodel) \
                and (still_3p_base not in do_not_remodel and still_3p_base not in tertiary_bases ):
            scratch_loop_regions[curr_region].append(still_3p_base)
            if still_3p_base in infasta.alignnum0:
                scratch_delete[curr_region].append(infasta.alignnum0.index(still_3p_base))
            still_3p_base += 1

        else:
            in_region = False
            for region in scratch_loop_regions:
                if still_3p_base in region:
                    still_3p_base += 1
                    in_region = True
                    break
            moving_base = in_region

    redundant_pair = False
    for pair_set in obligate_pairs:
        for k in range(0,len(pair_set),2):
            if pair_set[k] == still_5p_base:
                redundant_pair = True
                break
        for k in range(1,len(pair_set),2):
            if pair_set[k] == still_3p_base:
                redundant_pair = True
                break
        if redundant_pair:
            break

    if not redundant_pair:
        if still_5p_base > 0 and still_3p_base <= infasta.rnas[1].length:
            obligate_pairs[curr_region].append(still_5p_base)
            obligate_pairs[curr_region].append(still_3p_base)
        elif still_5p_base < 1:
            assert( still_3p_base <= infasta.rnas[1].length)
            obligate_pairs[curr_region].append( -999 ) # fake number -- signal below that this is not really a 'pair', but the 3' residue is a terminus.
            obligate_pairs[curr_region].append(still_3p_base)
        else:
            assert( still_5p_base > 0)
            obligate_pairs[curr_region].append(still_5p_base)
            obligate_pairs[curr_region].append( -999 ) # fake number -- signal below that this is not really a 'pair', but the 5' residue is a terminus.


    curr_region += 1
    scratch_loop_regions.append([])
    obligate_pairs.append([])
    scratch_delete.append([])

###############################################################################################################################
#Deletions.  Insertions trigger the nucleotides on either side to be remodeled.  Deletions at the ends are already gone.
prevdele = 0

deleted_bases = []
if len( infasta.deletions ) > 0: print "Scanning through potential delete positions: ", infasta.deletions
for i in range(len(infasta.deletions)):

    if infasta.deletions[i] not in deleted_bases:

	deleted_bases.append(infasta.deletions[i])

	scratch_delete[curr_region].append(infasta.deletions_in_template[i])

	j = i+1
	if j < len(infasta.deletions):
	    while scratch_delete[curr_region][-1] == infasta.deletions_in_template[j] - 1:
		scratch_delete[curr_region].append(infasta.deletions_in_template[j])
		j += 1
		if j >= len(infasta.deletions):
		    break

	scratch_loop_regions[curr_region].append(infasta.deletions[i])
	scratch_loop_regions[curr_region].append(infasta.deletions[i]+1)

	if infasta.deletions[i] in infasta.alignnum0:
	    scratch_delete[curr_region].append(infasta.alignnum0.index(infasta.deletions[i]))

	if infasta.deletions[i]+1 in infasta.alignnum0:
	    scratch_delete[curr_region].append(infasta.alignnum0.index(infasta.deletions[i]+1))

	moving_base = False
	still_5p_base = infasta.deletions[i]-1
	still_3p_base = infasta.deletions[i]+2


	if (still_5p_base in infasta.insertions or still_5p_base in infasta.loop_mutations or still_5p_base in infasta.deletions \
	or still_5p_base in infasta.insertions_plus or still_5p_base in infasta.insertions_minus or still_5p_base in infasta.deletions_plus or still_5p_base in force_remodel) \
	and still_5p_base not in do_not_remodel and still_5p_base not in tertiary_bases and still_5p_base > 0:
	    moving_base = True
	    scratch_loop_regions[curr_region].append(still_5p_base)
	    if still_5p_base in infasta.alignnum0:
		scratch_delete[curr_region].append(infasta.alignnum0.index(still_5p_base))
	    still_5p_base -= 1

	while moving_base and still_5p_base >= 0:
	    if (still_5p_base in infasta.insertions or still_5p_base in infasta.loop_mutations or still_5p_base in infasta.deletions \
	    or still_5p_base in infasta.insertions_plus or still_5p_base in infasta.insertions_minus or still_5p_base in infasta.deletions_plus or still_5p_base in force_remodel) \
	    and still_5p_base not in do_not_remodel and still_5p_base not in tertiary_bases:
		scratch_loop_regions[curr_region].append(still_5p_base)
		if still_5p_base in infasta.alignnum0:
		    scratch_delete[curr_region].append(infasta.alignnum0.index(still_5p_base))
		still_5p_base -= 1

	    else:
		in_region = False
		for region in scratch_loop_regions:
		    if still_5p_base in region:
			still_5p_base -= 1
			in_region = True
			break
		moving_base = in_region

	if (still_3p_base in infasta.insertions or still_3p_base in infasta.loop_mutations or still_3p_base in infasta.deletions \
	or still_3p_base in infasta.insertions_plus or still_3p_base in infasta.insertions_minus or still_3p_base in infasta.deletions_plus or still_3p_base in force_remodel) \
	and still_3p_base not in do_not_remodel and still_3p_base not in tertiary_bases and still_3p_base < infasta.rnas[1].length:
	    moving_base = True
	    scratch_loop_regions[curr_region].append(still_3p_base)
	    if still_3p_base in infasta.alignnum0:
		scratch_delete[curr_region].append(infasta.alignnum0.index(still_3p_base))
	    still_3p_base += 1
	while moving_base and still_3p_base <= infasta.rnas[1].length:
	    if (still_3p_base in infasta.insertions or still_3p_base in infasta.loop_mutations or still_3p_base in infasta.deletions \
	    or still_3p_base in infasta.insertions_plus or still_3p_base in infasta.insertions_minus or still_3p_base in infasta.deletions_plus or still_3p_base in force_remodel) \
	    and still_3p_base not in do_not_remodel and still_3p_base not in tertiary_bases:
		scratch_loop_regions[curr_region].append(still_3p_base)
		if still_3p_base in infasta.alignnum0:
		    scratch_delete[curr_region].append(infasta.alignnum0.index(still_3p_base))
		still_3p_base += 1
	    else:
		in_region = False
		for region in scratch_loop_regions:
		    if still_3p_base in region:
			still_3p_base += 1
			in_region = True
			break
		moving_base = in_region


	redundant_pair = False
	for pair_set in obligate_pairs:
	    for k in range(0,len(pair_set),2):
		if pair_set[k] == still_5p_base:
		    redundant_pair = True
		    break
	    for k in range(1,len(pair_set),2):
		if pair_set[k] == still_3p_base:
		    redundant_pair = True
		    break
	    if redundant_pair:
		break

	if still_5p_base > 0 and still_3p_base <= infasta.rnas[1].length and not redundant_pair:
	    obligate_pairs[curr_region].append(still_5p_base)
	    obligate_pairs[curr_region].append(still_3p_base)

	curr_region += 1
	scratch_loop_regions.append([])
	obligate_pairs.append([])
	scratch_delete.append([])


###############################################################################################################################
#Mutations do not necessarily remodel the whole loop.  Only the mutated bases are remodeled
i = 0
if len( infasta.loop_mutations ) > 0:
    print "##################################################################"
    print "Scanning through potential loop mutation positions: ", infasta.loop_mutations
    print "##################################################################"


# force_remodel used to be a different loop than infasta.loop_mutations -- now unifying with 'regular' loop mutations
all_remodel = infasta.loop_mutations + force_remodel

for i in range(len(all_remodel)):

    redundant_base = False
    for region in scratch_loop_regions:
	if infasta.loop_mutations[i] in region:
            # already slated for remodeling
	    redundant_base = True
	    break

    if not redundant_base and all_remodel[i] not in tertiary_bases:
	scratch_loop_regions[curr_region].append(all_remodel[i])
	scratch_delete[curr_region].append(infasta.alignnum0.index(all_remodel[i]))

        # determine loop boundaries
	end_mut = all_remodel[i]+1
	while end_mut in all_remodel and end_mut not in tertiary_bases:
	    scratch_loop_regions[curr_region].append(end_mut)
	    scratch_delete[curr_region].append(infasta.alignnum0.index(end_mut))
	    end_mut += 1

	moving_base = False
	still_5p_base = all_remodel[i]-1
	still_3p_base = end_mut

        # is loop boundary itself moveable? If so move it back
        # need to ask Matt about this -- why necessary?
	if (still_5p_base in infasta.insertions or still_5p_base in infasta.loop_mutations or still_5p_base in infasta.deletions \
                or still_5p_base in infasta.insertions_plus or still_5p_base in infasta.insertions_minus or still_5p_base in infasta.deletions_plus or still_5p_base in force_remodel) \
                and still_5p_base not in do_not_remodel and still_5p_base not in tertiary_bases:
	    moving_base = True
	    still_5p_base -= 1
	while moving_base and still_5p_base >= 0:
	    if (still_5p_base in infasta.insertions or still_5p_base in infasta.loop_mutations or still_5p_base in infasta.deletions \
	    or still_5p_base in infasta.insertions_plus or still_5p_base in infasta.insertions_minus or still_5p_base in infasta.deletions_plus or still_5p_base in force_remodel) \
	    and still_5p_base not in do_not_remodel and still_5p_base not in tertiary_bases:
		still_5p_base -= 1
	    else:
		in_region = False
		for region in scratch_loop_regions:
		    if still_5p_base in region:
			still_5p_base -= 1
			in_region = True
			break
		moving_base = in_region

	if (still_3p_base in infasta.insertions or still_3p_base in infasta.loop_mutations or still_3p_base in infasta.deletions \
                or still_3p_base in infasta.insertions_plus or still_3p_base in infasta.insertions_minus or still_3p_base in infasta.deletions_plus or still_3p_base in force_remodel) \
                and still_3p_base not in do_not_remodel and still_3p_base not in tertiary_bases:
	    moving_base = True
	    still_3p_base += 1
	while moving_base and still_3p_base <= infasta.rnas[1].length:
	    if (still_3p_base in infasta.insertions or still_3p_base in infasta.loop_mutations or still_3p_base in infasta.deletions \
	    or still_3p_base in infasta.insertions_plus or still_3p_base in infasta.insertions_minus or still_3p_base in infasta.deletions_plus or still_3p_base in force_remodel) \
	    and still_3p_base not in do_not_remodel and still_3p_base not in tertiary_bases:
		still_3p_base += 1
	    else:
		in_region = False
		for region in scratch_loop_regions:
		    if still_3p_base in region:
			still_3p_base += 1
			in_region = True
			break
		moving_base = in_region

	redundant_pair = False
	for pair_set in obligate_pairs:
	    for k in range(0,len(pair_set),2):
		if pair_set[k] == still_5p_base:
		    redundant_pair = True
		    break
	    for k in range(1,len(pair_set),2):
		if pair_set[k] == still_3p_base:
		    redundant_pair = True
		    break
	    if redundant_pair:
		break

	if still_5p_base > 0 and still_3p_base <= infasta.rnas[1].length and not redundant_pair:
	    obligate_pairs[curr_region].append(still_5p_base)
	    obligate_pairs[curr_region].append(still_3p_base)

	curr_region += 1
	scratch_loop_regions.append([])
	obligate_pairs.append([])
	scratch_delete.append([])

print "scratch_loop_regions:", scratch_loop_regions
print "obligate_pairs:      ", obligate_pairs
exit()

###############################################################################################################################
#Check if loops in the template have nucleotides within 13 A of each other
#If so, combine into a single region.  Many loops may be combined this way.
#


curr_region = 0
local_regions = []
local_regions.append([])
combined_regions = []
combined_regions.append([])
local_pairs = []
local_pairs.append([])
delete_from_pdb = []
delete_from_pdb.append([])

# it might be better to define 'motifs' based on desired secondary structure + crystallographic tertiary contacts
# then hash out which loops get merged. -- rhiju

COM_CUTOFF = 11.0 # centroid/centroid distance cutoff
for i in range(0,len(scratch_loop_regions)):
    for j in range(i+1, len(scratch_loop_regions)):
	combine = False
	for k in range(0, len(scratch_loop_regions[i])):
	    for m in range(0, len(scratch_loop_regions[j])):

		if abs(scratch_loop_regions[i][k] - scratch_loop_regions[j][m]) <= 2:
		    combine = True
		    print "Region",i,"combined with",j,"because of sequence closeness!"

		    break
		elif infasta.rnas[1].basepr[scratch_loop_regions[i][k]] in scratch_loop_regions[j]:
		    print "Region",i,"combined with",j,"because of base pairing!"

		    combine = True
		    break
		elif infasta.rnas[1].basepr[scratch_loop_regions[j][m]] in scratch_loop_regions[i]:
		    combine = True
		    print "Region",i,"combined with",j,"because of base pairing!"

		    break
		elif infasta.rnas[1].basepr[scratch_loop_regions[i][k]+1]>0  and \
                        infasta.rnas[1].basepr[scratch_loop_regions[i][k]+1]+1 in scratch_loop_regions[j]:
		    print "Region",i,"combined with",j,"because of nearby base pairing!"
		    combine = True
		    break
                elif  infasta.rnas[1].basepr[scratch_loop_regions[i][k]-1]>0 and \
                        infasta.rnas[1].basepr[scratch_loop_regions[i][k]-1]-1 in scratch_loop_regions[j] :
		    print "Region",i,"combined with",j,"because of nearby base pairing!"
		    combine = True
		    break
		elif infasta.rnas[1].basepr[scratch_loop_regions[j][m]-1] > 0 and  \
                        infasta.rnas[1].basepr[scratch_loop_regions[j][m]-1]-1 in scratch_loop_regions[i] :
		    combine = True
		    print "Region",i,"combined with",j,"because of nearby base pairing!"
		    break
                elif infasta.rnas[1].basepr[scratch_loop_regions[j][m]+1] > 0 and \
                        infasta.rnas[1].basepr[scratch_loop_regions[j][m]+1]+1 in scratch_loop_regions[i] :
		    combine = True
		    print "Region",i,"combined with",j,"because of nearby base pairing!"
		    break
		elif use_tertiary:
		    for contact in tertiary_restraints:
			if scratch_loop_regions[i][k] in contact and scratch_loop_regions[j][m] in contact:
			    combine = True
			    print "Region",i,"combined with",j,"because of tertiary restraints!"

			    break
		    if combine:
			break
		    for contact in insert_pair_restraints:
			if scratch_loop_regions[i][k] in contact and scratch_loop_regions[j][m] in contact:
			    combine = True
			    print "Region",i,"combined with",j,"because of tertiary restraints!"

			    break
		    if combine:
			break
		elif scratch_loop_regions[i][k] in infasta.alignnum0 and scratch_loop_regions[j][m] in infasta.alignnum0:
		    if infasta.rnas[0].ntcom[infasta.alignnum0.index(scratch_loop_regions[i][k])].distance(infasta.rnas[0].ntcom[infasta.alignnum0.index(scratch_loop_regions[j][m])]) < COM_CUTOFF:
			combine = True
			print "Region",i,"combined with",j,"because of distance in the template crystal structure!"

			break
	    if combine:
		break

	if combine:
	    if len(combined_regions[0]) == 0:
		combined_regions[0].append(i)
		combined_regions[0].append(j)
		if len(obligate_pairs[i]) > 0:
		    local_pairs[0].append(obligate_pairs[i][0])
		    local_pairs[0].append(obligate_pairs[i][1])
		if len(obligate_pairs[j]) > 0:
		    local_pairs[0].append(obligate_pairs[j][0])
		    local_pairs[0].append(obligate_pairs[j][1])

		for base in scratch_loop_regions[i]:
		    local_regions[0].append(base)
		for base in scratch_loop_regions[j]:
		    local_regions[0].append(base)
		for base in scratch_delete[i]:
		    delete_from_pdb[0].append(base)
		for base in scratch_delete[j]:
		    delete_from_pdb[0].append(base)
	    else:

		i_in = -1
		j_in = -1
		for region in range(0,len(combined_regions)):
		    if i in combined_regions[region]:
			i_in = region
		    if j in combined_regions[region]:
			j_in = region

		if j_in >= 0 and i_in >= 0 and j_in == i_in:
		    continue
		elif j_in < 0 and i_in >= 0:

		    combined_regions[i_in].append(j)

		    if len(obligate_pairs[j]) > 0:
			for base in obligate_pairs[j]:
			    local_pairs[i_in].append(base)
		    for base in scratch_loop_regions[j]:
			local_regions[i_in].append(base)
		    for base in scratch_delete[j]:
			delete_from_pdb[i_in].append(base)
		elif i_in < 0 and j_in >=0 :

		    combined_regions[j_in].append(i)
		    if len(obligate_pairs[i]) > 0:
			for base in obligate_pairs[i]:
			    local_pairs[j_in].append(base)
		    for base in scratch_loop_regions[i]:
			local_regions[j_in].append(base)
		    for base in scratch_delete[i]:
			delete_from_pdb[j_in].append(base)
		elif j_in >= 0 and i_in >= 0 and j_in != i_in:
		    if j_in < i_in:
			temp = j_in
			j_in = i_in
			i_in = temp
		    for region in combined_regions[j_in]:
			combined_regions[i_in].append(region)
		    combined_regions.pop(j_in)
		    if len(local_pairs[j_in]) > 0:
			for base in local_pairs[j_in]:
			    local_pairs[i_in].append(base)
		    local_pairs.pop(j_in)
		    for base in local_regions[j_in]:
			local_regions[i_in].append(base)
		    local_regions.pop(j_in)
		    for base in delete_from_pdb[j_in]:
			delete_from_pdb[i_in].append(base)
		    delete_from_pdb.pop(j_in)
		    curr_region -= 1
		else:
		    curr_region += 1

		    local_regions.append([])
		    combined_regions.append([])
		    local_pairs.append([])
		    delete_from_pdb.append([])

		    combined_regions[curr_region].append(i)
		    combined_regions[curr_region].append(j)

		    if len(obligate_pairs[i]) > 0:
			local_pairs[curr_region].append(obligate_pairs[i][0])
			local_pairs[curr_region].append(obligate_pairs[i][1])
		    if len(obligate_pairs[j]) > 0:
			local_pairs[curr_region].append(obligate_pairs[j][0])
			local_pairs[curr_region].append(obligate_pairs[j][1])
		    for n in range(0,len(scratch_loop_regions[i])):
			local_regions[curr_region].append(scratch_loop_regions[i][n])
		    for n in range(0,len(scratch_loop_regions[j])):
			local_regions[curr_region].append(scratch_loop_regions[j][n])

		    for base in scratch_delete[i]:
			delete_from_pdb[curr_region].append(base)
		    for base in scratch_delete[j]:
			delete_from_pdb[curr_region].append(base)


# remove the empty set at the end -- if its still there.
n = len( local_regions )
if len(local_regions[n-1]) == 0:
    del( local_regions[n-1] )
    del( local_pairs[n-1] )
    del( delete_from_pdb[n-1] )

used_regions = []

for combo in combined_regions:
    for i in combo:
	used_regions.append(i)

for i in range(len(scratch_loop_regions)):
    if i not in used_regions:
	if len(scratch_loop_regions[i]) > 0:
	    local_regions.append(scratch_loop_regions[i])
	    local_pairs.append(obligate_pairs[i])
	    delete_from_pdb.append(scratch_delete[i])
print "local_regions:", local_regions
print "local_pairs:  ", local_pairs
#exit()

###############################################################################################################################
# Each loop remodeling calculation on each local region will have a different sequence and potentially a different length.
# Compute the alignment between the sequence that will be remodeled in each calculation against the target sequence
#  and adjust the numbers of the bases to be remodeled and restrained to be accordingly.
final_regions = []
final_pairs = []
restrained_pairs = []
new_alignment = []
new_alignnum = []
for i in range(len(local_regions)):

    final_pairs.append([])
    final_regions.append([])
    restrained_pairs.append([])
    new_alignment.append([])
    new_alignnum.append([0])

    tempseq = ""
    targseq = ""

    wc_pair_complements = []
    for j in range(len(infasta.wc_pair_mutations)):
	wc_pair_complements.append(infasta.rnas[0].basepr[infasta.wc_pair_mutations[j]])

    for j in range(len(infasta.alignments[0])):

	if infasta.alignments[0][j] in ["a","c", "g", "u"]:
	    tempseq += str(infasta.alignments[0][j])
	if infasta.alignments[1][j] in ["a","c", "g", "u"]:
	    targseq += str(infasta.alignments[1][j])

	if (len(tempseq) <= five_prime_dele and infasta.alignments[1][j] not in ["a","c", "g", "u"]) \
	or (len(tempseq) >= three_prime_dele and infasta.alignments[1][j] not in ["a","c", "g", "u"]):
	    continue

        # new_alignment is the sequence of the pose during the loop modeling calculation. It is what will go into the .fasta file for rna_denovo.
	if (len(tempseq) in infasta.wc_pair_mutations or len(tempseq) in wc_pair_complements or len(tempseq) in template_tertiary_bases) and infasta.alignments[0][j] in ["a","c", "g", "u"]:
	    new_alignment[i].append(infasta.alignments[1][j]) # target sequence
	elif len(targseq) in local_regions[i] or len(targseq) in do_not_remodel:
	    new_alignment[i].append(infasta.alignments[1][j]) # target sequence
	elif (len(tempseq) > five_prime_dele or five_prime_dele == 0) and (len(tempseq) < three_prime_dele or three_prime_dele == 0 ):
	    new_alignment[i].append(infasta.alignments[0][j]) # template sequence
	elif (len(tempseq) <= five_prime_dele and infasta.alignments[1][j] in ["a","c", "g", "u"]) \
	or (len(tempseq) >= three_prime_dele and infasta.alignments[1][j] in ["a","c", "g", "u"]):
	    new_alignment[i].append(infasta.alignments[0][j]) # template sequence? what?

    targseq = ""
    newseq = ""

    for j in range(len(new_alignment[i])):
	if new_alignment[i][j] in ["a","c", "g", "u"]:
	    newseq += str(new_alignment[i][j])
	if infasta.alignments[1][j+five_prime_dele] in ["a","c", "g", "u"]:
	    new_alignnum[i].append(len(newseq))

    for j in range(len(local_regions[i])):
	final_regions[i].append(new_alignnum[i][local_regions[i][j]])
    #If there's an insertion in the target sequence compared to the template in a region that's not being remodeled, there will be repeat numbers here.
    j = 0
    while j  < len(final_regions[i]):
	if final_regions[i].count(final_regions[i][j]) > 1:
	    final_regions[i].pop(j)
	    continue
	else:
	    j += 1

    #I need to remove redundancy in local_regions, too, in anticipation of the combination step.
    j = 0
    while j  < len(local_regions[i]):
	if local_regions[i].count(local_regions[i][j]) > 1:
	    local_regions[i].pop(j)
	    continue
	else:
	    j += 1
    local_regions[i].sort()

    j = 0
    while j < len(local_pairs[i]):

	#if local_pairs[i][j] > 0 and local_pairs[i][j] < len(new_alignnum[i]) and local_pairs[i][j+1] < len(new_alignnum[i]):
        if local_pairs[i][j] < len(new_alignnum[i]) and local_pairs[i][j+1] < len(new_alignnum[i]):

            # 5' terminus
            if local_pairs[i][j] > 0:
                final_pairs[i].append(new_alignnum[i][local_pairs[i][j]])
            else:
                final_pairs[i].append( -999 )

            # 3' terminus
            if local_pairs[i][j+1] > 0:
                final_pairs[i].append(new_alignnum[i][local_pairs[i][j+1]])
            else:
                final_pairs[i].append( -999 )

	    if infasta.rnas[1].basepr[j] > 0:
		restrained_pairs[i].append(new_alignnum[i][j])
		restrained_pairs[i].append(new_alignnum[i][infasta.rnas[1].basepr[j]])
	    if infasta.rnas[1].basepr[j+1] > 0:
		restrained_pairs[i].append(new_alignnum[i][j+1])
		restrained_pairs[i].append(new_alignnum[i][infasta.rnas[1].basepr[j+1]])
	else:
	    local_pairs[i].pop(j+1)
	    local_pairs[i].pop(j)
	j+=2

cutpoints = []
for i in range(len(final_regions)):
    cutpoints.append([])

#print "Final regions: ", final_regions
#print "Final pairs: ", final_pairs
#exit()

###########################################################################################################################################################
# Will eventually need to add one extra nucleotide 5' and 3' of each obligate pair from the unchanged (or only WC pair mutation) part of the structure
# It's basically an impossibility for such an extra nucleotide to be a deletion in the template, but the calculation of the numbering of this nucleotide
# is done completley just in case.
final_pairs_plus = []
local_pairs_plus = []
for i in range(len(final_pairs)):
    final_pairs_plus.append([])
    local_pairs_plus.append([])

    for j in range(len(final_pairs[i])):
	if j % 2 == 0:
	    offset = 1
	    while final_pairs[i][j] - offset > 0:
		if (final_pairs[i][j]-offset) in new_alignnum[i]:
		    final_pairs_plus[i].append(final_pairs[i][j] - offset)
		    local_pairs_plus[i].append(new_alignnum[i].index(final_pairs[i][j] - offset))
		    break
		else:
		    offset -= 1

	else:
	    offset = 1
	    while final_pairs[i][j] + offset < new_alignnum[i][-1]:
		if (final_pairs[i][j]+offset) in new_alignnum[i]:
		    final_pairs_plus[i].append(final_pairs[i][j] + offset)
		    local_pairs_plus[i].append(new_alignnum[i].index(final_pairs[i][j] + offset))
		    break
		else:
		    offset += 1



##########################################################################################################################################################
# At long last, we've figured out which bases need to be remodeled, which need to be remodeled together,
#  which bases at the ends of loops need to be restrained, and which need to be deleted from the pdb.
# Time to write some input files and run the calculations.
##########################################################################################################################################################

for i in range(len(final_regions)):
    loopstr = ".loop"+str(i)

    """print "obligate_pairs",obligate_pairs
    print "final_pairs",final_pairs
    print "final_regions",final_regions"""
    """print "scratch_loop_regions",scratch_loop_regions

    print "local_regions",local_regions
    print "combined_regions",combined_regions
    print "scratch_delete",scratch_delete
    print "delete_from_pdb",delete_from_pdb
    print "infasta.loop_mutations",infasta.loop_mutations
    print "infasta.deletions",infasta.deletions
    print "infasta.insertions",infasta.insertions
    print "deleted_bases",deleted_bases
    print "local_pairs",local_pairs"""

    # rhiju to matt: why this offset? why not just use final_pairs_plus, defined above?
    for j in range(len(final_pairs[i])):
	while final_pairs[i][j] in final_regions[i]:
	    if j % 2 == 0:
		final_pairs[i][j] -= 1
	    else:
		final_pairs[i][j] += 1

    print
    print 'Remodeling: final_region --> ',  final_regions[i]
    print 'Remodeling: final_pair   --> ',  final_pairs[i]

    #First, delete remodeled residues from the PDB
    #current_pdb is the same as from the end of the base pair remodeling
    outpdb = open(prefix+loopstr+".cut.pdb","w")
    currnuc = five_prime_dele
    prevres = ""
    for line in open(current_pdb,"r"):
	column = line.split()
	if len(column) > 6:
	    if column[5] != prevres:
		    prevres = column[5]
		    currnuc += 1
	    if currnuc not in delete_from_pdb[i]:
		outpdb.write(line)
	else:
	    outpdb.write(line)
    outpdb.close()

    #Convert this PDB to a silent file
    #os.system("convert_pdb_to_silent_file.linuxgccrelease  -output_silent_file "+ \
    os.system("convert_pdb_to_silent_file."+rosetta_build+"  -output_silent_file "+ \
    prefix + loopstr + ".cut.out -database " + rosetta_root + "rosetta_database/ -s " + prefix + loopstr + \
    ".cut.pdb")



    #write main prm file
    outprm = open(prefix+loopstr+".prm","w")
    for j in range(0,len(final_pairs[i]),2):
        if final_pairs[i][j] < 1 or final_pairs[i][j+1] < 1: continue # this is a 'fake' final_pair, denoting a terminal fragment with only one boundary residue.
	if infasta.rnas[1].basepr[local_pairs[i][j]] == local_pairs[i][j+1]:
	    outprm.write("OBLIGATE PAIR " + str(final_pairs[i][j]) + " " + str(final_pairs[i][j+1]) + " W W A\n")
	else:
            # The X X X was dangerous as it could lead to jump-atom connections that were not at the bases.
            # This occasionally led to atom trees where fragment insertions went outside the loop! Horrible!
	    #outprm.write("OBLIGATE PAIR " + str(final_pairs[i][j]) + " " + str(final_pairs[i][j+1]) + " X X X\n")
	    outprm.write("OBLIGATE PAIR " + str(final_pairs[i][j]) + " " + str(final_pairs[i][j+1]) + " W W A\n")
	outprm.write("CUTPOINT_CLOSED " + str(final_pairs[i][j+1]-1) + "\n")
	cutpoints[i].append(final_pairs[i][j+1])

    paired_bases = []
    fixed_bases = []


    for j in range(len(final_regions[i])):
	if infasta.rnas[1].basepr[new_alignnum[i].index(final_regions[i][j])] > 0 and new_alignnum[i][infasta.rnas[1].basepr[new_alignnum[i].index(final_regions[i][j])]] in final_regions[i]:
	    paired_bases.append(final_regions[i][j])


	"""if infasta.rnas[1].basepr[new_alignnum[i].index(final_regions[i][j])] > 0 and new_alignnum[i].index(final_regions[i][j]) not in paired_bases \
	and infasta.rnas[1].basepr[new_alignnum[i].index(final_regions[i][j])] not in paired_bases and infasta.rnas[1].basepr[new_alignnum[i].index(final_regions[i][j])] > new_alignnum[i].index(final_regions[i][j]):
	    paired_bases.append(infasta.rnas[1].basepr[new_alignnum[i].index(final_regions[i][j])])
	    paired_bases.append(final_regions[i][j])
	    outprm.write("STEM PAIR " + str(final_regions[i][j]) + " " +str(new_alignnum[i][infasta.rnas[1].basepr[new_alignnum[i].index(final_regions[i][j])]]) + " W W A\n")"""

    if use_tertiary:
	for contact in insert_pair_restraints:

	    if contact[0] in local_regions[i] or contact[1] in local_regions[i]:
		outprm.write("OBLIGATE PAIR "+str(new_alignnum[i][contact[0]])+" "+str(new_alignnum[i][contact[1]])+" "+contact[2]+" "+contact[3]+" "+contact[4]+"\n")
		if contact[1] in local_regions[i]:
		    outprm.write("CUTPOINT_CLOSED " + str(new_alignnum[i][contact[1]]-1) + "\n")
		else:
		    outprm.write("CUTPOINT_CLOSED " + str(new_alignnum[i][contact[0]]-1) + "\n")
		    if contact[0] not in local_regions[i]:
			cutpoints[i].append(new_alignnum[i][contact[0]])

    paired_bases.sort()

    helix_chunk_str = ""
    print paired_bases
    os.system("rm -f "+prefix+loopstr+".helix*.out")

    helix_count = 0
    while len(paired_bases) > 0:

	fiveprimeseq = ""
	threeprimeseq = ""
	helix_chunk_5p_str = ""
	helix_chunk_3p_str = ""

	if new_alignnum[i][infasta.rnas[1].basepr[new_alignnum[i].index(paired_bases[0])]] not in paired_bases:
	    paired_bases.pop(0)
	    continue
	else:
	    helix_count += 1
	    k = 1
	    while k < len(paired_bases):
		if paired_bases[k] != paired_bases[k-1]+1:
		    break
		elif new_alignnum[i][infasta.rnas[1].basepr[new_alignnum[i].index(paired_bases[k])]] not in paired_bases:
		    break
		elif new_alignnum[i][infasta.rnas[1].basepr[new_alignnum[i].index(paired_bases[k])]] != new_alignnum[i][infasta.rnas[1].basepr[new_alignnum[i].index(paired_bases[k-1])]] - 1:
		    break
		else:
		    k+=1

            helix_chunk_5p_res = []
            helix_chunk_3p_res = []
	    for n in range(k):

		if paired_bases[n]-n-1 in final_regions[i] and infasta.rnas[1].basepr[new_alignnum[i].index(paired_bases[n]-n-1)] == 0:
		    fixed_bases.append(paired_bases[n])

		fiveprimeseq += infasta.rnas[1].seq[new_alignnum[i].index(paired_bases[n])]
		threeprimeseq = infasta.rnas[1].seq[infasta.rnas[1].basepr[new_alignnum[i].index(paired_bases[n])]] + threeprimeseq
		helix_chunk_5p_res.append( paired_bases[n] )
		helix_chunk_3p_res.append(new_alignnum[i][infasta.rnas[1].basepr[new_alignnum[i].index(paired_bases[n])]])
            helix_chunk_5p_str =  make_tag_with_dashes( helix_chunk_5p_res )
            helix_chunk_3p_str =  make_tag_with_dashes( reversed(helix_chunk_3p_res) )

	    if paired_bases[k-1] + 1 not in paired_bases and new_alignnum[i][infasta.rnas[1].basepr[new_alignnum[i].index(paired_bases[k-1])]] - 1 not in paired_bases:

		outprm.write("OBLIGATE PAIR "+str(paired_bases[k-1])+" "+str(new_alignnum[i][infasta.rnas[1].basepr[new_alignnum[i].index(paired_bases[k-1])]])+" W W A\n")
                # for now can we get away with not specifying cutpoints?
		#outprm.write("CUTPOINT_CLOSED "+str(new_alignnum[i][infasta.rnas[1].basepr[new_alignnum[i].index(paired_bases[k-1])]]-1) + "\n")


	    for n in range(k):
		paired_bases.pop(0)

	    helix_fasta = open(prefix+loopstr+".helix"+str(helix_count)+".fasta","w")
	    helix_fasta.write(infasta.labels[1])
	    helix_fasta.write(fiveprimeseq+threeprimeseq+"\n")
	    helix_fasta.close()

	    os.system("rna_helix."+rosetta_build+" -fasta " + prefix+loopstr+".helix"+str(helix_count)+".fasta -database "+rosetta_root+"rosetta_database/" \
	    +" -out:file:silent "+prefix+loopstr+".helix"+str(helix_count)+".out")

	    helix_chunk_str += helix_chunk_5p_str + " " + helix_chunk_3p_str + " "
	#done making helices

    insertstr = ""
    for base in final_regions[i]:
	if base not in fixed_bases:
            insertstr += str(base)+ " "

    outprm.write("ALLOW_INSERT_RES " + insertstr+"\n")
    outprm.close()

    #restrained_pairs = []


    """if(len(restrained_pairs[i]) > 1):
	for j in range(0, len(restrained_pairs[i]), 2):
	    outprm.write("STEM PAIR "+str(restrained_pairs[i][j])+ " " + str(restrained_pairs[i][j+1])+" W W A\n")"""


    #sys.exit()
        #write new fasta file with mutated sequence
    #Also, we will need to recompute the base numbers of the nucleotides that will be remodeled depending on the final sequence
    outfasta = open(prefix+loopstr+".fasta","w")
    outseq = ""

    for j in range(len(new_alignment[i])):
	if new_alignment[i][j] in ["a","c", "g", "u"]:
	    outseq += str(new_alignment[i][j])

    outfasta.write(infasta.labels[1])
    outfasta.write(outseq+"\n")
    outfasta.close()

    chunk_res = []
    for j in range(1, len(outseq)+1):
	if j not in final_regions[i]:
            chunk_res.append( j )
    chunkstr = make_tag_with_dashes( chunk_res )

    chunkstr += " " + helix_chunk_str

    helix_file_list = " "
    for j in range(1,helix_count+1):
	helix_file_list += prefix+loopstr+".helix"+str(j)+".out "

    avg_remodel = 0
    for region in final_regions:
	avg_remodel += len(region)
    avg_remodel /= len(final_regions)

    decoys_per_loop = max(2, int(math.ceil(float(numdecoys)**(1.0/float(len(final_regions))))))
    decoys_per_loop_total = decoys_per_loop * 5 # was 10x in matt's original script.
    numcycles = 10000*len(final_regions[i])


    #################################
    #Now actually mutate the loop:
    #################################
    command =  EXE+" -fasta "+prefix+loopstr+".fasta -params_file " \
        + prefix+loopstr+".prm -nstruct "+str(decoys_per_loop_total) +" -out::file::silent "+ \
        prefix+loopstr+".out -cycles "+str(numcycles)+" -minimize_rna -in:file:silent " + prefix + loopstr+ \
        ".cut.out "+helix_file_list+" -close_loops -chunk_res"+chunkstr
    command += " -close_loops_after_each_move "

    os.system( command )

    #Extract the structure with the lowest score
    os.system("extract_lowscore_decoys.py "+prefix+loopstr+".out "+str(decoys_per_loop))

    # Parenthetical note from Rhiju -- I think this script cannot currently handle modeling 5' and 3' terminal fragments simultaneously -- might be
    #  worth refactoring.

    for j in range(1, decoys_per_loop+1):
	currnuc = 0
	cleanpdb = open(prefix+loopstr+"."+str(j)+".pdb","w")
	for line in open(prefix+loopstr+".out."+str(j)+".pdb","r"):
	    column = line.split()
	    if len(column) > 6:
		if column[5] != prevres:
			prevres = column[5]
			currnuc += 1
		if (currnuc in final_regions[i] or currnuc in final_pairs[i] or currnuc in final_pairs_plus[i]) and "OV" not in column[2]:
		    cleanpdb.write(line)
	    else:
		cleanpdb.write(line)
	cleanpdb.close()

	os.system("rm -f "+prefix+loopstr+"."+str(j)+".out") # the -f prevents a warning
        #os.system("convert_pdb_to_silent_file.linuxgccrelease  -output_silent_file "+ \
	os.system("convert_pdb_to_silent_file."+rosetta_build+"  -output_silent_file "+ \
	prefix+loopstr+"."+str(j)+".out -database " + rosetta_root + "rosetta_database/ -s " + prefix+loopstr+"."+str(j)+".pdb")


#switching from final to local to get numbering right
changed_nucs = []
for i in range(len(local_regions)):
    for base in local_regions[i]:
	changed_nucs.append(base)

final_cutpoints = []
for i in range(len(final_regions)):
    final_cutpoints.append([])
    for j in range(len(cutpoints[i])):
	final_cutpoints[i].append(new_alignnum[i].index(cutpoints[i][j]))

all_cutpoints = []
for i in range(len(final_cutpoints)):
    for base in final_cutpoints[i]:
	all_cutpoints.append(base)


currnuc = 0
prevres = ""
nonloop_nucs = []
nonlooppdb = open(prefix+".nonloop.pdb","w")
for line in open(prefix+".loop0.out.1.pdb","r"):
    column = line.split()
    if len(column) > 6:
	if column[5] != prevres:
	    prevres = column[5]
	    currnuc += 1

	    if currnuc in new_alignnum[0]:

		if new_alignnum[0].index(currnuc) not in changed_nucs and new_alignnum[0].index(currnuc) not in all_cutpoints:
		    nonloop_nucs.append(new_alignnum[0].index(currnuc))

		    if new_alignnum[0].index(currnuc) < 10:
			nonlooppdb.write(line[:22]+"   "+str(new_alignnum[0].index(currnuc))+line[26:])
		    elif new_alignnum[0].index(currnuc) < 100:
			nonlooppdb.write(line[:22]+"  "+str(new_alignnum[0].index(currnuc))+line[26:])
		    elif new_alignnum[0].index(currnuc) < 1000:
			nonlooppdb.write(line[:22]+" "+str(new_alignnum[0].index(currnuc))+line[26:])
		    else:
			nonlooppdb.write(line[:22]+str(new_alignnum[0].index(currnuc))+line[26:])

	elif currnuc in new_alignnum[0]:
	    if new_alignnum[0].index(currnuc) not in changed_nucs and new_alignnum[0].index(currnuc) not in all_cutpoints:
		if new_alignnum[0].index(currnuc) < 10:
		    nonlooppdb.write(line[:22]+"   "+str(new_alignnum[0].index(currnuc))+line[26:])
		elif new_alignnum[0].index(currnuc) < 100:
		    nonlooppdb.write(line[:22]+"  "+str(new_alignnum[0].index(currnuc))+line[26:])
		elif new_alignnum[0].index(currnuc) < 1000:
		    nonlooppdb.write(line[:22]+" "+str(new_alignnum[0].index(currnuc))+line[26:])
		else:
		    nonlooppdb.write(line[:22]+str(new_alignnum[0].index(currnuc))+line[26:])
nonlooppdb.close()

os.system("rm "+prefix+".nonloop.out")
#os.system("convert_pdb_to_silent_file."+rosetta_build+" -output_silent_file "+ \
os.system("convert_pdb_to_silent_file."+rosetta_build+"  -output_silent_file "+ \
prefix+".nonloop.out -database " + rosetta_root + "rosetta_database/ -s " + prefix+".nonloop.pdb")

chunk_res = []
for i in range(len(nonloop_nucs)):
    chunk_res.append( nonloop_nucs[i] )
chunkstr = make_tag_with_dashes( chunk_res )

for i in range(len(local_regions)):
    changed_loop_nucs = local_regions[i] + local_pairs[i] + local_pairs_plus[i]
    changed_loop_nucs.sort()
    j = 0
    while j  < len(changed_loop_nucs):
	if changed_loop_nucs.count(changed_loop_nucs[j]) > 1 or  changed_loop_nucs[j] < 1:
	    changed_loop_nucs.pop(j)
	    continue
	else:
	    j += 1

    loop_chunk_res = []
    for j in range(len(changed_loop_nucs)):
	loop_chunk_res.append( changed_loop_nucs[j] )
    loop_chunkstr = " " + make_tag_with_dashes( loop_chunk_res )
    chunkstr += loop_chunkstr


current_decoy_by_loop = []
prev_decoy_by_loop = []

for i in range(len(local_regions)):
    current_decoy_by_loop.append(1) #decoys are indexed from one.
    prev_decoy_by_loop.append(0)


outfasta = open(prefix+".final.fasta","w")
outseq = ""

for j in range(1,len(infasta.rnas[1].seq)):
    outseq += infasta.rnas[1].seq[j]

outfasta.write(infasta.labels[1])
outfasta.write(outseq+"\n")
outfasta.close()


for i in range(numdecoys):
    print "Writing "+prefix+".decoy"+str(i)+".out..."
    print "Combining"
    for j in range(len(current_decoy_by_loop)):
	print "loop", j, "decoy", current_decoy_by_loop[j]
    print ""
    outpdb = open(prefix+".decoy"+str(i)+".pdb","w")

    file_str = prefix+".nonloop.out "
    for j in range(len(local_regions)):
	file_str += prefix+".loop" + str(j)+"."+str(current_decoy_by_loop[j])+".out "


    outprm = open(prefix+".final.prm","w")
    for j in range(len(local_regions)):
	for k in range(0,len(local_pairs[j]),2):
            if local_pairs[j][k] < 1 or local_pairs[j][k+1] < 1: continue # this is a 'fake' final_pair, denoting a terminal fragment with only one boundary residue.
	    if infasta.rnas[1].basepr[local_pairs[j][k]] == local_pairs[j][k+1]:
		outprm.write("STEM PAIR " + str(local_pairs[j][k]) + " " + str(local_pairs[j][k+1]) + " W W A\n")
	    else:
		outprm.write("OBLIGATE PAIR " + str(local_pairs[j][k]) + " " + str(local_pairs[j][k+1]) + " X X X\n")
	    outprm.write("CUTPOINT_CLOSED " + str(local_pairs[j][k+1]-1) + "\n")

	paired_bases = []
	fixed_bases = []


	for k in range(len(local_regions[j])):
	    if infasta.rnas[1].basepr[local_regions[j][k]] > 0 and infasta.rnas[1].basepr[local_regions[j][k]] in local_regions[j]:
		paired_bases.append(local_regions[j][k])


	paired_bases.sort()

	helix_chunk_str = ""

	helix_count = 0
	while len(paired_bases) > 0:

	    if infasta.rnas[1].basepr[paired_bases[0]] not in paired_bases:
		paired_bases.pop(0)
		continue
	    else:
		helix_count += 1
		k = 1
		while k < len(paired_bases):
		    if paired_bases[k] != paired_bases[k-1]+1:
			break
		    elif infasta.rnas[1].basepr[paired_bases[k]] not in paired_bases:
			break
		    elif infasta.rnas[1].basepr[paired_bases[k]] != infasta.rnas[1].basepr[paired_bases[k-1]] - 1:
			break
		    else:
			k+=1
		for n in range(k):
		    print paired_bases[n], n
		    if paired_bases[n]-n-1 in local_regions[j] and infasta.rnas[1].basepr[paired_bases[n]-n-1] == 0:
			fixed_bases.append(paired_bases[n])


		if paired_bases[k-1] + 1 not in paired_bases and infasta.rnas[1].basepr[paired_bases[k-1]] - 1 not in paired_bases:

		    outprm.write("OBLIGATE PAIR "+str(paired_bases[k-1])+" "+str(infasta.rnas[1].basepr[paired_bases[k-1]])+" W W A\n")
		    outprm.write("CUTPOINT_CLOSED "+str(infasta.rnas[1].basepr[paired_bases[k-1]]-1) + "\n")


		for n in range(k):
		    paired_bases.pop(0)

	insertstr = ""
	for base in local_regions[j]:
	    if base not in fixed_bases:
		insertstr += str(base) + " "

	outprm.write("ALLOW_INSERT_RES " + insertstr+"\n")

    if use_tertiary:
	for contact in insert_pair_restraints:

	    outprm.write("OBLIGATE PAIR "+str(contact[0])+" "+str(contact[1])+" "+contact[2]+" "+contact[3]+" "+contact[4]+"\n")
	    outprm.write("CUTPOINT_CLOSED " + str(contact[1]-1) + "\n")

    outprm.close()

    os.system("rm "+prefix+".decoy"+str(i)+".out")
    #
    os.system( EXE+" -fasta "+prefix+".final.fasta -nstruct 1 -out::file::silent "+ \
    prefix+".decoy"+str(i)+".out -cycles 0 -params_file " +prefix+".final.prm  -in:file:silent " + file_str \
    + " -close_loops -minimize_rna -chunk_res "+chunkstr)

    j = len(local_regions)-1
    while j >= 0:
	if current_decoy_by_loop[j] < decoys_per_loop:
	    current_decoy_by_loop[j] += 1
	    for k in range(j+1, len(local_regions)):
		current_decoy_by_loop[k] = 1

	    break
	else:
	    j -= 1


