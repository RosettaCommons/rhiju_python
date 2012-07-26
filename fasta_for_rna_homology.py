#!/usr/bin/python

#This class is for reading and processing fasta format sequence alignment of multiple RNAs.  It may also contain a dot-bracket format
#secondary structure common to all the sequences in the file.

import rna, sys

class Fasta_RNA:

    def __init__(self, fastafile):
	print "reading", fastafile
	align = ""
	self.rnas=[]
	self.alignments=[]
	self.labels=[]

	for line in open(fastafile, "r"):
	    if line[0] == ">":  #new alignment
		if len(align) > 1:
		    self.alignments.append(align)
		    align=""
		    self.rnas.append(rna.RNA())
		self.labels.append(line)
	    elif line[0] == "a" or line[0] == "g" or line[0] == "c" or line[0] == "u" or line[0] == "." or line[0] == "-" or line[0] == "(" or line[0] == "[":
		align += line

	if len(align) > 2:
	    self.alignments.append(align)
	    align=""

	for i in range(1, len(self.alignments)):
	    if len(self.alignments[i]) != len(self.alignments[0]):
		print "Error: alignments[",i,"] and alignments[0] are not the same length!"
		sys.exit()

    #Initialize the list of rna objects.  Give them sequences, structures, and find the loops.
    def init_rnas(self):

	#find the secondary structure, if present
	self.structurenum = -1

	for i in range(0, len(self.alignments)):
	    if "." in self.alignments[i] and "(" in self.alignments[i] and ")" in self.alignments[i]:
		self.structurenum = i

	#remove gaps, assign secondary structures
	for i in range(0, len(self.alignments)):
	    if i != self.structurenum:
		seq = ""
		struct = ""
		for j in range(0, len(self.alignments[i])):
		    if self.alignments[i][j] == "a" or self.alignments[i][j] == "g" or self.alignments[i][j] == "c" or self.alignments[i][j] == "u":
			seq += self.alignments[i][j]
			if self.structurenum >= 0:
			    struct += self.alignments[self.structurenum][j]

		self.rnas[i].readseq(seq)

		if self.structurenum >= 0:
		    self.rnas[i].readdb(struct, False)
		    self.rnas[i].findloops()

    #make two lists that list which nucleotide number in the 0th RNA corresponds to which number in the 1st, and vice versa
    def align_by_number(self, i, j):
	self.alignnum0 = []
	self.alignnum1 = []
	tempseq = ""
	targseq = ""
	for k in range(0,len(self.alignments[i])):
	    if self.alignments[i][k] in ["a","c", "g", "u"]:
		tempseq += self.alignments[i][k]
	    if self.alignments[j][k] in ["a","c", "g", "u"]:
		targseq += self.alignments[j][k]
	    if self.alignments[i][k] in ["a","c", "g", "u"]:
		self.alignnum0.append(len(targseq))
	    if self.alignments[j][k] in ["a","c", "g", "u"]:
		self.alignnum1.append(len(tempseq))
	    #set these to be indexed from 1 so that using .index() will return the correct numbers
	self.alignnum0.insert(0,0)
	self.alignnum1.insert(0,0)

    #Find the differences between two alignments.
    def finddifferences(self, i, j):
	self.wc_pair_mutations = []
	self.insertions = []
	self.deletions = []
	self.insertions_plus = []
	self.deletions_plus = []
	self.insertions_minus = []
	self.deletions_in_template = []
	self.loop_mutations = []
	self.template_to_check = []

	tempseq = ""
	targseq = ""

	for k in range(0,len(self.alignments[i])):

	    if self.alignments[i][k] in ["a","c", "g", "u"]:
		tempseq = tempseq + self.alignments[i][k]
	    if self.alignments[j][k] in ["a","c", "g", "u"]:
		targseq += str(self.alignments[j][k])
	    if self.alignments[i][k] != self.alignments[j][k] and self.alignments[i][k] in ["a","c", "g", "u"] and self.alignments[j][k] in ["a","c", "g", "u"] and self.rnas[i].basepr[len(tempseq)] > 0:
		if self.rnas[i].basepr[len(tempseq)] > len(tempseq):
		    self.wc_pair_mutations.append(len(tempseq))
		#this is the special case of a GU or GC pair changing to a GC or GU pair where the G is earlier in the sequence than the pyrimidine
		elif self.rnas[i].basepr[len(tempseq)] > 0 and self.rnas[i].basepr[len(tempseq)] not in self.wc_pair_mutations:
		    self.wc_pair_mutations.append(self.rnas[i].basepr[len(tempseq)])


	    elif self.alignments[i][k] in ["a","c", "g", "u"] and self.alignments[j][k] in [".","-"]:
		if k>0 and k < len(self.alignments[i])-1:
		    if self.rnas[i].basepr[len(tempseq)] > 0 and len(tempseq) < self.rnas[i].length and len(tempseq) > 1:
			if self.rnas[i].basepr[len(tempseq)+1] == self.rnas[i].basepr[len(tempseq)]-1 \
			and self.rnas[i].basepr[len(tempseq)-1] == self.rnas[i].basepr[len(tempseq)]+1 and self.alignments[j][k+1] in ["a","c", "g", "u"] \
			and self.alignments[j][k-1] in ["a","c", "g", "u"]:
			    print "Error. Base pair deletions from the middle of a helix are not supported at this time. Please reconsider your input alignment and see if the deletion might be more fitting at the end of a helix."
			    sys.exit()
		self.deletions.append(len(targseq))
		self.deletions_in_template.append(len(tempseq))
		self.deletions_plus.append(len(targseq)+1)
		self.template_to_check.append(len(tempseq))
		self.template_to_check.append(len(tempseq)+1)
	    elif self.alignments[j][k] in ["a","c", "g", "u"] and self.alignments[i][k] in [".","-"]:
		if k>0 and k < len(self.alignments[i])-1:
		    if self.rnas[j].basepr[len(targseq)] > 0 and len(targseq) < self.rnas[j].length and len(targseq) > 1:
			if self.rnas[j].basepr[len(targseq)] > 0 and self.rnas[j].basepr[len(targseq)+1] == self.rnas[j].basepr[len(targseq)]-1 \
			and self.rnas[j].basepr[len(targseq)-1] == self.rnas[j].basepr[len(targseq)]+1 and self.alignments[i][k+1] in ["a","c", "g", "u"] \
			and self.alignments[i][k-1] in ["a","c", "g", "u"]:
			    print "Error. Base pair insertions in the middle of a helix are not supported at this time. Please reconsider your input alignment and see if the insertion might be more fitting at the end of a helix."
			    print "Base pair between:", len(targseq), "and", self.rnas[j].basepr[len(targseq)], "in the target sequence, column", k+1, "in the alignment."
			    sys.exit()
		self.insertions.append(len(targseq))
		self.insertions_plus.append(len(targseq)-1)
		self.insertions_minus.append(len(targseq)+1)
		self.template_to_check.append(len(tempseq))
		self.template_to_check.append(len(tempseq)+1)
	    elif self.alignments[i][k] != self.alignments[j][k] and self.alignments[i][k] in ["a","c", "g", "u"] \
	    and self.alignments[j][k] in ["a","c", "g", "u"] and self.rnas[i].basepr[len(tempseq)] == 0:
		self.loop_mutations.append(len(targseq))
		self.template_to_check.append(len(tempseq))


