#!/usr/bin/python

import math

class Centerofmass(object):
    def __init__(self, xinit, yinit, zinit):
	self.x = xinit
	self.y = yinit
	self.z = zinit
    def add(self, xadd, yadd, zadd):
	self.x += xadd
	self.y += yadd
	self.z += zadd
    def divide(self, divisor):
	self.x /= divisor
	self.y /= divisor
	self.z /= divisor
    def distance(self, other_com):
	return math.sqrt( (self.x - other_com.x)*(self.x-other_com.x) +  (self.y - other_com.y)*(self.y-other_com.y) +(self.z - other_com.z)*(self.z-other_com.z))

atomicmasses = {"P": 30.97,"OP3": 16.0,"O1P": 16.0,"O2P": 16.0,"O5*": 16.0,"C5*": 12.01,"C4*": 12.01,"O4*": 16.0,"C3*": 12.01,"O3*": 16.0,"C1*": 12.01, \
                "C2*": 12.01,"O2*": 16.0,"N1": 14.007,"C2": 12.01,"N2": 14.007,"N3": 14.007,"C4": 12.01,"C5": 12.01,"C6": 12.01,"O6": 16.0,"N7": 14.007,\
                "C8": 12.01,"N9": 14.007,"1H5*": 1.008,"2H5*": 1.008,"H4*": 1.008 ,"H3*": 1.008 ,"H1*": 1.008,"1H2*": 1.008,"2HO*": 1.008,"H1": 1.008,\
                "1H2": 1.008,"2H2": 1.008,"H8": 1.008,"O2": 16.0,"N4": 14.007,"1H4": 1.008,"2H4": 1.008,"H5": 1.008,"H6": 1.008,"C4": 12.01,"O4": 16.0,\
                "H3": 1.008,"1H6": 1.008,"2H6": 1.008,"N6": 14.007}

class RNA:
    #read the sequence
    def readseq(self, sequence):
	self.seq=["x"]*(len(sequence)+1)
	self.length=len(sequence)
	self.basepr=[0]*(len(sequence)+1)
	self.ntcom=[]

	for i in range(0,len(sequence)+1):
	    self.ntcom.append(Centerofmass(0.0, 0.0, 0.0))

	self.ntmass=[0.0]*(len(sequence)+1)

	for i in range(0,len(sequence)):
	    if sequence[i] == "a" or sequence[i] == "c" or sequence[i] == "g" or sequence[i] == "u":
		self.seq[i+1] = sequence[i]
	    else:
		print "Error: invalid RNA sequence. Nucleotide", i+1, "=", sequence[i]
		break


    #Read a dot bracket notation string into the secondary structure array basepr[]
    def readdb(self, db, remove_noncanonical):
	if len(db) != self.length:
	    print "Error: len(db) = ",len(db),". Length of sequence is", self.length
	    return
	right = [False]*(self.length)
	#find right facing brackets
	for i in range(0,self.length):
	    if db[i] == "(":
		right[i] = True

	#pair with left facing brackets
	for i in range(0,self.length):
	    if db[i] == ")":
		j = i-1
		while right[j] == False and j >= 0:
		    j -= 1
		right[j] = False

		self.basepr[j+1] = i+1
		self.basepr[i+1] = j+1

	#Handle pseudoknots with alternate brackets
	#find right facing brackets
	right = [False]*(self.length)
	for i in range(0,self.length):
	    if db[i] == "[":
		right[i] = True

	#pair with left facing brackets
	for i in range(0,self.length):
	    if db[i] == "]":
		j = i-1
		while right[j] == False and j >= 0:
		    j -= 1
		right[j] = False

		self.basepr[j+1] = i+1
		self.basepr[i+1] = j+1

	#Check for non-canonical pairs:
	if remove_noncanonical:
	    for i in range(1,self.length+1):
		if self.basepr[i] > 0 and self.basepr[i] > i:
		    if not self.canpair(self.seq[i], self.seq[self.basepr[i]]):
			self.basepr[self.basepr[i]] = 0
			self.basepr[i] = 0

    #read a pdb file and compute the centers of mass of each base.
    def readcoms(self, pdbfile):
	currnuc = 0
	prevres = ""
	for line in open(pdbfile, "r"):
	    column = line.split()
	    if len(column) < 6:
		continue
	    if column[5] != prevres:
		prevres = column[5]
		currnuc += 1
	    if column[0] == "ATOM" and len(column) >= 11:
		x = float(column[6])*atomicmasses[column[2]]
		y = float(column[7])*atomicmasses[column[2]]
		z = float(column[8])*atomicmasses[column[2]]
		self.ntcom[int(currnuc)].add(x, y, z)
		self.ntmass[int(currnuc)] += atomicmasses[column[2]]
	for i in range(1, self.length+1):
	    self.ntcom[i].divide(self.ntmass[i])

    #find the loop nucleotides in an RNA
    def findloops(self):
	currentloop = 0
	self.loops = []
	self.loops.append([]) # Ask Matt: Why do we start with a null loop?
	for i in range(1, self.length+1):
	    if self.basepr[i] == 0:
		self.loops[currentloop].append(i)
	    elif self.basepr[i-1] == 0:
		currentloop = currentloop + 1
		self.loops.append([])
        print 'LOOPS: ', self.loops

    #Check for valid base pair
    def canpair(self, base1, base2):
	if base1 not in ["a","c","g","u"]:
	    return False
	if base2 not in ["a","c","g","u"]:
	    return False
	if base1 == "a" and base2 == "u":
	    return True
	elif base1 == "u" and base2 == "a":
	    return True
	elif base1 == "c" and base2 == "g":
	    return True
	elif base1 == "g" and base2 == "c":
	    return True
	elif base1 == "u" and base2 == "g":
	    return True
	elif base1 == "g" and base2 == "u":
	    return True
	else:
	    return False
