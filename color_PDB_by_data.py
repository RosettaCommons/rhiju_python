#!/usr/bin/python

from sys import argv

pdb = argv[1]
datafile = argv[2]

offset = 0
if len( argv ) > 3:
    offset = int( argv[3] )

lines = open( datafile ).readlines()
data= {}
avg_data = 0.0
num_data = 0
for line in lines:
    cols = line.split()
    resnum = int( cols[0] )
    dataval = float( cols[1] )
    data[ resnum ] = dataval
    avg_data += dataval
    num_data += 1
if num_data > 0 : avg_data /= num_data

lines = open( pdb ).readlines()

oldresnum = '   '
count = 0;

for line in lines:
    if (len(line)>20): # and (chainid == line[21]):
        line_edit = line
        if line[0:3] == 'TER':
            break
        elif (line[0:6] == 'HETATM') & (line[17:20]=='MSE'): #Selenomethionine
            line_edit = 'ATOM  '+line[6:17]+'MET'+line[20:]
            if (line_edit[12:14] == 'SE'):
                line_edit = line_edit[0:12]+' S'+line_edit[14:]
            if len(line_edit)>75:
                if (line_edit[76:78] == 'SE'):
                    line_edit = line_edit[0:76]+' S'+line_edit[78:]

        if line_edit[0:4] == 'ATOM':
            resnum = int(line_edit[23:26]) - offset

            B_value = avg_data
            if resnum in data.keys():
                B_value = data[ resnum ]

            line_edit = '%s%6.3f%s' % (line_edit[0:61],B_value,line_edit[66:])
            print line_edit[:-1]


