#!/usr/bin/python

import string
from sys import argv,stderr
from os import popen,system
from os.path import exists
from amino_acids import longer_names

assert( len(argv)>2)
pdbname = argv[1]

if (pdbname[-4:] != '.pdb'):
    pdbname += '.pdb'

outfile = pdbname

removechain = 0
if argv.count('-nochain'):
    pos = argv.index('-nochain')
    del( argv[ pos ] )
    removechain = 1

ignore_chain = 0
if argv.count('-ignorechain'):
    pos = argv.index('-ignorechain')
    del( argv[ pos ] )
    ignore_chain = 1

chainids = argv[2:]

netpdbname = '/net/pdb/' + pdbname[1:3] + '/' + pdbname
if not exists(netpdbname):
    netpdbname = pdbname

print 'Reading ... '+netpdbname

lines = open(netpdbname,'r').readlines()

outfile = string.lower( outfile )
outfile = outfile.replace( '.pdb', '_RNA.pdb');
outid = open( outfile, 'w')
print 'Writing ... '+outfile

#fastafile = pdbname[0:4]+chainid+'.pdb.fasta'
#fastaid = open( fastafile, 'w')
fastaid = stderr
#print 'Writing ... '+fastafile
fastaid.write('>'+pdbname+'\n');

oldresnum = '   '
count = 0;

for i in range( len( chainids ) ) :
    if chainids[i] == '_':
        chainids[i] = ' '

goodnames = [' rA',' rC',' rG',' rU']
for line in lines:
    if len(line)>5 and line[:6]=='ENDMDL':break #Its an NMR model.
    if len(line) <= 21:  continue
    if (line[21] in chainids or ignore_chain):
        line_edit = line

        if line[0:3] == 'TER':
            continue
        elif (line[0:6] == 'HETATM') & (line[17:20]=='MSE'): #Selenomethionine
            line_edit = 'ATOM  '+line[6:17]+'MET'+line[20:]
            if (line_edit[12:14] == 'SE'):
                line_edit = line_edit[0:12]+' S'+line_edit[14:]
            if len(line_edit)>75:
                if (line_edit[76:78] == 'SE'):
                    line_edit = line_edit[0:76]+' S'+line_edit[78:]
        elif (line[0:6] == 'HETATM') & (line[17:20]=='5BU'): #Selenomethionine
            line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
        elif (line[0:6] == 'HETATM') & (line[17:20]=='OMC'): #Selenomethionine
            line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
        elif (line[0:6] == 'HETATM') & (line[17:20]==' DC'): #Selenomethionine
            line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
        elif (line[0:6] == 'HETATM') & (line[17:20]=='CBR'): #Selenomethionine
            line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
        elif (line[0:6] == 'HETATM') & (line[17:20]=='CB2'):
            line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]

        #Don't save alternative conformations.
        if line[16] == 'B': continue;

        if line_edit[0:4] == 'ATOM':
            resnum = line_edit[23:26]
            if not resnum == oldresnum or line_edit[12:16] == ' P  ':
                #print "AAH ==> " ,  resnum, oldresnum, line_edit
                count = count + 1
                longname = line_edit[17:20]
                if longname == '  G':
                    longname = ' rG'
                elif longname == '  A':
                    longname =   ' rA'
                elif longname == '  C':
                    longname =   ' rC'
                elif longname == '  U':
                    longname =   ' rU'
                elif longname == 'G  ':
                    longname =   ' rG'
                elif longname == 'A  ':
                    longname =   ' rA'
                elif longname == 'C  ':
                    longname =   ' rC'
                elif longname == 'U  ':
                    longname =   ' rU'
                elif longname == ' DG':
                    longname = ' rG'
                elif longname == ' DA':
                    longname = ' rA'
                elif longname == ' DC':
                    longname = ' rC'
                elif longname == ' DT':
                    longname = ' rU'
                elif longname == 'GUA':
                    longname = ' rG'
                elif longname == 'ADE':
                    longname = ' rA'
                elif longname == 'CYT':
                    longname = ' rC'
                elif longname == 'URA':
                    longname = ' rU'
                else:
                    if longname not in goodnames:    continue

                if longer_names.has_key(longname):
                    fastaid.write( longer_names[longname] );
                else:
                    fastaid.write( 'X')
            oldresnum = resnum


            newnum = '%4d' % count
            line_edit = line_edit[0:16] + ' ' + longname + line_edit[20:22] + \
                        newnum + line_edit[26:]
            if removechain:
                line_edit = line_edit[0:21]+'  '+line_edit[23:]

            line_edit = line_edit.replace('\'','*')
            line_edit = line_edit.replace('OP1','O1P')
            line_edit = line_edit.replace('OP2','O2P')

            outid.write(line_edit)


fastaid.write('\n')
outid.close()
fastaid.close()
