#!/usr/bin/python

import string
from sys import stderr
from os import popen,system
from os.path import exists,dirname,basename,abspath

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
              ' rA': 'a', ' rC': 'c', ' rG': 'g', ' rU': 'u' }


# accepts a pdb file name, returns a string with pdb entries -- or None if there is an error.
def make_rna_rosetta_ready( pdbname, removechain=False, ignore_chain=True, chainids = [], no_renumber = False ):

    fastaid = stderr
    num_model = 0
    #max_model = 60 # for virus
    max_model = 0 # for virus

    # an old baker lab thing:
    #netpdbname = '/net/pdb/' + pdbname[1:3] + '/' + pdbname
    #if not exists(netpdbname):
    netpdbname = pdbname

    outstring = ''

    #print 'Reading ... '+netpdbname
    if not exists( netpdbname ):
        stderr.write( 'DOES NOT EXIST: %s\n' % netpdbname  )
        return None

    if ( netpdbname[-3:] == '.gz' ):
        lines = popen( 'gzcat '+netpdbname ).readlines()
    else:
        lines = open(netpdbname,'r').readlines()

    fastaid.write('>'+pdbname+'\n');

    oldresnum = '   '
    count = 0;

    for i in range( len( chainids ) ) :
        if chainids[i] == '_':
            chainids[i] = ' '

    goodnames = [' rA',' rC',' rG',' rU',' MG']
    #goodnames = [' rA',' rC',' rG',' rU']

    for line in lines:
        if len(line)>5 and line[:6]=='ENDMDL':
            num_model += 1
            if num_model > max_model:  break #Its an NMR model.
        if len(line) <= 21:  continue
        if (line[21] in chainids or ignore_chain):
            line_edit = line

            if line[0:3] == 'TER' and False:
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
            elif (line[0:6] == 'HETATM') & (line[17:20]==' MG'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+' MG'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='OMC'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='5MC'):
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]==' DC'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='CBR'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='CB2'):
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='2MG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='H2U'):
                line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='PSU'):
                line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='5MU'):
                line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='OMG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='7MG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='1MG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='GTP'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]==' YG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='1MA'):
                line_edit = 'ATOM  '+line[6:17]+'  A'+line[20:]

            #Don't save alternative conformations.
            if line[16] == 'B':
                continue;

            if line_edit[0:4] == 'ATOM':
                resnum = line_edit[23:26]
                if not resnum == oldresnum: #  or line_edit[12:16] == ' P  ':
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
                    elif longname == 'rA ':
                        longname =   ' rA'
                    elif longname == 'rC ':
                        longname =   ' rC'
                    elif longname == 'rU ':
                        longname =   ' rU'
                    elif longname == 'rG ':
                        longname =   ' rG'
                    elif longname == 'A  ':
                        longname =   ' rA'
                    elif longname == 'C  ':
                        longname =   ' rC'
                    elif longname == 'U  ':
                        longname =   ' rU'
                    elif longname == 'GUA':
                        longname = ' rG'
                    elif longname == 'ADE':
                        longname = ' rA'
                    elif longname == 'CYT':
                        longname = ' rC'
                    elif longname == 'URA':
                        longname = ' rU'
                    elif longname == 'URI':
                        longname = ' rU'
                    else:
                        if longname not in goodnames:    continue

                    if longer_names.has_key(longname):
                        fastaid.write( longer_names[longname] );
                    else:
                        fastaid.write( 'X')

                    #print "AAH ==> " ,  resnum, oldresnum, line_edit
                    count = count + 1

                oldresnum = resnum

                if not longname in goodnames:
                    continue

                newnum = '%4d' % count
                if no_renumber: newnum = '%4s' % resnum

                line_edit = line_edit[0:16] + ' ' + longname + line_edit[20:22] + \
                            newnum + line_edit[26:]
                if removechain:
                    line_edit = line_edit[0:21]+'  '+line_edit[23:]

                line_edit = line_edit.replace( 'HO2\'', '2HO*' )
                line_edit = line_edit.replace( 'HO5\'', '5HO*' )
                line_edit = line_edit.replace( 'H5\'\'', '2H5*' )

                line_edit = line_edit.replace('\'','*')
                line_edit = line_edit.replace('OP1','O1P')
                line_edit = line_edit.replace('OP2','O2P')

                outstring += line_edit

    fastaid.write('\n')

    return outstring
