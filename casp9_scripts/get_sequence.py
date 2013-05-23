#!/usr/bin/python

import string
from os import popen
from os.path import exists,basename
from amino_acids import longer_names

def get_sequence( pdbname, removechain = 0 ):

    netpdbname = pdbname
    assert( exists(netpdbname))
    #print 'Reading ... '+netpdbname

    lines = open(netpdbname,'r').readlines()

    oldresnum = '   '
    count = 0;
    fasta_line = ''
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
                resnum = line_edit[23:26]
                if not resnum == oldresnum:
                    count = count + 1
                    longname = line_edit[17:20]
                    if longer_names.has_key(longname):
                        fasta_line +=  longer_names[longname]
                    else:
                        fasta_line +=  'X'
                oldresnum = resnum

                newnum = '%3d' % count
                line_edit = line_edit[0:23] + newnum + line_edit[26:]
                if removechain:
                    line_edit = line_edit[0:21]+' '+line_edit[22:]

    return fasta_line


def get_annotated_sequence( outfile ):
    lines = popen( "head -n 20 "+outfile ).readlines()
    sequence = ''
    all_variant_types = {}
    for line in lines:
        if line.count( "ANNOTATED_SEQUENCE" ):
            # Need to figure out virtual res
            annotated_sequence = string.split( line[:-1] )[1]
            in_variant_type = 0
            count = -1
            for i in range( len( annotated_sequence) ):
                if ( annotated_sequence[i] == '[' ):
                    in_variant_type = 1
                    variant_type = ''
                elif ( annotated_sequence[i] == ']' ):
                    in_variant_type = 0
                    variant_types.append( variant_type )
                    all_variant_types[ count ] = variant_types
                elif in_variant_type:
                    variant_type += annotated_sequence[i]
                else:
                    count += 1
                    sequence += annotated_sequence[i]
                    variant_types = []
                    all_variant_types[ count ] = variant_types
            break
    return (sequence, all_variant_types )
