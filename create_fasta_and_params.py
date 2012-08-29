#!/usr/bin/python

from rna_server_conversions import prepare_fasta_and_params_file_from_sequence_and_secstruct
from sys import argv
from os import system
from parse_options import parse_options
from make_tag import make_tag, make_tag_with_dashes
import string

# I could make this a little bit smarter:
# if I look at input_res, I should be able to figure out obligate pairs "on the fly", right? This
#  could either happen here, or within Rosetta.

def Help():
    print
    print argv[0] + " 'sequence in quotes'  'secstruct in quotes'  [tag]"
    print
    print "example:  "+argv[0]+" 'aaa uuu'  '((( )))' test"
    print

if len(argv) < 3:
    Help()
    exit()


sequence = parse_options( argv, "sequence", "" )
secstruct = parse_options( argv, "secstruct", "")
working_res = parse_options( argv, "working_res", [-1] )
offset = parse_options( argv, "offset", 0 )
tag = parse_options( argv, 'tag', 'test' )
fasta = parse_options( argv, 'fasta', "" )
secstruct_file = parse_options( argv, 'secstruct_file', "" )
input_pdbs = parse_options( argv, 's', [""] )
input_res  = parse_options( argv, 'input_res', [-1] )
print argv
fixed_stems = parse_options( argv, 'fixed_stems', False )
print argv
native_pdb = parse_options( argv, 'native', "" )
cst_file = parse_options( argv, 'cst_file', "" )
cutpoint_closed = parse_options( argv, "cutpoint_closed", [-1] )
cutpoint_open = parse_options( argv, "cutpoint_open", [-1] )
extra_minimize_res = parse_options( argv, "extra_minimize_res", [-1] )
virtual_anchor = parse_options( argv, "virtual_anchor", [-1] )
obligate_pair = parse_options( argv, "obligate_pair", [-1] )
remove_pair = parse_options( argv, "remove_pair", [-1] )

print argv
#assert( len( argv ) == 1 )

extra_args = ""
if len( argv ) > 1: extra_args = string.join( argv[1:] )


def is_even( num ):
    return (  2 * (num/2) == num ) # even

def working_res_map( vector, working_res ):
    if len( working_res ) == 0: return vector
    working_vector = []
    if len( working_res ) > 0:
        for m in vector:
            if m in working_res: working_vector.append(  working_res.index( m )+1 )
    return working_vector

if len( fasta ) > 0 :
    assert( len(sequence) == 0 )
    sequence = open( fasta ).readlines()[1][:-1]

if len( secstruct_file ) > 0 :
    assert( len(secstruct) == 0 )
    secstruct = open( secstruct_file ).readlines()[0][:-1]


if len( secstruct ) == 0:  sectruct = len(sequence) * '.'

assert( len( sequence ) == len( secstruct ))
assert( secstruct.count('(') == secstruct.count(')') )
assert( secstruct.count('[') == secstruct.count(']') )


assert( is_even( len(remove_pair) ) )
for m in remove_pair:
    pos = m - offset - 1
    assert( pos > -1 and pos < len( secstruct) and secstruct[ pos ] != '.' )
    secstruct = secstruct[:pos] + '.' + secstruct[ pos+1:]

if len(working_res) <= 1:
    working_sequence = sequence
    working_secstruct = secstruct
else:# slice out the residues.
    working_sequence = ''
    working_secstruct = ''
    for i in range(len(working_res)):
        m = working_res[ i ]
        if i > 0 and m > working_res[ i-1 ] + 1:
            working_sequence  += ' '
            working_secstruct += ' '
        working_sequence  += sequence[  m-1-offset ]
        working_secstruct += secstruct[ m-1-offset ]

print working_sequence
print working_secstruct

assert( working_secstruct.count('(') == working_secstruct.count(')') )
assert( working_secstruct.count('[') == working_secstruct.count(']') )

( fasta_file_outstring, params_file_outstring ) = prepare_fasta_and_params_file_from_sequence_and_secstruct( working_sequence, working_secstruct, fixed_stems )

working_cst_file = ""
cst_file_outstring = ""
if len( cst_file ) > 0: working_cst_file = tag+"_"+cst_file

if len( working_res ) > 0:
    working_input_res = []
    for m in input_res:
        if m not in working_res :
            print m, 'not in working_res: ', working_res
            exit()
        i = working_res.index( m )
        working_input_res.append( i+1 )
    #print
    #print

    # also create a constraints file that will help close chains across 2-5 residue gaps...
    # this is a slight hack, but if it works, might be worth putting a term into Rosetta, as
    # well as automated handling of "working_res"
    cst_gaps = []
    for m in range( len(working_res)-1 ):
        gap = working_res[m+1] - working_res[m]
        if ( gap > 1 and gap < 6 ): cst_gaps.append( [ m+1, m+2, gap, working_res[m], working_res[m+1] ] )
    #print "GAPS TO APPLY CST: ", cst_gaps

    if len( cst_gaps ) > 0:

        if len( working_cst_file ) == 0:  working_cst_file = tag+'_closegaps.cst'
        cst_file_outstring += "[ atompairs ]\n"
        for cst_gap in cst_gaps:
            stdev = 10.0
            gap_length = cst_gap[2]
            max_dist = gap_length * 5.0 + 4
            bonus = 200.0
            cst_file_outstring +=  " O3* %d  C5* %d   FADE %6.3f  %6.3f  %6.3f %6.3f %6.3f \n" % \
                ( cst_gap[0], cst_gap[1],  -stdev, max_dist, stdev, -1*bonus, bonus)

else:
    working_input_res = input_res

if len( cutpoint_closed ) > 0:
    cutpoint_closed = working_res_map( cutpoint_closed, working_res )
    params_file_outstring += "CUTPOINT_CLOSED "+make_tag( cutpoint_closed )+ "\n"

if len( cutpoint_open ) > 0:
    cutpoint_open = working_res_map( cutpoint_open, working_res )
    params_file_outstring += "CUTPOINT_OPEN "+make_tag( cutpoint_open )+ "\n"

if len( virtual_anchor ) > 0:
    virtual_anchor = working_res_map( virtual_anchor, working_res )
    params_file_outstring += "VIRTUAL_ANCHOR "+make_tag( virtual_anchor )+ "\n"


assert( is_even( len(obligate_pair) ) )
if len( obligate_pair ) > 0:
    for m in range( len( obligate_pair)/2 ):
        pos1 = obligate_pair[ 2*m ]
        pos2 = obligate_pair[ 2*m+1 ]
        if len( working_res ) > 0:
            if pos1 not in working_res: continue
            if pos2 not in working_res: continue
            pos1 = working_res.index( pos1 ) + 1
            pos2 = working_res.index( pos2 ) + 1
        params_file_outstring += "OBLIGATE PAIR %d %d W W A \n" % (pos1, pos2)

if len(cst_file) > 0:  # also have input data...
    lines = open( cst_file ).readlines()

    if len( working_res ) > 0:
        for line in lines: # for now, just handle coordinate constraints, later be smart about atompairs
            if len( line ) == 0: continue
            if line[0] == "[":
                assert( line.split()[1] == "coordinates" )
                cst_file_outstring += line
            else:
                cols = line.split()
                if len( cols ) < 4: continue
                if int(cols[1]) not in working_res: continue
                assert( int(cols[3]) in working_res )
                newline = "%s %d %s %d  %s\n" % (cols[0],  working_res.index(int(cols[1]))+1, cols[2],  working_res.index( int(cols[3]) )+1, string.join( cols[4:] ) )
                cst_file_outstring += newline
    else:
        for line in lines:  cst_file_outstring += lines


for i in range(len(input_res)):
    if input_res[i] in input_res[:i]: print 'WARNING -- double input_res?', input_res[i]


fasta_file = tag+'.fasta'
print
print 'Writing to fasta file: ', fasta_file
fid = open( fasta_file, 'w' )
#print fasta_file_outstring
fid.write( fasta_file_outstring )
fid.close()

params_file = tag+'.params'
print 'Writing to params file: ', params_file
fid = open( params_file, 'w' )
#print params_file_outstring
fid.write( params_file_outstring )
fid.close()

if len( working_cst_file ) > 0 :
    cst_file = tag+'.cst'
    print 'Writing to cst file: ', working_cst_file
    fid = open( working_cst_file, 'w' )
    #print cst_file_outstring
    fid.write( cst_file_outstring )
    fid.close()

if ( len(native_pdb) > 0 and len( working_res ) > 0):
    command = "pdbslice.py " + native_pdb + " -subset"
    for m in working_res: command += " %d" % m
    command += " "+tag+"_"
    system( command )
    working_native_pdb = "%s_%s" % (tag,native_pdb)
    print "Writing native to:", working_native_pdb
else:
    working_native_pdb = ""

#########################################
print
print "Sample command line: "

command = "rna_denovo.macosgccrelease  -nstruct 500 -params_file %s -fasta %s  -out:file:silent %s.out  -include_neighbor_base_stacks -minimize_rna" % (params_file, fasta_file, tag )

if len( native_pdb ) > 0:
    if len( working_native_pdb ) > 0:
        command += " -native %s " % working_native_pdb
    else:
        command += " -native %s " % native_pdb

if len( extra_minimize_res ) > 0:
    extra_minimize_res = working_res_map( extra_minimize_res, working_res )
    command += " -extra_minimize_res " + make_tag_with_dashes( extra_minimize_res )

if len( input_pdbs ) > 0:
    command += " -s"
    for pdb in input_pdbs: command += " "+pdb
    assert( len( input_res ) > 0 )

if len( input_res ) > 0:
    command += " -input_res " +make_tag_with_dashes( working_input_res )

if len( working_cst_file ) > 0:
    command += " -cst_file " + working_cst_file

command += ' ' + extra_args

print command

readme = "README_FARFAR"
print "outputting command line to: ", readme
fid = open( readme, 'w' )
fid.write( command + "\n" )
fid.close()

