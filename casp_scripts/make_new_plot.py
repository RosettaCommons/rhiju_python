#!/usr/bin/python

## this is the worst script in the world

import string
from glob import glob
from os.path import exists
from math import floor,log10,log,exp
from operator import add
from os import chdir,system,popen,getcwd
import sys
from random import random
from math import sqrt

#print 'this script is under construction. tell phil you want to use it\n'
#sys.exit()

MAX_RMS = 15
COURIER_FONT = 12
ss_color = {'H':5,'E':6,'L':7}

if len(sys.argv)<2:
    print '\n'+'-'*75
    print 'usage:',sys.argv[0],' {-gs} {-e} {-ss <ss-prefix>} {-loff} {-f <frag-file>} <contact-file1> {<contact-file2> ... }\n\n\n'
    print '-gs gives grey-scale output (default is rainbow colors)'
    print '-e specifies .eps format (better for combining multiple plots)'
    print '-ss allows you to pass a prefix (ss-prefix) to which the program'
    print '  will add .rdb, .phd, .psipred_ss2, and .jufo_ss in looking for'
    print '  secondary structure information'
    print '-loff turns log scale off'
    print '-f will read the ss-info from the designated frag-file\n'+'-'*75+'\n\n'
    sys.exit()

######### parse args ######################################################
arguments = sys.argv[1:]

# defaults:
file_tag = ''
FASTA = 0    ## should we show the sequence on the plot?
PPO_PLOT = 0 ##           add color line showing phi-psi-omega in native?
ss_base = '' ## base for finding ss-files
frag_file = ''
SSBOND = 0  ## show ss_bond info
GREY_SCALE = 0
RED_SCALE = 0
BG_WHITE = 0
format = 'ps'
SKIP_NEARBY = 5 ## min sequence separation
SS_ONLY = 0
VERBOSE = 0
RESCALE_COLOR = 1
ss_1d = []
SHOW_DIAGONALS = 1
JUST_CONTACTS = 0
ZOOM = 0
LOFF = 0

## check for alternate tags #################

tag_list = ['NS','NC','NC2',
            'DS','DC','DC2',
            'MS',
            'SR','ST']

tag_converter = {'CC':['DC'],
                 'SS':['DS']}  ## mapping from each new tag to list of old tags
narg = len(arguments)
for i in range(narg):
    pos = narg-1-i
    a = arguments[pos]
    if a[0] == '-' and a[1:] in tag_list:
        old_tag = a[1:]
        new_tag = arguments[pos+1]
        if not tag_converter.has_key( new_tag ):tag_converter[new_tag] = []
        tag_converter[new_tag].append( old_tag )
        if old_tag != new_tag: tag_converter[old_tag] = ['XXXXX']
        file_tag = file_tag+'.%s:%s'%(old_tag,new_tag)
        del arguments[pos]
        del arguments[pos]

        if old_tag == 'DC' and new_tag == 'BP':
            tag_converter['MS'] = ['XX']

################################################################

# get the extra 1d scores
extra_1D_info = {}
extra_1D_info_tag_list = [] ## to preserve the order!
while '-1D' in arguments:
    pos = arguments.index('-1D')
    new_tag = arguments[pos+1]
    extra_1D_info[new_tag] = {}
    extra_1D_info_tag_list.append(new_tag)
    del arguments[pos]
    del arguments[pos]

if '-fasta' in arguments:
    pos = arguments.index('-fasta')
    fasta_file = arguments[pos+1]
    if (fasta_file[-3:] == '.gz'):
        sequence= string.join(map(lambda x:string.split(x)[0], popen('gzip -dc '+fasta_file).readlines()[1:]),'')
    else:
        sequence = string.join(map(lambda x:string.split(x)[0], open(fasta_file,'r').readlines()[1:]),'')

    FASTA = 1
    ss_color = {'H':29, 'E':6, 'L':7} ## change the helix color
    del arguments[pos]
    del arguments[pos]

if '-big_bin' in arguments: ## add a color line showing the phi-psi-omega's
    PPO_PLOT = 1
    pos = arguments.index('-big_bin')
    rosetta_pdb = arguments[pos+1]
    assert exists(rosetta_pdb)
    del arguments[pos]
    del arguments[pos]


if '-ss' in arguments:
    pos = arguments.index('-ss')
    ss_base = arguments[pos+1]
    del arguments[pos]
    del arguments[pos]

while '-ss_1d' in arguments:
    pos = arguments.index('-ss_1d')
    ss_1d.append(arguments[pos+1])
    del arguments[pos]
    del arguments[pos]

if '-loff' in arguments:
    pos = arguments.index('-loff')
    LOFF = 1
    del arguments[pos]


if '-f' in arguments:
    pos = arguments.index('-f')
    frag_file = arguments[pos+1]
    del arguments[pos]
    del arguments[pos]

## for disulfides
if '-ssbond' in arguments:
    pos = arguments.index('-ssbond')
    ssbond_file = arguments[pos+1]
    del arguments[pos]
    del arguments[pos]

    lines = map(string.split,open(ssbond_file,'r').readlines())
    disulfides = []
    for line in lines:
        if len(line) == 6 and line[1] == 'CEN':
            pair = map(lambda x:int(x)-1,[line[0],line[2]]) ## numbering starts at 0
            pair.sort()
            print 'disulfide between %d and %d'%(pair[0]+1,pair[1]+1)
            disulfides.append(pair)

    SSBOND = 1

### for handling colors #############################
if '-gs' in arguments:
    pos = arguments.index('-gs')
    del arguments[pos]
    GREY_SCALE = 1

if '-rs' in arguments:
    pos = arguments.index('-rs')
    del arguments[pos]
    RED_SCALE = 1

if '-bg_white' in arguments:
    del arguments[arguments.index('-bg_white')]
    BG_WHITE=1

if '-no_rescale' in arguments:
    del arguments[arguments.index('-no_rescale')]
    RESCALE_COLOR = 0


####### controlling the format #################
if '-L' in arguments:
    pos = arguments.index('-L')
    format = arguments[pos+1]
    del arguments[pos]
    del arguments[pos]
elif '-e' in arguments:
    pos = arguments.index('-e')
    del arguments[pos]
    format = 'eps'
elif '-eps' in arguments:
    pos = arguments.index('-eps')
    del arguments[pos]
    format = 'eps'
elif '-png' in arguments:
    del arguments[ arguments.index('-png') ]
    format = 'png'

## other display parameters

if '-no_diagonals' in arguments:
    del arguments[ arguments.index('-no_diagonals')]
    SHOW_DIAGONALS = 0

if '-no_skip' in arguments:
    del arguments[ arguments.index('-no_skip')]
    SKIP_NEARBY = 0

if '-ss_only' in arguments:
    del arguments[arguments.index('-ss_only')]
    SS_ONLY = 1

if '-v' in arguments:
    del arguments[arguments.index('-v')]
    VERBOSE = 1

if '-just_contacts' in arguments:
    del arguments[ arguments.index('-just_contacts') ]
    JUST_CONTACTS = 1

if '-zoom' in arguments:
    pos = arguments.index('-zoom')
    ZOOM = float( arguments[pos+1] )
    del arguments[pos]
    del arguments[pos]

data_files = arguments


##########################################################################
##########################################################################

def Convert_tag(tag):
    if tag_converter.has_key(tag):
        return tag_converter[tag]
    else:
        return [tag]

def In_range(angle):
    while angle>180: angle = angle-360
    while angle<=-180:angle = angle+360
    return angle

def pp_class(ppo): ## E G A B O
    ppo = ( In_range( ppo[0]), In_range(ppo[1]), In_range(ppo[2]))
    assert -180<=ppo[0]<=180 and -180<=ppo[1]<=180 and -180<=ppo[1]<=180

    if abs( ppo[2] ) < 90: return 'O'
    elif ppo[0]>=0:
        if -100< ppo[1] <= 100:return 'G'
        else: return 'E'
    elif -125 < ppo[1] <= 50: return 'A'
    else: return 'B'

def Read_rosetta_pdb(file):
    lines = map(string.split,
                popen('grep "complete" %s -A10000'%file).readlines()[1:])

    nat_bb = ''
    for line in lines:
        nat_bb = nat_bb + pp_class ( map(float,line[2:5]))
    return nat_bb
##     PP_color = {'A':0.5, #green
##                 'B':0.0, #blue
##                 'G':1.0, #red
##                 'E':0.25, #light blue
##                 'O':0.75} #orange

##     colors = []
##     for line in lines:
##         colors.append ( PP_color[ pp_class ( map(float,line[2:5])) ] )
##     return colors


def Fig(xpos,ypos,marked_color,color,width,marked,out,height=0.0):
    if color == 7 and marked == 0:
        return
    if height==0.0:height=width
    ypos = 10000-ypos
    if marked:
        out.write(string.join(map(str,['2 2 0 1',marked_color,color,
                                       '50 0 20 0.000 0 0 -1 0 0 5\n',
                                       '\t',xpos,ypos,
                                       width+xpos,ypos,
                                       width+xpos,height+ypos,
                                       xpos,height+ypos,
                                       xpos,ypos]))+'\n')
    else:
        out.write(string.join(map(str,['2 2 0 0',color,color,
                                       '50 0 20 0.000 0 0 -1 0 0 5\n',
                                       '\t',xpos,ypos,
                                       width+xpos,ypos,
                                       width+xpos,height+ypos,
                                       xpos,height+ypos,
                                       xpos,ypos]))+'\n')
def Line(xpos,ypos,width,out):
    ypos = 10000-ypos
    out.write(string.join(map(str,['2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 0 0 2\n',
                                   xpos,ypos,xpos+width,ypos]))+'\n')
    return

def RMS100(rmsd,length):
    return rmsd / (1 + 0.5 * (log ( length/100.0)))

def Text(xpos,ypos,label,out,font=0,size=12):
    ypos = 10000-ypos
    out.write(string.join(map(str,['4 0 0 50 0 %d %d 0.0000 4 0 0'%(font,size),
                                   xpos,ypos,label+'\\001']))+'\n')
##     out.write(string.join(map(str,['4 0 0 50 0 %d %d 0.0000 4 195 135'%(font,size),
##                                    xpos,ypos,label+'\\001']))+'\n')
    return

def Color(fraction,LOG):
    if LOFF:
        EXPONENT = 1.0
    else:
        if LOG:
            EXPONENT = log10(2)
        else:
            EXPONENT = 1.0

    if fraction < 0.00001:
        score = fraction
    else:
        score = exp( EXPONENT * log(fraction))

    color = 32 + int(floor(2*16*16*(score/1.00001)))
    return color


def Color_RMSD(rmsd):
    if rmsd > MAX_RMS:
        fraction = 0.0
    else:
        fraction = (float(MAX_RMS) - rmsd) / MAX_RMS
    return Color(fraction, 0) ## linear interpolation between 0 and 15

def Read_phd(phd_file):

    if phd_file[-3:] == '.gz':
        lines = map(string.split,popen('gzip -dc '+phd_file).readlines())
    else:
        lines = map(string.split,open(phd_file,'r').readlines())
    e = []
    h = []
    l = []
    seq= ''

    for line in lines:
        if len(line) == 2 and line[0] == 'AA' and line[1][0] == '|':
            seq = seq + line[1][1:-1]
        elif len(line) == 1 and line[0][:5] == 'prH-|':
            h = h + map(lambda x:float(x)/10,list(line[0][5:]))
        elif len(line) == 1 and line[0][:5] == 'prE-|':
            e = e + map(lambda x:float(x)/10,list(line[0][5:]))
        elif len(line) == 1 and line[0][:5] == 'prL-|':
            l = l + map(lambda x:float(x)/10,list(line[0][5:]))

    seq = seq[:len(h)]

    L = len(seq)
    if len(h) != L or len(e) != L or len(l) != L:
        sys.stderr.write('WARNING: error reading phd file: seq:%d e:%d h:%d l:%d file: %s\n'\
                         %(L,len(e),len(h),len(l),phd_file))
        return [[],[],[],'']
    return [e,h,l,seq]

def Read_fragments(fragment_file):
# 1st five lines of frag file:
# position:            1 neighbors:          200
#
# 1di1 A   229 S L -147.015  136.320  177.891
# 1di1 A   230 A H  -61.103  -28.146  179.117
# 1di1 A   231 V H  -66.832  -46.265  179.757

    base = string.split(fragment_file,'/')[-1]
    if fragment_file[-3:] == '.gz':
        data = popen('gzip -dc '+fragment_file)
        base = base[:-3]
    else:
        data = open(fragment_file,'r')

    size = int(base[-14:-12])

    sys.stderr.write('Reading fragment file: %s size= %d\n'\
                     %(fragment_file,size))

    line = data.readline()
    prev = (-1,-1)
    ss_count = {}
    ppo_count = {}
    while line:
        l = string.split(line)
        line = data.readline()
        assert len(l) == 4 and l[0] == 'position:'

        window = int(l[1])-1 ## numbering starts at 0
        nbrs = int(l[3])
        for i in range(size):
            pos = window + i
            if not ss_count.has_key( pos ):
                ss_count[ pos ] = {'H':0,'E':0,'L':0}
                ppo_count[ pos ] = {'A':0,'B':0,'G':0,'E':0,'O':0}

        for n in range(nbrs):
            for i in range(size):
                line = data.readline()
                l = string.split(line)
                pos = window + i
                ss = l[4]
                ppo = pp_class( map(float, l[5:8] ) )
                ss_count[pos][ss] += 1
                ppo_count[pos][ppo] += 1
            line = data.readline() ## extra blank line

        line = data.readline()
    data.close()

    L = len(ss_count.keys())

    l = []
    e = []
    h = []
    abgeo = ( [], [], [], [], [] )
    ABGEO = 'ABGEO'
    for i in range(L):
        total = reduce(add,ss_count[i].values())
        if total>0:
            e.append( float(ss_count[i]['E'])/total)
            h.append( float(ss_count[i]['H'])/total)
            l.append( float(ss_count[i]['L'])/total)
            for j in range(len(ABGEO)):
                count = ppo_count[i][ABGEO[j]]
                frac = float( count )/total
                if not (frac>= 0 and frac<=1):
                    print i,count,total,frac
                abgeo[j].append( frac )

    frag_ss = (e,h,l)

    return frag_ss,abgeo


##     while line:
##         if len(line)>20 and line[-5] == 'F' and line[-10] == 'P':
##             position = int(line[-9:-6])
##             fragment = int(line[-4:-1])
##             if fragment <= 25:
##                 if (position,fragment) != prev: ## first in fragment
##                     prev = (position,fragment)
##                     count = 0
##                     if fragment == 1: ## initialize counts
##                         if not position%10:
##                             sys.stderr.write('read frag file: position %d\n'%position)

##                         for i in range(9):
##                             fpos = position+i-1 ## numbering starts at 0
##                             if fpos not in ss_count.keys():
##                                 ss_count[fpos] = {'H':0,'E':0,'L':0}

##                 pos = position + count - 1## numbering starts at 0
##                 count = count+1

##                 ss = line[16]
##                 ss_count[pos][ss] = ss_count[pos][ss] + 1
##         line = data.readline()
##     data.close()

def Read_prof(prof_file):
    if (prof_file[-3:] == '.gz'):
        lines = popen('gzip -dc '+prof_file).readlines()
    else:
        lines = open(prof_file,'r').readlines()

    while lines and lines[0][0] == '#':del lines[0]
    del lines[:2]
    lines = map(string.split,lines)
    seq = ''
    e = []
    h = []
    l = []
    for line in lines:
        if len(line) != 5:continue
        seq = seq + line[1]
        e.append(float(line[2]))
        h.append(float(line[3]))
        l.append(float(line[4]))
    return [e,h,l,seq]

def Read_jufo(jufo_file):
    if (jufo_file[-3:]=='.gz'):
        lines = map(string.split,popen('gzip -dc '+jufo_file).readlines())
    else:
        lines = map(string.split,open(jufo_file,'r').readlines())
    l = []
    e = []
    h = []
    seq = ''
    for line in lines:
        if len(line)== 6:
            seq = seq + line[1]
            l.append(float(line[3]))
            h.append(float(line[4]))
            e.append(float(line[5]))
    return [e,h,l,seq]

def Read_jones(jones_file):
    if (jones_file[-3:] == '.gz'):
        lines = map(string.split, popen('gzip -dc '+jones_file).readlines())
    else:
        lines = map(string.split,open(jones_file,'r').readlines())
    l = []
    e = []
    h = []
    seq = ''
    for line in lines:
        if len(line)== 6:
            seq = seq + line[1]
            l.append(float(line[3]))
            h.append(float(line[4]))
            e.append(float(line[5]))
    return [e,h,l,seq]

def Read_rdb(rdb_file):
    if (rdb_file[-3:] == '.gz'):
        lines = popen('gzip -dc '+rdb_file).readlines()
    else:
        lines = open(rdb_file,'r').readlines()
    while lines and lines[0][0] == '#':del lines[0]
    del lines[:2]
    lines = map(string.split,lines)
    seq = ''
    e = []
    h = []
    l = []
    for line in lines:
        if len(line) != 5:continue
        seq = seq + line[1]
        e.append(float(line[2]))
        h.append(float(line[3]))
        l.append(float(line[4]))
    return [e,h,l,seq]

def Read_ss_info(ss_base,ss_1d):
    ss_info = {}

    ## read any 1d files
    for file in ss_1d:
        if not exists(file):
            print 'missing ss_1d file:',file
            continue
        print 'reading ss_1d file:',file

        if (file[-3:] == '.gz'):
            ss = string.split(popen('gzip -dc '+file).readline())[0]
        else:
            ss = string.split(open(file,'r').readline())[0]
        e = []
        h = []
        l = []
        for i in range(len(ss)):
            if ss[i] == 'E':
                e.append(1.0)
                h.append(0.0)
                l.append(0.0)
            elif ss[i] == 'H':
                e.append(0.0)
                h.append(1.0)
                l.append(0.0)
            else:
                e.append(0.0)
                h.append(0.0)
                l.append(1.0)
        ss_info[string.split(file,'/')[-1]] = [e,h,l,''] ## no sequence info

    if ss_base:
        prof_file = ss_base+'.prof_rdb'
        if exists(prof_file):
            print 'reading prof file:',prof_file
            ss_info['prof'] = Read_prof(prof_file)
        else:
            prof_file += '.gz'
            if exists(prof_file):
                print 'reading prof file:',prof_file
                ss_info['prof'] = Read_prof(prof_file)
            else:
                print 'couldnt find prof file:',prof_file

        rdb_file = ss_base+'.rdb'
        if exists(rdb_file):
            print 'reading rdb file:',rdb_file
            ss_info['rdb'] = Read_rdb(rdb_file)
        else:
            rdb_file += '.gz'
            if exists(rdb_file):
                print 'reading rdb file:',rdb_file
                ss_info['rdb'] = Read_rdb(rdb_file)
            else:
                print 'couldnt find rdb file:',rdb_file

        jones_file = ss_base+'.psipred_ss2'

        if exists(jones_file):
            print 'reading jones file:',jones_file
            ss_info['jon'] = Read_jones(jones_file)
        else:
            jones_file += '.gz'
            if exists(jones_file):
                print 'reading jones file:',jones_file
                ss_info['jones'] = Read_jones(jones_file)
            else:
                new_jones_file = string.join(string.split(jones_file[:-3],'/')[:-1],'/')+'/psipred_ss2'
                if exists(new_jones_file):
                    print 'WARNING: couldnt find',jones_file
                    print 'WARNING: Using',new_jones_file,'instead'
                    ss_info['jon'] = Read_jones(new_jones_file)
                else:
                    print 'couldnt find jones file:',jones_file

        jufo_file = ss_base+'.jufo_ss'
        if exists(jufo_file):
            print 'reading jufo file:',jufo_file
            ss_info['juf'] = Read_jufo(jufo_file)
        else:
            jufo_file += '.gz'
            if exists(jufo_file):
                print 'reading jufo file:',jufo_file
                ss_info['jufo'] = Read_jufo(jufo_file)
            else:
                jufo_file = ss_base+'.jufo_1D_ss'
                if exists(jufo_file):
                    print 'reading jufo file:',jufo_file
                    ss_info['juf'] = Read_jufo(jufo_file)
                else:
                    print 'couldnt find jufo file:',jufo_file

        phd_file = ss_base+'.phd'
        if exists(phd_file):
            print 'reading phd file:',phd_file
            ss_info['phd'] = Read_phd(phd_file)
        else:
            phd_file += '.gz'
            if exists(phd_file):
                print 'reading phd file:',phd_file
                ss_info['phd'] = Read_phd(phd_file)
            else:
                print 'couldnt find phd file:',phd_file

    return ss_info

def PPO_plot( decoy_bb, pos1, pos2, width, height, out,
              text_list, text_over, text_down):
    ## get native values
    nat_bb = Read_rosetta_pdb( rosetta_pdb )

    boxed = []
    ABGEO = 'ABGEO'

    if GREY_SCALE:
        MARKED_COLOR = 4
    else:
        MARKED_COLOR = 0

    for i in range(5): ## ABGEO
        for j in range(len( decoy_bb[i] ) ):
            color = Color( decoy_bb[i][j] ,0) ## not log
            if nat_bb[j] != ABGEO[i]:
                Fig(pos1+width*j, pos2, color,color,width,0,out,height)
            else:
                boxed.append( [pos1+width*j, pos2, color ] )
        text_list.append( [text_over, pos2+text_down, ABGEO[i]+':  '] )
        pos2 = pos2 + height

    for b in boxed:
        x = b[0]
        y = b[1]
        color = b[2]
        Fig(x,y,MARKED_COLOR,color,width,1,out,height)

def Color_line(values,pos1,pos2,width,height,out):
    for i in range(len(values)):
        color = Color(values[i],0) ## not logarithmic
        if 1 or values[i]>0:
            Fig(pos1+width*i, pos2, color,color,width,0,out,height)
    return


def Contact_plot(ps_file,contact_file,info,format,LOG=RESCALE_COLOR):
#    ss_color = {'H':5,'E':6,'L':7}


    if GREY_SCALE:
        MARKED_COLOR = 4
    else:
        MARKED_COLOR = 0

    SSBOND_COLOR = 1

    if ss_base or ss_1d:
        ss_info = Read_ss_info(ss_base,ss_1d)
    else:
        ss_info = {}

    if frag_file:
        frag_ss,frag_ppo = Read_fragments(frag_file)
    else:
        frag_ss = ()
        frag_ppo = ()

    if VERBOSE:
        print 'reading file'
    lines = map(string.split,open(contact_file,'r').readlines())
    if exists(contact_file[:-8]+"server_contacts"):
        lines2 = map(string.split,open(contact_file[:-8]+"server_contacts",'r').readlines())

        spps = ps_file.split('.')
        spps.reverse()
        if (spps[1] == 'contacts'):
            spps[1] = 'server_contacts'
        spps.reverse()
        ps_file = ".".join(spps)


    if VERBOSE:
        print 'done reading file'


    native_ss = ''
    decoy_ss = [ [], [], [] ] ############## e, h, l
    decoy_bb = [ [], [], [], [], [] ] ## ABGEO
    ssl = ['E','H','L']

    NC = {}
    NC2 = {}
    DC = {}
    DC2 = {}
    MS = {}
    SR = {}
    ST = {}

    if VERBOSE:
        print 'parsing file'
    for line in lines:
        for tag in Convert_tag( line[0] ):

            if tag == 'NS':
                native_ss = line[1]

            elif tag == 'DS':
                pos = int(line[1])
                assert pos == len(decoy_ss[0])
                for i in range(3):
                    ss = line[2*i+2][0]
                    decoy_ss[ ssl.index(ss) ].append(float(line[2*i+3]))

            elif tag == 'DBB':
                pos = int(line[1])
                ABGEO = ['A','B','G','E','O']
                assert pos == len(decoy_bb[0])
                for i in range(5):
                    bb = line[2*i+2][0]
                    assert ABGEO.count(bb)
                    decoy_bb[ ABGEO.index(bb) ].append(float(line[2*i+3]))

            elif tag == 'NC':
                NC [ (int(line[1]),int(line[2])) ] = 1

            elif tag == 'NC2':
                NC2 [ (int(line[1]),int(line[2])) ] = 1

            elif tag == 'DC':
                #print contact_file[:-8]+"server_contacts"
                if (exists(contact_file[:-8]+'server_contacts')):
                    if (len(lines2) > 0):
                        lin2 = lines2.pop()
                        assert 0 <= float(lin2[3]) <= 1
                        DC [(int(lin2[0]),int(lin2[1]))] = float(lin2[3])
                else:
                    assert 0<= float(line[3]) <= 1
                    DC [(int(line[1]),int(line[2]))] = float(line[3])

            elif tag == 'DC2':
                assert 0<= float(line[3]) <= 1
                DC2 [(int(line[1]),int(line[2]))] = float(line[3])

            elif tag == 'MS':
                assert 0<= float(line[3]) <= 1
                MS [(int(line[2]),int(line[1]))] = float(line[3]) ## reverse i,j

            elif tag == 'SR':
                SR [(int(line[1]),int(line[2]))] = float(line[3])

            elif tag == 'ST':
                ST [(int(line[1]),int(line[2]))] = float(line[3])

            elif extra_1D_info.has_key( tag ):
                pos = int(line[1])
                f = float(line[2])
                extra_1D_info[ tag ] [ pos ] = f

    if VERBOSE:
        print 'done parsing file'
    L = len(decoy_ss[0])

    if native_ss:
        assert len(native_ss) == L
    else:
        native_ss = 'L'*L


    fig_file = '/tmp/phil_junk'+str(random())+'.fig'

    sys.stderr.write('Making %s\n'%ps_file)
    out = open(fig_file,'w')
    out.write('#FIG 3.2\nLandscape\nCenter\nInches\nLetter  \n100.00\nSingle\n-2\n1200 2\n')

    ###########################  define new colors
    l = map(str,range(10))+['a','b','c','d','e','f']

    counter = 32
    if GREY_SCALE:
        for i in range(512):
            red = 255 - i/2
            green = 255 - i/2
            blue = 255 - i/2
            out.write(string.join(map(str,[0,counter,
                                           '#'+l[red/16]+l[red%16]+l[green/16]+\
                                           l[green%16]+l[blue/16]+l[blue%16]]))+'\n')
            counter = counter + 1
    elif RED_SCALE:
        for i in range(512):
            red = 255
            green = 255 - i/2
            blue = 255 - i/2
            out.write(string.join(map(str,[0,counter,
                                           '#'+l[red/16]+l[red%16]+l[green/16]+\
                                           l[green%16]+l[blue/16]+l[blue%16]]))+'\n')
            counter = counter + 1
    else:
        for i in range(256):
            if 1:
##             if LOG:
                green = min(200,(200*i)/128)
##                 green = min(200,2*i)
                blue = min(255,510-2*i)
                red = 0
            else:
                green = 255- i/2
                red = 255-i/2
                blue = 255
            out.write(string.join(map(str,[0,counter,
                                           '#'+l[red/16]+l[red%16]+l[green/16]+\
                                           l[green%16]+l[blue/16]+l[blue%16]]))+'\n')
            counter = counter + 1

        for i in range(256):
            if 1:
##             if LOG:
                red = min(255,2*i)
                green = max(0,min(200,400 - (200*i)/128))
##                 green = min(200,510-2*i)
                blue = 0
            else:
                green = 128-i/2
                red = 128-i/2
                blue = 255
            out.write(string.join(map(str,[0,counter,
                                           '#'+l[red/16]+l[red%16]+l[green/16]+\
                                           l[green%16]+l[blue/16]+l[blue%16]]))+'\n')
            counter = counter + 1


    width = (12500/(2*L + 5))

    if ZOOM:
        width = int( ZOOM * width )

    offset = width*(L+5)

    to_fig = []
    disulf_to_fig = []

    if not SS_ONLY:
        for i in range(L):
            for j in range(L): ## i<j is below the diagonal -- it's flipped!!!


                if abs(j-i) < SKIP_NEARBY: continue
                k = (i,j)

                pos1 = i*width
                pos2 = j*width

                if MS.has_key(k):
                    ## maxsub above the diagonal
                    color = Color(MS[k],LOG)

                    Fig(pos2,pos1,color,color,width,0,out)

                else:
                    if DC.has_key(k):
                        f = DC[k]
                    else:
                        f = 0.0
                    if BG_WHITE and f==0.0:color=7
                    else: color = Color(f,LOG)
                    if SSBOND and [min(i,j),max(i,j)] in disulfides:
                        disulf_to_fig.append([pos2,pos1,SSBOND_COLOR,color,width,1,out])
                    elif NC.has_key(k):
                        to_fig.append([pos2,pos1,MARKED_COLOR,color,width,1,out])
                    else:
                        Fig(pos2,pos1,color,color,width,0,out)

                if DC2: ## for a second contact-type plot
                    if DC2.has_key(k):
                        f = DC2[k]
                    else:
                        f = 0.0
                    if BG_WHITE and f==0.0: color = 7
                    else: color = Color(f,LOG)
                    if NC2.has_key(k):
                        to_fig.append([pos2+offset,pos1,MARKED_COLOR,color,width,1,out])
                    else:
                        Fig(pos2+offset,pos1,color,color,width,0,out)

                elif ST: ## subclustering
                    ## subcluster rmsd above diagonal, offset
                    if SR.has_key(k):
                        color = Color_RMSD ( SR[k])
                        Fig (pos1+offset,pos2,color,color,width,0,out)

                    if ST.has_key(k):
                        color = Color_RMSD (ST[k])
                        Fig (pos2+offset,pos1,color,color,width,0,out)


                ## sequence position marks along the top:
                if (k[1] == L-1 and not k[0]%10):
                    Fig(pos1,pos2+width,0,0,width,0,out)
                elif (k[1] == L-1 and not k[0]%5):
                    Fig(pos1,pos2+width,0,0,width/2,0,out)


        for tf in to_fig+disulf_to_fig:
            Fig(tf[0],tf[1],tf[2],tf[3],tf[4],tf[5],tf[6])


    text_list = [[0,(L+2)*width,info]]

    if not JUST_CONTACTS:
        ## make secondary structure plot ####################33

        #calculate deviations of decoys from predictions
        dev = [ [0.5]*L, [0.5]*L ] # 0=E,1=H
        ks =ss_info.keys()
        for type in ks:
            if (len(ss_info[type][0]) !=L):
                sys.stderr.write('WARNING: prediction/sequence mismatch: %s\n'%type);
                del ss_info[type]
        if ss_info:
            for pos in range(L):
                for ss in range(2):
                    pred = 0.0
                    for type in ss_info.keys():
                        pred = pred + ss_info[type][ss][pos]
                    pred = pred/len(ss_info.keys())
                    dev[ss] [pos] = max(0.0,min(1.0, 0.5 + (decoy_ss[ss][pos]-pred)))

        ## smooth dev
        for i in range(2):
            for pos in range(L):
                if pos==0:
                    dev[i][pos] = ( dev[i][pos]+dev[i][pos+1])/2
                elif pos==L-1:
                    dev[i][pos] = ( dev[i][pos]+dev[i][pos-1])/2
                else:
                    dev[i][pos] = ( dev[i][pos]+dev[i][pos-1]+dev[i][pos+1])/3

        assert len(decoy_ss[0]) == len(native_ss)
        native = [[],[],[]]
        for i in range(L):
            if native_ss[i]=='E':
                native[0].append(1.0)
                native[1].append(0.0)
                native[2].append(0.0)
            elif native_ss[i]=='H':
                native[0].append(0.0)
                native[1].append(1.0)
                native[2].append(0.0)
            elif native_ss[i]=='L':
                native[0].append(0.0)
                native[1].append(0.0)
                native[2].append(1.0)

        if frag_ss:
            if len(frag_ss[0]) != L:
                sys.stderr.write('WARNING: fragment file length mismatch\n')
                frag_ss = ()
                frag_ppo = ()

        pos2 = L*width+600

        text_down = -125
        text_over = -500
        line_up = 3
        DW=200
        PW=150
        for ss in range(2): ## 0=strand, 1=helix
            ss_name = 'EH'[ss]
            ## decoys:
            text_list.append([text_over,pos2+text_down,'dec: '+ss_name])
            Color_line(decoy_ss[ss],0,pos2,width,DW,out)
            Line(0,pos2+line_up,width*L,out)
            if DC2:
                Color_line(decoy_ss[ss],offset,pos2,width,DW,out)
                Line(offset,pos2+line_up,width*L,out)

            pos2 = pos2+PW+10

            ## fragments
            if frag_ss:
                text_list.append([text_over,pos2+text_down,'frg: '+ss_name])
                Color_line(frag_ss[ss],0,pos2,width,PW,out)
                Line(0,pos2+line_up,width*L,out)
                pos2 = pos2+PW+10


            ## predictions
            for type in ss_info.keys():
                if len(ss_info[type][ss]) == L:
                    text_list.append([text_over,pos2+text_down,type+': '+ss_name])
                    Color_line(ss_info[type][ss],0,pos2,width,PW,out)
                    pos2 = pos2+PW+10

            ## native
            Line(0,pos2-PW-10+line_up,width*L,out)
            text_list.append([text_over,pos2+text_down,'nat: '+ss_name])
            Color_line(native[ss],0,pos2,width,PW,out)
            if DC2:
                Color_line(native[ss],offset,pos2,width,PW,out)
            pos2 = pos2+DW+50

        if ss_info:
            Color_line(dev[0],0,pos2,width,DW,out)
            text_list.append([text_over,pos2+text_down,'dev: E'])
            pos2 = pos2+DW
            Color_line(dev[1],0,pos2,width,DW,out)
            text_list.append([text_over,pos2+text_down,'dev: H'])
        if PPO_PLOT:
            PPO_plot( decoy_bb, 0, pos2, width, DW, out, text_list, text_over,
                      text_down )
            pos2 = pos2 + 5*DW
            if frag_ppo:
                PPO_plot( frag_ppo, 0, pos2, width, DW, out, text_list, text_over,
                          text_down )
                pos2 = pos2 + 5*DW

            #pp_colors = Read_rosetta_pdb(rosetta_pdb) [1:-1] ## ignore 1st and last
            #pos2 = pos2+DW
            #Color_line(pp_colors,width,pos2,width,DW,out)
            #text_list.append([text_over,pos2+text_down,'ppo:  '])
        for tag in extra_1D_info_tag_list:
            values = []
            extra = extra_1D_info[tag]
            for pos in range(L):
                if extra.has_key(pos):
                    values.append( extra[pos] )
            if len(values) != L:
                print 'problem with extra_1D_info for tag:',tag
                continue
            pos2 = pos2+DW
            Color_line(values,0,pos2,width,DW,out)
            text_list.append([text_over,pos2+text_down,tag])


        ## consensus decoy ss
        pred_ss = ''
        for pos in range(L):
            if decoy_ss[0][pos] > decoy_ss[1][pos] and decoy_ss[0][pos] > decoy_ss[2][pos]:
                pred_ss = pred_ss + 'E'
            elif decoy_ss[1][pos] > decoy_ss[2][pos]:
                pred_ss = pred_ss + 'H'
            else:
                pred_ss = pred_ss + 'L'

        assert len(pred_ss) == L


        if not SS_ONLY:
            for i in range(len(pred_ss)):
                ss = pred_ss[i]

                ## boundaries of first square:
                Fig(i*width,(L+1)*width,ss_color[ss],ss_color[ss],width,0,out)
                Fig(-1*width,i*width,ss_color[ss],ss_color[ss],width,0,out)
                Fig(i*width,-1*width,ss_color[native_ss[i]],ss_color[native_ss[i]],
                    width,0,out)
                Fig((L+1)*width,i*width,ss_color[native_ss[i]],
                    ss_color[native_ss[i]],width,0,out)

                ## diagonals:
                if SHOW_DIAGONALS:
                    Fig(i*width,i*width,ss_color[ss],ss_color[ss],width,0,out)
                    if ST or DC2:
                        Fig(i*width+offset,i*width,ss_color[ss],ss_color[ss],width,0,out)

                if DC2: ## boundaries of second square
                    Fig(offset+i*width,(L+1)*width,ss_color[ss],
                        ss_color[ss],width,0,out)
                    Fig(offset-1*width,i*width,ss_color[ss],
                        ss_color[ss],width,0,out)
                    Fig(offset+i*width,-1*width,ss_color[native_ss[i]],
                        ss_color[native_ss[i]],width,0,out)
                    Fig(offset+(L+1)*width,i*width,ss_color[native_ss[i]],
                        ss_color[native_ss[i]],width,0,out)

                if FASTA: ## add sequence
                    size = min(12,max(5,width/10))
                    text_list.append( [i*width,(i-1)*width,sequence[i],COURIER_FONT,size])
                    text_list.append( [i*width,L*width,sequence[i],COURIER_FONT,size])
                    text_list.append( [(L+1)*width,(i-1)*width,sequence[i],COURIER_FONT,size])
                    text_list.append( [i*width,-2*width,sequence[i],COURIER_FONT,size])
                    text_list.append( [-1*width,(i-1)*width,sequence[i],COURIER_FONT,size])

    if not SS_ONLY:
        ## make key #######################3
        show = [1000,500,100,50,10,5,1]
        labels = {}
        for s in show:
            cf = float(s)/1000
            labels[Color(cf,LOG)] = str(s)


        w2 = 5000/512
        for i in range(512):
            color = 32+i

            pos1 = -150
            pos2 = i*w2

            for j in range(10):
                Fig(pos1+j*w2,pos2,color,color,w2,0,out)

            if color in labels.keys():
                for j in range(10):
                    Fig(pos1-1-j*w2,pos2,color,color,w2,0,out)

                text_list.append([pos1-500,pos2,labels[color]])

        ## make key for contacts and maxsub #######################3
        show = [1000,500,100,50,10,5,1]
        labels = {}
        for s in show:
            cf = float(s)/1000
            labels[Color(cf,LOG)] = str(s)


        w2 = 5000/512
        for i in range(512):
            color = 32+i

            pos1 = -150
            pos2 = i*w2

            for j in range(10):
                Fig(pos1+j*w2,pos2,color,color,w2,0,out)

            if color in labels.keys():
                for j in range(10):
                    Fig(pos1-1-j*w2,pos2,color,color,w2,0,out)

                text_list.append([pos1-500,pos2,labels[color]])


    if not SS_ONLY and (ST and not DC2): ## make key for rmsd100 #######################3

        lengths = [25,50,100,150]
        rmsds = range(1,15)
        labels  = {}
        for r in rmsds:
            for l in lengths:
                if RMS100(r,l) >= MAX_RMS:
                    continue
                labels [ Color_RMSD (RMS100(r,l)) ] = [r,l]

        for l in lengths:
            pos1 = offset + 10*width - 10*w2
            pos2 = L*width + 2*width

            up = lengths.index(l)
            text_list.append([pos1,pos2+10*w2+200*up,str(l)])



        w2 = 5000/512
        for i in range(512):
            color = 32+i

            pos1 = offset + 10*width + (511-i)*w2
            pos2 = L*width + 2*width

            for j in range(10):
                Fig(pos1,pos2-j*w2,color,color,w2,0,out)

            if color in labels.keys():
                for j in range(10):
                    Fig(pos1,pos2+1+j*w2,color,color,w2,0,out)

                up = lengths.index(labels[color][1])
                text_list.append([pos1,pos2+10*w2+200*up,str(labels[color][0])])


    ## write all the text
    for tl in text_list:
        if len(tl) == 3:
            Text(tl[0],tl[1],tl[2],out)
        else:
            font = tl[3]
            size = tl[4]
            Text(tl[0],tl[1],tl[2],out,font,size)


    out.close()

    command = 'fig2dev -L %s -l 0 %s %s'\
              %(format,fig_file,ps_file)
    if VERBOSE:
        print command
    system(command)
    if VERBOSE:
        print 'done:',command

    #print 'fig_file:',fig_file
    system('rm '+fig_file)

    return


for file in data_files:
    if file[0] == '/':
        info = file
    else:
        info = getcwd()
        if info[-1] != '/':info = info +'/'
        info = info + file

    Contact_plot(file+file_tag+'.'+format,file,info,format)
