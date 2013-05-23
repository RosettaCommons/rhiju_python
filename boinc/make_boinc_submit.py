#!/usr/bin/python
##
## make boinc submit files for homologs
##

from phil import *
import sys
from os.path import dirname, basename, abspath
import string

NSTRUCT = 30
QUEUE = 10

select_mode = 0

mode='-shortrelax'
if argv.count(mode):
    commandlinedescription = 'CASP7_ABRELAX_SHORTRELAX_SAVE_ALL_OUT'
    args = '-output_silent_gz -silent -increase_cycles 10  -new_centroid_packing -abrelax  -output_chi_silent -stringent_relax -vary_omega -omega_weight 0.5 -farlx -ex1 -ex2 -termini -short_range_hb_weight 0.50 -long_range_hb_weight 1.0 -no_filters -rg_reweight 0.5 -rsd_wt_helix 0.5 -rsd_wt_loop 0.5 -output_all -accept_all -nstruct %d -relax_score_filter' % NSTRUCT
    descriptionfile = '/work/rhiju/CASP7/casp7.description.shorter.txt';
    select_mode = 1
    pos = argv.index(mode); del(argv[pos])

mode='-abrelax'
if argv.count(mode):
    commandlinedescription = 'CASP7_ABRELAX_SAVE_ALL_OUT'
    args = '-output_silent_gz -silent -increase_cycles 10  -new_centroid_packing -abrelax  -output_chi_silent -stringent_relax -vary_omega -omega_weight 0.5 -farlx -ex1 -ex2 -termini -short_range_hb_weight 0.50 -long_range_hb_weight 1.0 -no_filters -rg_reweight 0.5 -rsd_wt_helix 0.5 -rsd_wt_loop 0.5 -output_all -accept_all -do_farlx_checkpointing -nstruct %d -relax_score_filter' % NSTRUCT
    descriptionfile = '/work/rhiju/CASP7/casp7.description.shorter.txt';
    select_mode = 1
    pos = argv.index(mode); del(argv[pos])

mode='-abinitio'
if argv.count(mode):
    commandlinedescription = 'CASP7_ABINITIO_SAVE_ALL_OUT'
    args = '-output_silent_gz -silent -increase_cycles 10  -new_centroid_packing -termini -rg_reweight 0.5 -rsd_wt_helix 0.5 -rsd_wt_loop 0.5 -output_all -accept_all -nstruct %d ' % NSTRUCT
    descriptionfile = '/work/rhiju/CASP7/casp7.description.shorter.txt';
    select_mode = 1
    pos = argv.index(mode); del(argv[pos])

jump_code = 0
mode='-jumpcode'
if argv.count(mode):
    commandlinedescription = 'CASP7_JUMPABINITIO_SAVE_ALL_OUT'
#    args = '-output_silent_gz -silent -increase_cycles 10  -new_centroid_packing -termini -rg_reweight 0.5 -rsd_wt_helix 0.5 -rsd_wt_loop 0.5 -output_all -accept_all -nstruct %d -jumping -sheet_from_barcode -filter_jumps_with_bonus -ignore_sspair_barcode_in_score -handedness_score -close_chainbreaks ' % NSTRUCT
    args = '-output_silent_gz -silent -increase_cycles 10  -new_centroid_packing -termini -rg_reweight 0.5 -rsd_wt_helix 0.5 -rsd_wt_loop 0.5 -output_all -accept_all -nstruct %d -jumping -sheet_from_barcode -filter_jumps_with_bonus -ignore_sspair_barcode_in_score -handedness_score ' % NSTRUCT
    descriptionfile = '/work/rhiju/CASP7/casp7.description.shorter.txt';
    select_mode = 1
    pos = argv.index(mode); del(argv[pos])
    jump_code = 1


jump1 = 0
mode='-jump1'
if argv.count(mode):
    commandlinedescription = 'CASP7_JUMP_ONEPAIR_SAVE_ALL_OUT'
    args = '-output_silent_gz -silent -increase_cycles 10  -new_centroid_packing -termini -rg_reweight 0.5 -rsd_wt_helix 0.5 -rsd_wt_loop 0.5 -output_all -accept_all -nstruct %d -jumping -sheet1 2 ' % NSTRUCT
    descriptionfile = '/work/rhiju/CASP7/casp7.description.shorter.txt';
    select_mode = 1
    pos = argv.index(mode); del(argv[pos])

jump2 = 0
mode='-jump2'
if argv.count(mode):
    commandlinedescription = 'CASP7_JUMP_TWOPAIR_SAVE_ALL_OUT'
    args = '-output_silent_gz -silent -increase_cycles 10  -new_centroid_packing -termini -rg_reweight 0.5 -rsd_wt_helix 0.5 -rsd_wt_loop 0.5 -output_all -accept_all -nstruct %d -jumping -sheet1 2 -sheet2 2 ' % NSTRUCT
    descriptionfile = '/work/rhiju/CASP7/casp7.description.shorter.txt';
    select_mode = 1
    pos = argv.index(mode); del(argv[pos])

mode='-jumprelax'
if argv.count(mode):
    commandlinedescription = 'CASP7_JUMPRELAX_SAVE_ALL_OUT'
    args = '-output_silent_gz -silent -increase_cycles 10  -new_centroid_packing -output_chi_silent -stringent_relax -vary_omega -omega_weight 0.5 -farlx -ex1 -ex2 -termini -short_range_hb_weight 0.50 -long_range_hb_weight 1.0 -no_filters -rg_reweight 0.5 -rsd_wt_helix 0.5 -rsd_wt_loop 0.5 -output_all -accept_all -nstruct %d -jumping -sheet_from_barcode -filter_jumps_with_bonus -ignore_sspair_barcode_in_score -handedness_score -jump_relax  -relax_score_filter' % NSTRUCT
    descriptionfile = '/work/rhiju/CASP7/casp7.description.shorter.txt';
    jump_code = 1
    select_mode = 1
    pos = argv.index(mode); del(argv[pos])


mode='-jump1relax'
if argv.count(mode):
    commandlinedescription = 'CASP7_JUMP_ONEPAIR_RELAX_SAVE_ALL_OUT'
    args = '-output_silent_gz -silent -increase_cycles 10  -new_centroid_packing -output_chi_silent -stringent_relax -vary_omega -omega_weight 0.5 -farlx -ex1 -ex2 -termini -short_range_hb_weight 0.50 -long_range_hb_weight 1.0 -no_filters -rg_reweight 0.5 -rsd_wt_helix 0.5 -rsd_wt_loop 0.5 -output_all -accept_all -nstruct %d -jumping -sheet1 2 -handedness_score -jump_relax  -relax_score_filter' % NSTRUCT
    descriptionfile = '/work/rhiju/CASP7/casp7.description.shorter.txt';
    select_mode = 1
    pos = argv.index(mode); del(argv[pos])

put_in_filters = 0
if argv.count('-filter'):
    i = argv.index('-filter')
    del argv[i]
    put_in_filters = 1

if argv.count('-no_termini'):
    i = argv.index('-no_termini')
    del argv[i]
    arguments = arguments.replace('-termini','')

put_in_barcode = 0
if argv.count('-barcode'):
    i = argv.index('-barcode')
    del argv[i]
    barcode_file_first = abspath(argv[i])
    del argv[i]
    put_in_barcode = 1
    defined_barcode_prefix = 0
    commandlinedescription += '_BARCODE'

put_in_constraint = 0
if argv.count('-constraint_file'):
    i = argv.index('-constraint_file')
    del argv[i]
    constraint_file_first = abspath(argv[i])
    del argv[i]
    put_in_constraint = 1
    defined_constraint_prefix = 0
    commandlinedescription += '_CONSTRAINT'


if argv.count('-no_relax_score_filter'):
    i = argv.index('-no_relax_score_filter')
    del argv[i]
    args = args.replace('-relax_score_filter','')

put_in_contact = 0
if argv.count('-contact_file'):
    i = argv.index('-contact_file')
    del argv[i]
    contact_file_first = abspath(argv[i])
    del argv[i]
    put_in_contact = 1
    defined_contact_prefix = 0
    commandlinedescription += '_CONTACT'
    args  += ' -score_contact_weight 0.5 -score_contact_readindist -score_contact_threshold 0.0 -score_contact_flag'

if argv.count('-score_contact_weight'):
    i = argv.index('-score_contact_weight')
    del argv[i]
    score_contact_weight = float(argv[i])
    del argv[i]
    args = args.replace('-score_contact_weight 0.5','-score_contact_weight %3.1f' % score_contact_weight)


user_tag = ''
if argv.count('-tag'):
    i = argv.index('-tag')
    del argv[i]
    user_tag = argv[i]
    del argv[i]
    commandlinedescription = commandlinedescription.replace('CASP7_','CASP7_'+user_tag)

put_in_pairing_file = 0
if argv.count('-pairing_file'):
    i = argv.index('-pairing_file')
    del argv[i]
    pairing_file_first = abspath(argv[i])
    del argv[i]
    put_in_pairing_file = 1
    defined_pairing_file_prefix = 0

if jump_code:
    assert(put_in_barcode)
    assert(put_in_pairing_file)

if argv.count('-queue'):
    i = argv.index('-queue')
    del argv[i]
    QUEUE = int(argv[i])
    del argv[i]

set_priority = 0
if argv.count('-priority'):
    i = argv.index('-priority')
    del argv[i]
    PRIORITY = int(argv[i])
    del argv[i]
    set_priority = 1

if argv.count('-increase_cycles'):
    i = argv.index('-increase_cycles')
    del argv[i]
    increase_cycles = float(argv[i])
    del argv[i]
    args = args.replace('-increase_cycles 10','-increase_cycles %3.1f' % increase_cycles)

filter1_defined = 0
if argv.count('-filter1'):
    i = argv.index('-filter1')
    del argv[i]
    filter1_input = float(argv[i])
    del argv[i]
    filter1_defined = 1

fasta_list = argv[1:]

outfile = sys.stdout;


fasta_list.sort()

for fasta in fasta_list:
    prefix_defined = 0
    gzipped = 0

    frag_dir = dirname(abspath(fasta))+'/'
    fasta = basename(fasta)

    if fasta[-3:] == '.gz' :
        fasta = fasta[:-3]
        gzipped = 1

    assert(gzipped == 1)

    fivelettercode = fasta[-11:-6]
    fourlettercode = fasta[-11:-7]
    chaincode = fasta[-7:-6]

    prefix = fasta[:-11]
    if len(prefix)>0:
       prefix_defined = 1

    if put_in_barcode and not defined_barcode_prefix:
        defined_barcode_prefix = 1
        barcode_prefix_start = barcode_file_first.find(prefix)
        barcode_prefix_end = barcode_prefix_start + len(prefix)
        assert(barcode_prefix_start >= 0)


    if put_in_constraint and not defined_constraint_prefix:
        defined_constraint_prefix = 1
        constraint_prefix_start = constraint_file_first.find(prefix)
        constraint_prefix_end = constraint_prefix_start + len(prefix)
        assert(constraint_prefix_start >= 0)

    if put_in_contact and not defined_contact_prefix:
        defined_contact_prefix = 1
        contact_prefix_start = contact_file_first.find(prefix)
        contact_prefix_end = contact_prefix_start + len(prefix)
        assert(contact_prefix_start >= 0)


    if put_in_pairing_file and not defined_pairing_file_prefix:
        defined_pairing_file_prefix = 1
        pairing_file_prefix_start = pairing_file_first.find(prefix)
        pairing_file_prefix_end = pairing_file_prefix_start + len(prefix)
        assert(pairing_file_prefix_start >= 0)


    outfile.write('name = '+fivelettercode+'_'+commandlinedescription+'_'+prefix+'\n')
    if prefix_defined:
        outfile.write('description = homolog '+prefix+' for '+fivelettercode+' with '+commandlinedescription+' command line\n')
    else:
        outfile.write('description = '+fivelettercode+' with '+commandlinedescription+' command line\n')

    outfile.write('inputfiles = ')

    # Here come the input files. Run asserts to make sure they all exist!
    for suffix in ['fasta.gz','psipred_ss2.gz']:
        new_file = '%s/%s%s.%s'%(frag_dir,prefix,fivelettercode,suffix)
        #print new_file
        assert(exists(new_file))
        outfile.write(new_file+';')

    new_frag_tag = "_05.200_v1_3"
    for m in [3,9]:
        new_file = '%s/%saa%s%02d%s.gz'%(frag_dir,'boinc_'+prefix,fivelettercode,m,new_frag_tag)
        assert(exists(new_file))
        outfile.write(new_file+';')

    for suffix in ['pdb.gz']:
        new_file = '%s/%s.%s'%(frag_dir,fourlettercode,suffix)
        if (exists(new_file)):
            outfile.write(new_file+';')

    if (put_in_barcode):
        barcode_file = barcode_file_first[:barcode_prefix_start] + prefix + barcode_file_first[barcode_prefix_end:]
        assert( exists(barcode_file))
        outfile.write(barcode_file+';')

    if (put_in_constraint):
        constraint_file = constraint_file_first[:constraint_prefix_start] + prefix + constraint_file_first[constraint_prefix_end:]
        assert( exists(constraint_file))
        outfile.write(constraint_file+';')

    if (put_in_contact):
        contact_file = contact_file_first[:contact_prefix_start] + prefix + contact_file_first[contact_prefix_end:]
        assert( exists(contact_file))
        outfile.write(contact_file+';')



    if (put_in_pairing_file):
        pairing_file = pairing_file_first[:pairing_file_prefix_start] + prefix + pairing_file_first[pairing_file_prefix_end:]
        assert( exists(pairing_file))
        outfile.write(pairing_file+';')


#    if prefix == 'hom001_':
#        new_file = '/work/rhiju/CASP7/bioinfo/outfiles/%s.3djury.out.gz'%(fivelettercode)
#        assert (exists(new_file))
#        outfile.write(new_file+';')

    outfile.write('%s' % descriptionfile)
    outfile.write('\n')

    if prefix_defined:
        args_withprefix = args+' -protein_name_prefix '+prefix+' -frags_name_prefix '+'boinc_'+prefix
    else:
        args_withprefix = args+' -frags_name_prefix boinc_'


    assert( exists(descriptionfile))
    args_withprefix += ' -description_file '+ basename(descriptionfile)


    if prefix == 'hom001_':
        new_file = '%s.3djury.out'%(fivelettercode)
        args_withprefix += ' -server_models ' + new_file

    if (shortrelax):
        args_withprefix += ' -filter1 -9999 -filter2 -9999'

    if (put_in_barcode):
        barcode_file = basename(barcode_file)
        if barcode_file[-3:] == '.gz':
            barcode_file = barcode_file[:-3]
        args_withprefix += ' -barcode_mode 3 -barcode_file %s -output_flavor' % barcode_file

    if (put_in_contact):
        contact_file = basename(contact_file)
        if contact_file[-3:] == '.gz':
            contact_file = contact_file[:-3]
#        args_withprefix += ' -score_contact_file %s -score_contact_weight 0.5 -score_contact_calpha -score_contact_readindist -score_contact_threshold 0.0 -score_contact_flag' % contact_file
        args_withprefix += ' -score_contact_file %s ' % contact_file


    if (put_in_pairing_file):
        pairing_file = basename(pairing_file)
        if pairing_file[-3:] == '.gz':
            pairing_file = pairing_file[:-3]
        args_withprefix += ' -pairing_file %s' % pairing_file


    if (put_in_filters):
        filterfile = 'filters/'+ prefix+fivelettercode+'filters.txt'
        if not exists(filterfile):
            filterfile =  prefix+fivelettercode+'filters.txt'
        if exists(filterfile):
            filters = open( filterfile ) .readlines()
            filter1 = int(float(string.split(filters[0])[0]))
            filter2 = int(float(string.split(filters[0])[1]))
            if filter1_defined:
                filter1 = filter1_input
            args_withprefix += ' -filter1 %d  -filter2 %d' % (filter1, filter2)
        else:
            print '************  '+filterfile+' NOT FOUND'
            sys.exit()

    outfile.write('arguments = xx '+fourlettercode+' '+chaincode+' '+args_withprefix+'\n')

    resultfile = 'xx'+fourlettercode+'.out.gz'
    outfile.write('resultfiles = '+resultfile+'\n')

    if set_priority:
        outfile.write('priority = %d\n' % PRIORITY)

    outfile.write('Queue = %d\n\n' % QUEUE)

outfile.close()
