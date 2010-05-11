#!/usr/bin/python

def cluster_check( cluster_in ):
    clusterlist = [ 'syd','niau','seth','bes','hapy','apep','gebb','ptah','yah','isis','yah','maat','nut','fin','dig','biox2','biox2_scratch','ade','ade.stanford.edu','steele','steele_scratch','tg-condor','tg-condor_scratch','abe','ncsa','abe_scratch','ade_scratch','vanlang','kwipapat','lovejoy' ];

    cluster = cluster_in
    if cluster not in clusterlist:
        print 'Hey, '+cluster+' is not a known cluster.'
        cluster = ''

    cluster_dir = ''

    if cluster == 'biox2': cluster = 'biox2.stanford.edu'
    if cluster == 'ade': cluster = 'ade.stanford.edu'
    if cluster == 'steele': cluster = 'dasr@tg-steele.purdue.teragrid.org'
    if cluster == 'tg-condor': cluster ='dasr@tg-condor.purdue.teragrid.org'
    if cluster == 'abe': cluster = 'rdas@login-abe.ncsa.teragrid.org'
    if cluster == 'ncsa': cluster ='rdas@tg-login.ncsa.teragrid.org'

    if cluster == 'steele_scratch':
        cluster = 'dasr@tg-steele.purdue.teragrid.org'
        cluster_dir = '/scratch/scratch95/d/dasr/'

    if cluster == 'tg-condor_scratch':
        cluster = 'dasr@tg-condor.purdue.teragrid.org'
        cluster_dir = '/scratch/scratch95/d/dasr/'

    if cluster == 'biox2_scratch':
        cluster = 'biox2.stanford.edu'
        cluster_dir = '/scratch/users/rhiju/'

    if cluster == 'ade_scratch':
        cluster = 'ade.stanford.edu'
        cluster_dir = '/scr/rhiju/'

    if cluster == 'abe_scratch':
        cluster = 'rdas@login-abe.ncsa.teragrid.org'
        cluster_dir = '/scratch/users/rdas/'

    if cluster == 'vanlang':  cluster = 'vanlang@biox2.stanford.edu'
    if cluster == 'kwipapat':  cluster = 'kwipapat@biox2.stanford.edu'
    if cluster == 'lovejoy':  cluster = 'lovejoy@biox2.stanford.edu'

    return (cluster,cluster_dir)

