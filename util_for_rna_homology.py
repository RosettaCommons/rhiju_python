#!/usr/bin/python

#This class is for useful functions, io, etc. for rna_homology
#

# rhiju -- put this in a util py or something.
def make_tag_with_dashes( int_vector ):
    tag = ''

    start_res = int_vector[0]
    for i in range( 1, len(int_vector)+1 ):
        if i==len( int_vector)  or  int_vector[i] != int_vector[i-1]+1:

            stop_res = int_vector[i-1]
            if stop_res > start_res:
                tag += ' %d-%d' % (start_res, stop_res )
            else:
                tag += ' %d' % (stop_res )

            if ( i < len( int_vector) ): start_res = int_vector[i]

    return tag
