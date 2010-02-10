#!/usr/bin/python

################################################################
def parse_options( argv, tag, default):
    value = default
    if argv.count( "-"+tag ):
        pos = argv.index( "-"+tag )

        if ( default == 0 and
             ( pos == (len( argv )-1) or
               argv[ pos+1 ][0] == '-' ) ): # Just a boolean
            value = 1
        elif( isinstance( default, list ) ):
            value = []
            offset = 1
            while ( (pos + offset) < len (argv) and not (argv[ pos+offset ][0] == '-') ):
                if isinstance( default[0], int ):
                    value.append( int( argv[ pos+offset ] ) )
                elif isinstance( default[0], float ):
                    value.append( float( argv[ pos+offset ] ) )
                else:
                    value.append( argv[ pos+offset ] )
                offset += 1
        elif isinstance( default, int ):
            value = int( argv[ pos + 1 ] )
        elif isinstance( default, float ):
            value = float( argv[ pos + 1 ] )
        else:
            value = argv[ pos + 1 ]
    else:
        if isinstance( default, list ) and isinstance( default[0], int ): value = []
    return value

