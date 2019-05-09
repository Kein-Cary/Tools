#!/usr/bin/env python2

import sqlcl
import os
import sys
import string

infile_dir = sys.argv[1]
oufile_dir = sys.argv[2]

infile_list = os.listdir( infile_dir )
url = sqlcl.default_url
fmt = 'csv'

for f in infile_list:
    infile = './' + os.path.join( infile_dir, f )
    oufile = './' + os.path.join( oufile_dir, f ) + '.' + fmt

    print "request %s ..."%infile
    qry = open( infile ).read()
    file = sqlcl.query( qry, url, fmt )
    line = file.readline()

    print "save %s to %s..."%( infile, oufile )
    oufile_fd = open( oufile, 'w' )
    sqlcl.write_header( oufile_fd, "#", url, qry )

    if line.startswith( 'ERROR' ):
        oufile_fd.write( 'ERROR' )
        oufile_fd.close()
        continue

    while line:
        oufile_fd.write( string.rstrip(line) + os.linesep )
        line = file.readline()
