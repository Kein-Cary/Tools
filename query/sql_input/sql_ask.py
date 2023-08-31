#!/usr/bin/env python2
import h5py
import numpy as np
import astropy.io.fits as fits
import sqlcl

import pandas as pd 
from io import StringIO

#url = "http://skyserver.sdss.org/dr12/en/tools/search/x_sql.aspx"
#url = 'http://cas.sdss.org/dr7/en/tools/search/sql.asp'

def sdss_sql(url_link, query_str, query_fmt, out_file, id_print = False ):

    if id_print == True:
        print( query_str )

    file = sqlcl.query( query_str, url_link, query_fmt )

    fd = open( out_file, 'w' )
    #sqlcl.write_header( fd, "#", url, qry ) #also take the header and query infomation

    lines = file.readlines()

    for l in lines:
        dtl = l.decode('utf-8')
        fd.write( dtl )
    fd.close()

    return

