# this file use for query source catalogue
import h5py
import numpy as np
import astropy.io.fits as fits

import mechanize
import pandas as pd 
from io import StringIO

url = 'http://skyserver.sdss.org/dr12/en/tools/search/sql.aspx'
'''
load = '/home/xkchen/mywork/ICL/data/test_data/'
file = 'frame-r-ra36.455-dec-5.896-redshift0.233.fits'
# comparation 
'''
r_select = 0.11286
ra = 36.455
dec = -5.896

c_ra0 = str(ra - r_select)
c_dec0 = str(dec - r_select)
c_ra1 = str(ra + r_select)
c_dec1 = str(dec + r_select)

data_set = """
SELECT ALL
   p.ra,p.dec,p.u,p.g,p.r,p.i,p.z,p.type,
   s.class, s.z as redshift, s.class
FROM PhotoObj AS p
   JOIN SpecObj AS s ON s.class = 'STAR'
WHERE
   p.ra BETWEEN %s AND %s
   AND p.dec BETWEEN %s AND %s
   AND p.type = 6
""" % (c_ra0, c_ra1, c_dec0, c_dec1)

br = mechanize.Browser()
resp = br.open(url)
resp.info()
br.select_form(name = "sql")
br['cmd'] = data_set
br['format'] = ['csv']
response = br.submit()

doc = open('/home/xkchen/mywork/ICL/query/SDSS_SQL_data_out.txt', 'w') # all the data is in a line
print(response.get_data(), file = doc)
doc.close()

# change the data from one line to a table
d = open( '/home/xkchen/mywork/ICL/query/SDSS_SQL_data_out.txt' ).readlines()
d = d[0] 
d = d.split( '\\n' )
f = open( '/home/xkchen/mywork/ICL/query/SDSS_SQL_data_out.txt', 'w' )
for l in d:
    f.write( l + '\n' )
f.close()
