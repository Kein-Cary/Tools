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
        p.psffwhm_u, p.psffwhm_g, p.psffwhm_r, p.psffwhm_i, p.psffwhm_z,
        p.petroR90_u, p.petroR90_g, p.petroR90_r, p.petroR90_i, p.petroR90_z,

        p.deVRad_u, p.deVRad_g, p.deVRad_r, p.deVRad_i, p.deVRad_z,
        p.deVAB_u, p.deVAB_g, p.deVAB_r, p.deVAB_i, p.deVAB_z,
        p.deVPhi_u, p.deVPhi_g, p.deVPhi_r, p.deVPhi_i, p.deVPhi_z,

        p.expRad_u, p.expRad_g, p.expRad_r, p.expRad_i, p.expRad_z,
        p.expAB_u, p.expAB_g, p.expAB_r, p.expAB_i, p.expAB_z,
        p.expPhi_u, p.expPhi_g, p.expPhi_r, p.expPhi_i, p.expPhi_z
FROM PhotoObj AS p
WHERE
   p.ra BETWEEN %s AND %s
   AND p.dec BETWEEN %s AND %s
   AND p.type = 6
ORDER by p.r
""" % (c_ra0, c_ra1, c_dec0, c_dec1)

br = mechanize.Browser()
resp = br.open(url)
resp.info()

br.select_form(name = "sql")
br['cmd'] = data_set
br['format'] = ['csv']
response = br.submit()
s = str(response.get_data(), encoding = 'utf-8')
doc = open('/home/xkchen/mywork/ICL/query/SDSS_SQL_data.txt', 'w') # all the data is in a line
print(s, file = doc)
doc.close()
'''
# change the data from one line to a table
d = open( '/home/xkchen/mywork/ICL/query/SDSS_SQL_data.txt' ).readlines()
d = d[0] 
d = d.split( '\\n' )
f = open( '/home/xkchen/mywork/ICL/query/SDSS_SQL_data.txt', 'w' )
for l in d:
    f.write( l + '\n' )
f.close()
'''
