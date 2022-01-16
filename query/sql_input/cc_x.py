#!/usr/bin/env python2
import h5py
import numpy as np
import astropy.io.fits as fits
import sqlcl

import pandas as pds
from io import StringIO
from sql_ask import sdss_sql

url = "http://skyserver.sdss.org/dr12/en/tools/search/x_sql.aspx"
#url = 'http://cas.sdss.org/dr7/en/tools/search/sql.asp'

fmt = 'csv'

r_select = 0.1128 #(1025*0.396/3600)

def put_sql( z_set, ra_set, dec_set, ):

    Nz = len( z_set )

    for q in range(Nz):

        ra_g = ra_set[q]
        dec_g = dec_set[q]

        cord_z = z_set[q]
        cord_ra = ra_set[q]
        cord_dec = dec_set[q]  

        c_ra0 = str(ra_g - r_select)
        c_dec0 = str(dec_g - r_select)
        c_ra1 = str(ra_g + r_select)
        c_dec1 = str(dec_g + r_select)

        qry = """
        SELECT ALL
        p.ra,p.dec,p.u,p.g,p.r,p.i,p.z,p.type,  
        p.probPSF,
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

        out_file = './ra%.3f_dec%.3f_z%.3f.csv' % (cord_ra, cord_dec, cord_z)

        sdss_sql(url, qry, fmt, out_file, id_print = False,)

        if q == 5:
            break

    return

with h5py.File('./sdss_sql_img_catalog.h5', 'r') as f:
    catalogue = np.array(f['a'])
z = catalogue[0]
Ra = catalogue[1]
Dec = catalogue[2]

'''
with h5py.File('./sample_catalog.h5', 'r') as f:
    head_record = np.array(f['a'])
z_cod = head_record[0]
ra_cod = head_record[1]
dec_cod = head_record[2]
'''

put_sql( z, Ra, Dec, )

