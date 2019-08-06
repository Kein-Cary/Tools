#!/usr/bin/env python2
import h5py
import numpy as np
import astropy.io.fits as fits
import sqlcl

#import mechanize
import pandas as pd 
#from io import StringIO
with h5py.File('./sdss_sql_img_catalog.h5') as f:
    catalogue = np.array(f['a'])
z = catalogue[0]
Ra = catalogue[1]
Dec = catalogue[2]

with h5py.File('./sample_catalog.h5') as f:
    head_record = np.array(f['a'])
z_cod = head_record[0]
ra_cod = head_record[1]
dec_cod = head_record[2]

#url = sqlcl.default_url
url = "http://skyserver.sdss.org/dr12/en/tools/search/x_sql.aspx"
fmt = 'csv'

r_select = 0.1128 #(1025*0.396/3600)
Nz = len(z)
def sdss_sql():

    #i = -1 
    for q in range(Nz):
        '''
        i = i+1
        if ( i % 10 != 0 ):
            continue
        '''
        time = z[q]
        ra = Ra[q]
        dec = Dec[q]

        c_ra0 = str(ra - r_select)
        c_dec0 = str(dec - r_select)
        c_ra1 = str(ra + r_select)
        c_dec1 = str(dec + r_select)

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

        cord_z = z_cod[q]
        cord_ra = ra_cod[q]
        cord_dec = dec_cod[q]
        print( qry )
        file = sqlcl.query( qry, url, fmt ) 

        #fd = open( '/home/xkchen/mywork/ICL/data/star_catalog/ra%.3f_dec%.3f_z%.3f.csv'%(cord_ra, cord_dec, cord_z), 'w' )
        fd = open( './ra%.3f_dec%.3f_z%.3f.csv'%(cord_ra, cord_dec, cord_z), 'w' )
        sqlcl.write_header( fd, "#", url, qry )
        
        lines = file.readlines()
        for l in lines:
            fd.write( l )
        fd.close()

    return

def main():
    sdss_sql()

if __name__ == "__main__":
    main()

