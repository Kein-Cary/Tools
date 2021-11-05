import h5py
import numpy as np
import astropy.io.fits as fits
from astropy.table import Table

import mechanize
import pandas as pds
from io import StringIO

import astropy.units as U
import astropy.constants as C
from astropy import cosmology as apcy
from astropy.coordinates import SkyCoord
from hscquery import HSC


band = ['r', 'g', 'i']

### === ### query HSC catalog
## divide sample by i-band Cmodel_mag
'''
# i_mag_arr = np.arange( 14, 24 ) # brighter than 24
i_mag_arr = np.arange( 24, 29, 0.5 ) # fainter than 24

N_subs = len( i_mag_arr )

for jj in range( 3, N_subs ):

    if jj == 0:
        cc_mag_0 = 0.
        cc_mag_1 = i_mag_arr[ jj ]

    else:
        cc_mag_0 = i_mag_arr[ jj - 1 ]
        cc_mag_1 = i_mag_arr[ jj ]

    sql_str = """
    SELECT
        object_id, ra, dec,
        a_g, a_r, a_i,
        g_cmodel_mag, r_cmodel_mag, i_cmodel_mag,
        specz_redshift

    FROM
        pdr3_wide.forced
        LEFT JOIN pdr3_wide.specz AS specz USING (object_id)
    WHERE
        isprimary
        AND g_cmodel_mag >0
            AND NOT g_cmodel_mag = 'NaN'
            AND NOT g_cmodel_mag = 'inf'
        AND r_cmodel_mag >0
            AND NOT r_cmodel_mag = 'NaN'
            AND NOT r_cmodel_mag = 'inf'
        AND i_cmodel_mag >0
            AND NOT i_cmodel_mag = 'NaN'
            AND NOT i_cmodel_mag = 'inf'

        AND i_extendedness_value = 1

        AND NOT i_extendedness_flag
        AND NOT g_cmodel_flag
        AND NOT r_cmodel_flag
        AND NOT i_cmodel_flag
        AND i_cmodel_mag >= %.1f
        AND i_cmodel_mag <= %.1f

    """ % ( cc_mag_0, cc_mag_1 )

    print( sql_str )

    output_format = 'csv'
    out_file = '/home/xkchen/figs/hsc_cat_sql/HSC_catalog/HSC-dr3_galaxy_gri-cat_i-cmag-up-%.1f.csv' % i_mag_arr[ jj ]
    delet_job = True

    h_sql = HSC( survey = 'wide' )
    data = h_sql.send_query( sql_str, output_file = out_file, delete_job = delet_job)

raise
'''

"""
## merger catalog into fits file
i_mag_0 = np.arange(14, 25)
i_mag_1 = np.arange(24.5, 29, 0.5)
i_mag = np.r_[ i_mag_0, i_mag_1 ]
N_sub = len( i_mag )

dat = pds.read_csv('/home/xkchen/data/HSC_data/HSC_catalog/' + 
                    'HSC-dr3_galaxy_gri-cat_i-cmag-up-14.0.csv', skiprows = 3)
keys = dat.columns
N_ks = len( keys )
pas = keys[0].replace('# ', '')
pt_lis = [ pas ]

for jj in range( 1, N_ks ):
    pt_lis.append( keys[ jj ] )

tmp_arr = []
for jj in range( N_ks ):
    tmp_arr.append( np.array([]) )

#. merger catalog
for jj in range( N_sub ):

    dat = pds.read_csv('/home/xkchen/data/HSC_data/HSC_catalog/' + 
                        'HSC-dr3_galaxy_gri-cat_i-cmag-up-%.1f.csv' % i_mag[ jj ], skiprows = 3)

    for nn in range( N_ks ):
        sub_arr = np.array( dat[ keys[nn] ] )
        tmp_arr[nn] = np.r_[ tmp_arr[nn], sub_arr ]

sub_tab = Table( tmp_arr, names = pt_lis )
sub_tab.write('/home/xkchen/data/HSC_data/HSC_gri-band_query-cat.fits', format = 'fits')
"""

### === ### match to HSC
"""
# hsc_dat = fits.open( '/home/xkchen/data/HSC_data/HSC_full_PDR3.fits' ) #.. HSC observation in all five bands
hsc_dat = fits.open( '/home/xkchen/data/HSC_data/HSC_gri-band_query-cat.fits' ) #.. HSC observation in gri bands

HSC_ID = hsc_dat[1].data['object_id']
HSC_ra, HSC_dec = np.array( hsc_dat[1].data['ra'] ), np.array( hsc_dat[1].data['dec'] )
HSC_r_cmag, HSC_i_cmag = np.array( hsc_dat[1].data['r_cmodel_mag'] ), np.array( hsc_dat[1].data['i_cmodel_mag'] )
af_r, af_g, af_i = np.array( hsc_dat[1].data['a_r'] ), np.array( hsc_dat[1].data['a_g'] ), np.array( hsc_dat[1].data['a_i'] )

hsc_coord = SkyCoord( HSC_ra * U.deg, HSC_dec * U.deg )

band = ['r', 'g', 'i']
load = '/home/xkchen/fig_tmp/'

#. match the BCG position at z_ref
for ll in range( 3 ):

    band_str = band[ ll ]

    #. BCG position at z_ref 
    ref_dat = pds.read_csv( load + 'pkoffset_cat/' + 
                            'low_BCG_star-Mass_%s-band_photo-z-match_pk-offset_BCG-pos_cat_z-ref.csv' % band_str )

    ref_ra_0, ref_dec_0, ref_z_0 = np.array( ref_dat['ra'] ), np.array( ref_dat['dec'] ), np.array( ref_dat['z'] )
    ref_bcgx_0, ref_bcgy_0 = np.array( ref_dat['bcg_x'] ), np.array( ref_dat['bcg_y'] )


    ref_dat = pds.read_csv( load + 'pkoffset_cat/' + 
                            'high_BCG_star-Mass_%s-band_photo-z-match_pk-offset_BCG-pos_cat_z-ref.csv' % band_str )

    ref_ra_1, ref_dec_1, ref_z_1 = np.array( ref_dat['ra'] ), np.array( ref_dat['dec'] ), np.array( ref_dat['z'] )
    ref_bcgx_1, ref_bcgy_1 = np.array( ref_dat['bcg_x'] ), np.array( ref_dat['bcg_y'] )  

    print('low_N = ', len(ref_ra_0) )
    print('high_N = ', len(ref_ra_1) )

    ref_ra = np.r_[ ref_ra_0, ref_ra_1 ]
    ref_dec = np.r_[ ref_dec_0, ref_dec_1 ]
    ref_z = np.r_[ ref_z_0, ref_z_1 ]
    ref_bcgx = np.r_[ ref_bcgx_0, ref_bcgx_1 ]
    ref_bcgy = np.r_[ ref_bcgy_0, ref_bcgy_1 ]

    ref_coord = SkyCoord( ra = ref_ra * U.deg, dec = ref_dec * U.deg,)

    idx, d2d, d3d = ref_coord.match_to_catalog_sky( hsc_coord )
    id_lim = d2d.value < 2.7e-4

    #. matched information in HSC
    mp_ra, mp_dec = HSC_ra[ idx[ id_lim ] ], HSC_dec[ idx[ id_lim ] ]
    mp_ID, mp_i_cmag = HSC_ID[ idx[ id_lim ] ], HSC_i_cmag[ idx[ id_lim ] ]
    mp_ar, mp_ag, mp_ai = af_r[ idx[ id_lim ] ], af_g[ idx[ id_lim ] ], af_i[ idx[ id_lim ] ]

    keys = [ 'ra', 'dec', 'obj_ID', 'i_cmag', 'a_r', 'a_g', 'a_i' ]
    values = [ mp_ra, mp_dec, mp_ID, mp_i_cmag, mp_ar, mp_ag, mp_ai ]
    fill = dict( zip( keys, values ) )
    data = pds.DataFrame( fill )
    data.to_csv( '/home/xkchen/%s-band_HSC-matched_cat.csv' % band_str )


    lim_ra, lim_dec, lim_z = ref_ra[ id_lim ], ref_dec[ id_lim ], ref_z[ id_lim ]

    keys = [ 'ra', 'dec', 'z' ]
    values = [ lim_ra, lim_dec, lim_z ]
    fill = dict( zip( keys, values ) )
    data = pds.DataFrame( fill )
    data.to_csv( '/home/xkchen/SDSS_%s-band_sql_cat.csv' % band_str )

    print('done ! ')

"""


### === ### query the light profile in HSC
def hsc_aper_flux_sql( ra_set, dec_set, set_IDs, out_file, band_info):

    Nz = len( ra_set )

    for jj in range( Nz ):

        ra_g, dec_g = ra_set[ jj ], dec_set[ jj ]
        ID_g = set_IDs[ jj ]

        sql_str = """
            SELECT
                object_id, ra, dec,
                a_g, a_r, a_i,
                g_cmodel_mag, r_cmodel_mag, i_cmodel_mag,
                g_apertureflux_10_flux, r_apertureflux_10_flux, i_apertureflux_10_flux,
                g_apertureflux_10_fluxerr, r_apertureflux_10_fluxerr, i_apertureflux_10_fluxerr,

                g_apertureflux_15_flux, r_apertureflux_15_flux, i_apertureflux_15_flux,
                g_apertureflux_15_fluxerr, r_apertureflux_15_fluxerr, i_apertureflux_15_fluxerr,

                g_apertureflux_20_flux, r_apertureflux_20_flux, i_apertureflux_20_flux,
                g_apertureflux_20_fluxerr, r_apertureflux_20_fluxerr, i_apertureflux_20_fluxerr,

                g_apertureflux_30_flux, r_apertureflux_30_flux, i_apertureflux_30_flux,
                g_apertureflux_30_fluxerr, r_apertureflux_30_fluxerr, i_apertureflux_30_fluxerr,

                g_apertureflux_40_flux, r_apertureflux_40_flux, i_apertureflux_40_flux,
                g_apertureflux_40_fluxerr, r_apertureflux_40_fluxerr, i_apertureflux_40_fluxerr,

                g_apertureflux_57_flux, r_apertureflux_57_flux, i_apertureflux_57_flux,
                g_apertureflux_57_fluxerr, r_apertureflux_57_fluxerr, i_apertureflux_57_fluxerr,

                g_apertureflux_84_flux, r_apertureflux_84_flux, i_apertureflux_84_flux,
                g_apertureflux_84_fluxerr, r_apertureflux_84_fluxerr, i_apertureflux_84_fluxerr,

                g_apertureflux_118_flux, r_apertureflux_118_flux, i_apertureflux_118_flux,
                g_apertureflux_118_fluxerr, r_apertureflux_118_fluxerr, i_apertureflux_118_fluxerr,

                g_apertureflux_168_flux, r_apertureflux_168_flux, i_apertureflux_168_flux,
                g_apertureflux_168_fluxerr, r_apertureflux_168_fluxerr, i_apertureflux_168_fluxerr,

                g_apertureflux_235_flux, r_apertureflux_235_flux, i_apertureflux_235_flux,
                g_apertureflux_235_fluxerr, r_apertureflux_235_fluxerr, i_apertureflux_235_fluxerr

            FROM
                pdr3_wide.forced
                LEFT JOIN pdr3_wide.forced3 USING (object_id)
            WHERE
                isprimary
                AND object_id = %d
        """ % ID_g

        output_format = 'csv'
        delet_job = True

        h_sql = HSC( survey = 'wide' )
        data = h_sql.send_query( sql_str, output_file = out_file % (ra_g, dec_g, band_info), delete_job = delet_job)

    return


#. build the my_catalog used for HSC flux array matching
for kk in range( 3 ):

    dat = pds.read_csv('/home/xkchen/figs/hsc_cat_sql/%s-band_HSC-matched_cat.csv' % band[ kk ] )
    sub_ra, sub_dec = np.array( dat['ra'] ), np.array( dat['dec'] )
    sub_IDs = np.array( dat['obj_ID'] )

    out_files = '/home/xkchen/figs/hsc_cat_sql/HSC_aperF_cat/hsc_aper-flux_ra%.3f_dec%.3f_%s-band.csv'
    band_str = band[ kk ]

    Ns = len( sub_ra )

    #... query for BCGs
    # hsc_aper_flux_sql( sub_ra, sub_dec, sub_IDs, out_files, band_str )

    #... hsc query str list
    # doc = open( '/home/xkchen/%s-band_BCG-cat.txt' % band[ kk ], 'w')

    # for jj in range( Ns ):
    #     _str_line = '(%d, %.8f, %.8f),' % ( sub_IDs[ jj ], sub_ra[ jj ], sub_dec[ jj ] )
    #     print( _str_line, file = doc )

    # doc.close()

    doc = open( '/home/xkchen/%s-band_read_test.txt' % band_str, 'w')
    print( 'band = ', band_str, file = doc )
    for jj in range( Ns ):

        ra_g, dec_g = sub_ra[jj], sub_dec[jj]

        try:
            cat = pds.read_csv( out_files % (ra_g, dec_g, band_str ), skiprows = 3)
            _jj_mag = np.array( cat['i_cmodel_mag'] )[0]

        except:
            print('err_dex = ', jj, file = doc)
            print('ra = , dec = ', ra_g, dec_g, file = doc)

    doc.close()
