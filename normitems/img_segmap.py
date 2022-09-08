import h5py
import numpy as np
import pandas as pds
import astropy.wcs as awc
import subprocess as subpro
import astropy.io.ascii as asc
import astropy.io.fits as fits

import astropy.units as U
import astropy.constants as C


###... for any input image and catalog mask
def simp_finder_func( d_file, out_file, config_file, params_file, tmp_file = None):

	"""
	d_file : path where image data saved (include file-name structure: '/xxx/xxx/xxx.fits')
	out_file : save sources information (detected by SExTractor)

	config_file : the '.sex' file of sextractor, setting for source detection
	params_file : the '.param' file of sextractor, setting for output parameters (i.e. position, Kron radius)

	tmp_file : ".fits" file, image use to find objects ( in case d_file is '.fits.bz2' or others non-fits files)
	"""

	param_A = config_file
	out_param = params_file

	img_data = fits.open( d_file )
	img = img_data[0].data
	head = img_data[0].header

	if tmp_file is not None:

		#. need to change the d_file in format of fits
		hdu = fits.PrimaryHDU()
		hdu.data = img
		hdu.header = head
		hdu.writeto( tmp_file, overwrite = True)
		file_source = tmp_file

	else:
		file_source = '~/pre_detect.fits'

	##.
	# cmd = 'sex '+ file_source + ' -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s' % (param_A, out_file, out_param)
	cmd = 'source-extractor '+ file_source + ' -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s' % (param_A, out_file, out_param)

	a = subpro.Popen(cmd, shell = True)
	a.wait()

	return

def array_mask_func( img_arr, cen_x, cen_y, cen_ar, cen_br, cen_chi, tag_arr = None):
	"""
	img_arr : input data array, on which masks will be applied on
	---------
	cen_x, cen_y, cen_ar, cen_br, cen_chi : mask catalog, record the location of objects for given array
	including position (cen_x, cen_y), major and minor axis (cen_ar, cen_br), position angle (cen_chi)								
	---------
	tag_arr : location of objects which will be remained after masking, including position, shape
	tag_x, tag_y, tag_a, tag_b, tag_phi = tag_arr[:]

	"""
	if tag_arr is not None:

		tag_x, tag_y, tag_a, tag_b, tag_phi = tag_arr[:]

		ef1 = ( tag_x - cen_x ) * np.cos( tag_phi ) + ( tag_y - cen_y ) * np.sin( tag_phi )
		ef2 = ( tag_y - cen_y ) * np.cos( tag_phi ) - ( tag_x - cen_x ) * np.sin( tag_phi )
		er = ef1**2 / tag_a**2 + ef2**2 / tag_b**2
		idx = er < 1

		if np.sum( idx ) >= 1:

			id_bcg = np.where( idx == True )[0]

		if np.sum( idx ) == 0:
			id_bcg = np.array( [] )

	else:
		id_bcg = np.array( [] )


	major = cen_ar / 2
	minor = cen_br / 2
	senior = np.sqrt( major**2 - minor**2 )

	Numb = len( major )

	mask_path = np.ones( (img_arr.shape[0], img_arr.shape[1]), dtype = np.float32)
	ox = np.linspace(0, img_arr.shape[1] - 1, img_arr.shape[1])
	oy = np.linspace(0, img_arr.shape[0] - 1, img_arr.shape[0])
	basic_coord = np.array( np.meshgrid(ox, oy) )


	# masking 'galaxies'
	for k in range( Numb ):

		xc = cen_x[k]
		yc = cen_y[k]

		lr = major[k]
		sr = minor[k]
		cr = senior[k]
		chi = cen_chi[k] * np.pi / 180

		if k in id_bcg:
			continue

		else:
			set_r = np.int( np.ceil(1.2 * lr) )
			la0 = np.max( [np.int(xc - set_r), 0])
			la1 = np.min( [np.int(xc + set_r + 1), img_arr.shape[1] ] )
			lb0 = np.max( [np.int(yc - set_r), 0] ) 
			lb1 = np.min( [np.int(yc + set_r + 1), img_arr.shape[0] ] )

			df1 = (basic_coord[0,:][lb0: lb1, la0: la1] - xc)* np.cos(chi) + (basic_coord[1,:][lb0: lb1, la0: la1] - yc)* np.sin(chi)
			df2 = (basic_coord[1,:][lb0: lb1, la0: la1] - yc)* np.cos(chi) - (basic_coord[0,:][lb0: lb1, la0: la1] - xc)* np.sin(chi)
			fr = df1**2 / lr**2 + df2**2 / sr**2
			jx = fr <= 1

			iu = np.where(jx == True)
			iv = np.ones((jx.shape[0], jx.shape[1]), dtype = np.float32)
			iv[iu] = np.nan
			mask_path[lb0: lb1, la0: la1] = mask_path[lb0: lb1, la0: la1] * iv

	mask_img = mask_path * img_arr

	return mask_img

