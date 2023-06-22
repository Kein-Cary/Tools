import h5py
import numpy as np
import pandas as pds

#.
import scipy.stats as sts
from skimage import measure


##. basic source detection
def groups_find_func(img_data, threshold, pont_num = None):
	"""
	img_data : img will be use to find point groups

	threshold : the surface brightness limitation,
	general it is correlated with the 

	pont_num : point number limitation for groups
	[may also set the points number limit for the groups selection]
	"""
	identi = img_data > threshold
	copy_arr = img_data + 0.

	copy_arr[identi] = 1
	copy_arr[identi == False] = 0
	Ny, Nx = img_data.shape[0], img_data.shape[1]

	coord_x = []
	coord_y = []
	source_npix = []

	for kk in range(Ny):

		for tt in range(Nx):

			sub_n = 0
			sub_x = []
			sub_y = []

			if copy_arr[kk, tt] > 0:

				sub_x.append(tt)
				sub_y.append(kk)

				copy_arr[kk, tt] = 0
				tmp_len = 0
				sub_n += 1

				while sub_n != tmp_len:

					cont_n = tmp_len + 0
					tmp_len = sub_n + 0

					for nn in range(cont_n, sub_n):

						sor_x, sor_y = sub_x[nn], sub_y[nn]

						for pp in range(-1, 2):
							for qq in range(-1, 2):

								da0 = (sor_x + qq >= 0) & (sor_x + qq < Nx)
								da1 = (sor_y + pp >= 0) & (sor_y + pp < Ny)

								if (pp == 0) & (qq == 0):
									continue

								if (da0 & da1):
									if copy_arr[sor_y + pp, sor_x + qq] > 0:
										sub_x.append(sor_x + qq)
										sub_y.append(sor_y + pp)
										sub_n += 1
										copy_arr[sor_y + pp, sor_x + qq] = 0

				## record the source info.
				coord_x.append(sub_x)
				coord_y.append(sub_y)
				source_npix.append(sub_n)

			else:
				continue

	##.
	if pont_num is None:
		return source_npix, coord_x, coord_y

	else:
		id_nx = np.array( source_npix ) >= pont_num
		mp_dex = np.where( id_nx )[0]

		cc_Nobj = np.sum( id_nx )

		cc_coord_x = []
		cc_coord_y = []
		cc_n_pix = []

		for dd in range( cc_Nobj ):
			_dd_x = coord_x[ mp_dex[ dd ] ]
			_dd_y = coord_y[ mp_dex[ dd ] ]

			cc_coord_x.append( _dd_x )
			cc_coord_y.append( _dd_y )
			cc_n_pix.append( len( _dd_x ) )

		return cc_n_pix, cc_coord_x, cc_coord_y


def segmap_obj_func( img_arr, coord_x, coord_y, min_Npix = 1):
	"""
	img_arr : the image used to find objects
	coord_x, coord_y : location of sources pixels
	min_Npix : the minium pixel number ~ ( for an object)
	"""

	N_xx = np.array( [len(ll) for ll in coord_x] )

	id_mpx = N_xx >= min_Npix
	mp_dex = np.where( id_mpx )[0]

	##.
	N_obj = np.sum( id_mpx )
	cat_x = np.zeros( N_obj,)
	cat_y = np.zeros( N_obj,)
	cat_a = np.zeros( N_obj,)
	cat_b = np.zeros( N_obj,)
	cat_PA = np.zeros( N_obj,)

	for dd in range( N_obj ):

		_dd_x = np.array( coord_x[ mp_dex[ dd ] ] )
		_dd_y = np.array( coord_y[ mp_dex[ dd ] ] )

		_dd_flux = img_arr[ _dd_y, _dd_x ]

		cen_x = np.nansum( _dd_flux * _dd_x ) / np.nansum( _dd_flux )
		cen_y = np.nansum( _dd_flux * _dd_y ) / np.nansum( _dd_flux )

		X2 = np.nansum( _dd_flux * _dd_x**2 ) / np.nansum( _dd_flux ) - cen_x**2
		Y2 = np.nansum( _dd_flux * _dd_y**2 ) / np.nansum( _dd_flux ) - cen_y**2
		XY2 = np.nansum( _dd_flux * _dd_x * _dd_y ) / np.nansum( _dd_flux ) - cen_x * cen_y

		theta2 = 0.5 * np.arctan( 2 * XY2 / (X2**2 - Y2**2) )

		##. shape parameter
		mf0 = ( (X2**2 - Y2**2) / 2 )**2
		A22 = (X2**2 + Y2**2) / 2 + np.sqrt( mf0 + XY2**2 )
		B22 = (X2**2 + Y2**2) / 2 - np.sqrt( mf0 + XY2**2 )

		A2 = np.sqrt( A22 )
		B2 = np.sqrt( B22 )

		##. elongation and ellipticity
		elong = A2 / B2
		ellip = 1 - B2 / A2

		##.
		cat_x[dd] = cen_x
		cat_y[dd] = cen_y
		cat_a[dd] = A2
		cat_b[dd] = B2
		cat_PA[dd] = theta2

	return N_obj, cat_x, cat_y, cat_a, cat_b, cat_PA


def main():

	### === test for segmentation

	import time
	import matplotlib as mpl
	import matplotlib.pyplot as plt

	from astropy.modeling import models
	import matplotlib.patches as mpathes


	imshape = (300, 400)

	y, x = np.indices(imshape)

	# Generate random source model list
	rng = np.random.default_rng(0)
	nsrc = 100
	model_params = [
		dict(amplitude = rng.uniform(.5, 1),
			x_mean = rng.uniform(0, imshape[1] - 1),
			y_mean = rng.uniform(0, imshape[0] - 1),
			x_stddev = rng.uniform(2, 6),
			y_stddev = rng.uniform(2, 6),
			theta = rng.uniform(0, 2 * np.pi))
		for _ in range(nsrc)]

	model_list = [models.Gaussian2D(**kwargs) for kwargs in model_params]

	# Render models to image using bounding boxes
	bb_image = np.zeros(imshape)

	for model in model_list:
		model.render(bb_image)


	# Render models to image using full evaluation
	full_image = np.zeros(imshape)

	for model in model_list:
		model.bounding_box = None
		model.render(full_image)

	Nx, Ny = 400, 300
	BG_arr = np.random.random( (Ny, Nx),)

	input_img = full_image + BG_arr

	lim_x = np.median( input_img ) + 1.5 * np.std( input_img )

	source_n, coord_x, coord_y = groups_find_func( input_img, lim_x )

	loop_n = len(source_n)

	N_obj, cat_x, cat_y, cat_a, cat_b, cat_PA = segmap_obj_func( input_img, coord_x, coord_y, min_Npix = 5)


	fig = plt.figure( figsize = (11,11) )
	ax0 = fig.add_axes( [0.10, 0.53, 0.40, 0.42] )
	ax1 = fig.add_axes( [0.51, 0.53, 0.40, 0.42] )
	ax2 = fig.add_axes( [0.10, 0.11, 0.40, 0.42] )

	ax0.imshow( input_img, origin = 'lower', cmap = 'Greys', vmin = -3, vmax = 3,)

	ax1.imshow( input_img, origin = 'lower', cmap = 'Greys', vmin = -3, vmax = 3,)

	for ll in range( loop_n ):
		ax1.scatter(coord_x[ll], coord_y[ll], marker = 'o', s = 10, color = mpl.cm.rainbow(ll / loop_n),)

	ax2.imshow( input_img, origin = 'lower', cmap = 'Greys', vmin = -3, vmax = 3,)

	for ll in range( N_obj ):
		rect = mpathes.Ellipse( xy = (cat_x[ll], cat_y[ll]), width = cat_a[ll], height = cat_b[ll], 
								angle = cat_PA[ll] + 90, ec = 'r', fc = 'none', ls = '-',)
		ax2.add_patch( rect )

	plt.savefig('/home/xkchen/obj_limit_test.png', dpi = 300)
	plt.close()

	raise


if __name__ == "__main__":
	main()

