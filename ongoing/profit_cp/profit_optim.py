"""
use for Pyprofit fitting process, 
Taken from 'https://github.com/ICRAR/pyprofit/blob/master/examples/profit_optim.py'
and for our test, we set no mask images and segmentation by default.

PS: 

use total luminosity as correction for amplitude of model profiles, thus we need input
an magnitude ~ (not Brightness) as initial luminosity fitting

for a given image, we also need to input an initial angle ~ (i.e. along the major-axies)
for luminosity distribution fitting

Radius during fitting is in unit of pixel number!


PS:

PSF input ~ (give the magnitude of PSF and location on images)
modellist = list(
			pointsource = list(
				xcen = c(34,10,85),
				ycen = c(74,64,13),
				mag = c(10,13,16) )
			)

Sky input ~ (a constant background)

modellist = list(
			sky = list(
				bg = 3e-12 )
			)


PS:

combined model case ~ (just list all model for profile modeling)

modellist = list(
	sersic = list(
		xcen = c(50, 50),
		ycen = c(50, 50),
		mag = c(12, 13),
		re = c(14, 5),
		nser = c(1, 8),
		ang = c(46, 80),
		axrat = c(0.4, 0.6),
		box = c(0,-0.3)
	),

	ferrer = list(
		xcen = 50,
		ycen = 50,
		mag = 14,
		rout = 12,
		a = 0.3,
		b = 1.5,
		ang = 130,
		axrat = 0.2,
		box = 0.5
	),

	pointsource = list(
		xcen
		= c(34,10,85),
		ycen
		= c(74,64,13),
		mag = c(18,15,16)
	),

	sky = list(
		bg = 3e-12
	)
)

"""

import math
import pyprofit

from scipy import signal
from scipy import stats
from scipy import fftpack

import numpy as np
import time


class Data(object):
	pass


def to_pyprofit_image(params, data, use_mask = True ):

	# merge, un-sigma all, un-log some
	sigmas_tofit = data.sigmas[data.tofit]
	allparams = data.model0.copy()
	allparams[data.tofit] = params * sigmas_tofit
	allparams[data.tolog] = 10**allparams[data.tolog]

	fields = ['xcen','ycen','mag','re','nser','ang','axrat','box']
	s1params = [x for i,x in enumerate(allparams) if i%2 == 0]
	s2params = [x for i,x in enumerate(allparams) if i%2 != 0]

	if hasattr(data, 'psf') and len(data.psf) > 0:
		fields.append('convolve')
		s1params.append(True)
		s2params.append(True)

	sparams = [{name: val for name, val in zip(fields, params)} for params in (s1params, s2params)]
	if data.verbose:
		print(sparams)

	profit_model = {'width':  data.image.shape[1],
					'height': data.image.shape[0],
					'magzero': data.magzero,
					'psf': data.psf,
					'profiles': {'sersic': sparams}
					}
	if use_mask:
		profit_model['calcmask'] = data.calcregion
	image, _ = pyprofit.make_model( profit_model )

	return allparams, np.array(image)


def profit_like_model(params, data):

	# Get the priors sum
	priorsum = 0
	sigmas_tofit = data.sigmas[data.tofit]

	for i, p in enumerate(data.priors):
		priorsum += p(data.init[i] - params[i] * sigmas_tofit[i])

	# Calculate the new model
	allparams, modelim = to_pyprofit_image(params, data)

	# Scale and stuff
	scaledata = (data.image[data.region] - modelim[data.region]) / data.sigim[data.region]
	variance = scaledata.var()
	dof = 2 * variance / ( variance - 1 )
	dof = max( min( dof, float('inf') ) ,0 )

	ll = np.sum(stats.t.logpdf(scaledata, dof))
	lp = ll + priorsum
	lp = -lp

	if data.verbose:
		print( lp, { name: val for name, val in zip(data.names, allparams) } )
	return lp

def profit_setup_data(magzero, image, mask, sigim, segim, psf, names, model0, tofit, 
						tolog, sigmas, priors, lowers, uppers):

	im_w, im_h = image.shape
	psf_w, psf_h = psf.shape

	##. the 'mask' setting is for segimentation region of galaxy image during fitting
	region = np.zeros(image.shape, dtype = bool)

	region[(im_w - psf_w) // 2:(im_w + psf_w) // 2][(im_h - psf_h) // 2:(im_h + psf_h) // 2] = True
	region[(im_w - psf_w) // 2:(im_w + psf_w) // 2][(im_h - psf_h) // 2:(im_h + psf_h) // 2] = True

	segim_center_pix = segim[ int(math.ceil(im_w / 2.)) ][int( math.ceil(im_h / 2.) ) ]
	region[ segim == segim_center_pix ] = True

	##. original version
	psf[ psf<0 ] = 0.
	calcregion = signal.convolve2d( region.copy(), psf+1, mode = 'same')
	calcregion = calcregion > 0


	# ##. use FFT for convolution
	# psf[ psf<0 ] = 0.

	# sz = ( im_w - psf_w, im_h - psf_h )
	# cp_psf = np.pad( psf, ( ( (sz[0] + 1 ) // 2, sz[0] // 2), ( (sz[1] + 1 ) // 2, sz[1] // 2) ), 
	# 				'constant', constant_values = 0.,)

	# fft_psf = fftpack.ifftshift( cp_psf + 1)
	# fft_img = fftpack.fft2( region )
	# fft_psf_1 = fftpack.fft2( fft_psf )
	# calcregion = np.real( fftpack.ifft2( fft_img * fft_psf_1 ) )

	# calcregion = calcregion > 0

	#.
	data = Data()

	data.magzero = magzero
	data.names = names
	data.model0 = model0
	data.tolog = np.logical_and(tolog, tofit)
	data.tofit = tofit
	data.sigmas = sigmas
	data.image = image
	data.sigim = sigim
	data.psf = psf
	data.priors = priors[tofit]
	data.region = region
	data.calcregion = calcregion
	data.verbose = False

	# copy initial parameters
	# log some, /sigma all, filter
	data.init = data.model0.copy()
	data.init[tolog] = np.log10( data.model0[ tolog ] )
	data.init = data.init / sigmas
	data.init = data.init[tofit]

	# Boundaries are scaled by sigma values as well
	data.bounds = np.array( list( zip( lowers / sigmas, uppers / sigmas ) ) )[ tofit ]

	return data
