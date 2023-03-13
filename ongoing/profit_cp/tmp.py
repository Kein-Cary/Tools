import math
import pyprofit

from scipy import signal
from scipy import stats
from scipy import fftpack

import numpy as np


class Data(object):
	pass


def cp_to_pyprofit_image(params, data, model_str, use_mask = True):

	# merge, un-sigma all, un-log some
	sigmas_tofit = data.sigmas[data.tofit]
	allparams = data.model0.copy()
	allparams[data.tofit] = params * sigmas_tofit
	allparams[data.tolog] = 10**allparams[data.tolog]

	## model list and parameter selection

	# moffat profile
	if model_str == 'moffat':
		fields = ['xcen', 'ycen', 'mag', 'fwhm', 'con', 'ang', 'axrat', 'box']

	# Broken exponential profile
	if model_str == 'brokenexp':
		fields = ['xcen', 'ycen', 'mag', 'h1', 'h2', 'rb', 'ang', 'a', 'axrat', 'box']

	# Modified Ferrer's profile
	if model_str == 'ferrer':
		fields = ['xcen', 'ycen', 'mag', 'rout', 'a', 'b', 'ang', 'axrat', 'box']

	# Modified empirical King profile
	if model_str == 'king':
		fields = ['xcen', 'ycen', 'mag', 'rc', 'rt', 'a', 'ang', 'axrat', 'box']

	# sersic model
	if model_str == 'sersic':
		fields = ['xcen', 'ycen', 'mag', 're', 'nser', 'ang', 'axrat', 'box']

	# core-sersic profile
	if model_str == 'coresersic':
		fields = ['xcen', 'ycen', 'mag', 'rb', 're', 'nser', 'a', 'b', 'ang', 'axrat', 'box']

	##
	s1params = [ x for i,x in enumerate(allparams) if i%2 == 0 ]
	s2params = [ x for i,x in enumerate(allparams) if i%2 != 0 ]

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
					'profiles': { model_str: sparams}
					}

	if use_mask:
		profit_model['calcmask'] = data.calcregion
	image, _ = pyprofit.make_model(profit_model)

	return allparams, np.array(image)


def profit_like_model(params, data, model_str):

	# Get the priors sum
	priorsum = 0
	sigmas_tofit = data.sigmas[data.tofit]
	for i, p in enumerate(data.priors):
		priorsum += p(data.init[i] - params[i]*sigmas_tofit[i])

	# Calculate the new model
	allparams, modelim = cp_to_pyprofit_image( params, data, model_str )

	# Scale and stuff
	scaledata = (data.image[data.region] - modelim[data.region])/data.sigim[data.region]
	variance = scaledata.var()
	dof = 2*variance/(variance-1)
	dof = max(min(dof,float('inf')),0)

	ll = np.sum(stats.t.logpdf(scaledata, dof))
	lp = ll + priorsum
	lp = -lp

	if data.verbose:
		print(lp, {name: val for name, val in zip(data.names, allparams)})
	return lp


def profit_setup_data(magzero,
					  image, mask, sigim, segim, psf,
					  names, model0, tofit, tolog, sigmas, priors, lowers, uppers):

	im_w, im_h = image.shape
	psf_w, psf_h = psf.shape

	##.
	# All the center containing the PSF is considered, as well as the section of the 
	#                                                    image containing the galaxy
	# the 'mask' setting is for segimentation region of galaxy image during fitting

	region = np.zeros(image.shape, dtype = bool)
	region[(im_w - psf_w) // 2:(im_w + psf_w) // 2][(im_h - psf_h) // 2:(im_h + psf_h) // 2] = True
	segim_center_pix = segim[ int(math.ceil(im_w / 2.)) ][int( math.ceil(im_h / 2.) ) ]
	region[ segim == segim_center_pix ] = True

	##. PSF and galaxy image convolution ~ (for small array, i.e. 200*200)
	psf[ psf<0 ] = 0
	calcregion = signal.convolve2d( region.copy(), psf+1, mode = 'same')
	calcregion = calcregion > 0

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
	data.bounds = np.array( list( zip( lowers/sigmas,uppers/sigmas ) ) )[ tofit ]

	return data

