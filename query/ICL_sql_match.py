import scipy.stats as sts
import numpy as np
import astropy.io.fits as fits

import pandas as pd 
import astropy.wcs as awc
import astropy.units as U
from astropy import cosmology as apcy

import handy.scatter as hsc
import matplotlib as mpl
import matplotlib.pyplot as plt
## 
Test_model = apcy.Planck15.clone(H0 = 67.74, Om0 = 0.311)
H0 = Test_model.H0.value
h = H0/100
Omega_m = Test_model.Om0
Omega_lambda = 1.-Omega_m
Omega_k = 1.- (Omega_lambda + Omega_m)

pixel = 0.396
rad2asec = U.rad.to(U.arcsec)

load = '/home/xkchen/mywork/ICL/data/test_data/'
file = 'frame-r-ra36.455-dec-5.896-redshift0.233.fits'

data = fits.open(load + file)
img = data[0].data
wcs = awc.WCS(data[0].header)
x_side = data[0].data.shape[1]
y_side = data[0].data.shape[0]

cat = pd.read_csv('/home/xkchen/mywork/ICL/query/SDSS_SQL_data.txt', skiprows = 1)
ra = np.array(cat['ra'])
dec = np.array(cat['dec'])
mag = np.array(cat['r'])

'''
R0 = np.array(cat['psffwhm_r'])
iu = R0 >= 0
R = 4.25*R0[iu]/pixel
Ra = ra[iu]
Dec = dec[iu]
Mag = mag[iu]
'''
x, y = wcs.all_world2pix(Ra*U.deg, Dec*U.deg, 1)
ia = (x >= 0) & (x <= x_side)
ib = (y >= 0) & (y <= y_side)
ie = (mag <= 20)
ic = ia & ib & ie
comx = x[ic]
comy = y[ic]
comr = 2*1.5/pixel

R_ph = rad2asec/(Test_model.angular_diameter_distance(z = 0.233).value)
R_p = R_ph/pixel
cra = 36.455
cdec = -5.896
cenx, ceny = wcs.all_world2pix(cra*U.deg, cdec*U.deg, 1)

Numb = len(cr)
mask_B = np.ones((img.shape[0], img.shape[1]), dtype = np.float)
ox = np.linspace(0,2047,2048)
oy = np.linspace(0,1488,1489)
basic_coord = np.array(np.meshgrid(ox,oy))

# circle
for k in range(Numb):
	xc = comx[k]
	yc = comy[k]
	idr = np.sqrt((xc - basic_coord[0,:])**2 + (yc - basic_coord[1,:])**2)/cr[k]
	jx = idr <= 1
	jx = (-1)*jx+1
	mask_B = mask_B*jx

mirro_B = mask_B*img
plt.imshow(mirro_B, cmap = 'Greys', origin = 'lower', norm = mpl.colors.LogNorm())
hsc.circles(cenx, ceny, s = R_p, fc = '', ec = 'b', )
hsc.circles(comx, comy, s = cr, fc = '', ec = 'r', lw = 1)
plt.xlim(0, 2048)
plt.ylim(0, 1489)
plt.savefig('mask_test.png', dpi = 600)
plt.show()
raise
plt.hist(mag, histtype = 'step', color = 'b', label = 'all star')
plt.hist(mag[ic], histtype = 'step', color = 'r', label = 'mask sample')
plt.legend(loc = 2)
plt.xlabel('$Mag_r$')
plt.ylabel('Number')
plt.savefig('star-mag-distribution.png')
plt.close()
