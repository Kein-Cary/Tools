import h5py
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

import astropy.io.fits as fits
import astropy.wcs as awc
import astropy.units as U
import astropy.constants as C

load = '/home/xkchen/mywork/ICL/data/test_data/'
file = 'frame-u-ra203.834-dec41.001-redshift0.228.fits'
sql = '/home/xkchen/tmp/ra203.834_dec41.001_z0.228.csv'

data = fits.open(load + file)
img = data[0].data
wcs = awc.WCS(data[0].header)
x_side = data[0].data.shape[1]
y_side = data[0].data.shape[0]

cata_log = pd.read_csv(sql, skiprows = 15)
ra = np.array(cata_log['ra'][0:-1])
dec = np.array(cata_log['dec'][0:-1])
x, y = wcs.all_world2pix(ra*U.deg, dec*U.deg, 1)

ia = (x >= 1e-3) & (x <= 2047)
ib = (y >= 1e-3) & (y <= 1488)
ic = ia & ib
com_x = x[ic]
com_y = y[ic]

plt.imshow(img, cmap = 'Greys', vmin = 1e-3, origin = 'lower', norm = mpl.colors.LogNorm())
plt.scatter(com_x, com_y, s = 5, marker = 'o', c = 'r', alpha = 0.35)
plt.savefig('/home/xkchen/tmp/match_test.png', dpi = 600)
plt.show()
