import matplotlib as mpl
import matplotlib.pyplot as plt

import h5py
import numpy as np
import pandas as pds

def groups_find_func(img_data, threshold, pont_num = None):
	"""
	img_data : img will be use to find point groups
	threshold : the surface brightness limitation
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
	source_n = []

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
				source_n.append(sub_n)

			else:
				continue
	return source_n, coord_x, coord_y

with h5py.File('test_over-sb.h5', 'r') as f:
	over_sb = np.array(f['a'])

lim_sb = 5.

source_n, coord_x, coord_y = groups_find_func(over_sb, lim_sb)

plt.figure()
ax = plt.subplot(111)
ax.set_title('groups based on over_SB img [$\\Delta $ > %.1f]' % lim_sb)
tf = ax.imshow(over_sb, origin = 'lower', cmap = 'seismic', vmin = -5, vmax = 5,)
plt.colorbar(tf, ax = ax, fraction = 0.035, pad = 0.01, label = '$\\Delta = (\\mu_{P} - \\mu_{C})$ / $\\sigma_{C}$')

for mm in range( len(source_n) ):
	tmp_x = np.array(coord_x[mm])
	tmp_y = np.array(coord_y[mm])
	ax.scatter(tmp_x, tmp_y, s = 10, color = mpl.cm.hsv(mm / len(source_n) ), marker = 'o', label = 'group %d' % mm)
ax.legend(loc = 'left center',)

plt.savefig('groups_test_%.1f-sigma.png' % lim_sb, dpi = 300)
plt.close()
