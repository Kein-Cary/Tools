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

def main():

	####### test part
	import time
	import matplotlib as mpl
	import matplotlib.pyplot as plt

	Nx, Ny = 15, 10
	A = np.ones((Ny, Nx), dtype = np.float)
	for nn in range(Ny):
		for mm in range(Nx):
			pr = np.random.random()
			if pr > 0.5:
				pass
			else:
				A[nn, mm] = 5 * pr

	lim_x = 1.25
	source_n, coord_x, coord_y = groups_find_func(A, lim_x,)
	loop_n = len(source_n)

	plt.figure()
	plt.imshow(A, origin = 'lower', cmap = 'rainbow', vmin = -3, vmax = 3,)
	for ll in range( loop_n ):
		plt.scatter(coord_x[ll], coord_y[ll], marker = 'o', s = 10, color = mpl.cm.rainbow(ll / loop_n),)
	plt.show()

	raise

if __name__ == "__main__":
	main()
