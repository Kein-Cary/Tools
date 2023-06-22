"""
test for rotation of 3D points
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import matplotlib.patches as mpathes

import h5py
import numpy as np
import pandas as pds


### === func.s
def point_rotate_2D_func( x_arr, y_arr, x0, y0, angle ):
	"""
	x_arr, y_arr : points need to rotate
	angle : angle of rotation, in unit of 'rad'
	x0, y0 : the center of rotation
	"""

	R = [ 	[ np.cos( angle ), -np.sin( angle ), x0 * (1 - np.cos( angle ) ) + y0 * np.sin( angle ) ], 
			[ np.sin( angle ), np.cos( angle ),  y0 * (1 - np.cos( angle ) ) - x0 * np.sin( angle ) ], 
			[ 0,               0,                1 ] ]

	R = np.array( R )

	if type( x_arr ) == np.ndarray:
		P = np.array( [ x_arr, y_arr, [1] * x_arr.size ] )

	else:
		P = np.array( [ x_arr, y_arr ] )
		P = np.r_[ P, [1] * x_arr.size ]

	##. x, y, 1
	out_P = np.dot( R, P )

	return out_P


def point_rotate_3D_func( x, y, z, x0, y0, z0, alpha, beta, gamma ):
	"""
	x, y, z : points need to rotate
	x0, y0, z0 : the center of rotation
	alpha, beta, gamma : angle of rotation, in unit of 'rad'
	"""
	Rx = [ 	[ 1,  0,  0,  0 ], 
			[ 0, np.cos(alpha), -np.sin(alpha), y0 * (1 - np.cos(alpha) ) + z0 * np.sin(alpha) ], 
			[ 0, np.sin(alpha), np.cos(alpha), z0 * (1 - np.cos(alpha) ) - y0 * np.sin(alpha) ], 
			[ 0,  0,  0,  1 ] ]

	Ry = [ 	[ np.cos(beta),  0,  np.sin(beta),  x0 * (1 - np.cos(beta) ) - z0 * np.sin(beta) ], 
			[ 0,  1,  0,  0 ], 
			[ -np.sin(beta),  0,  np.cos(beta),  z0 * (1 - np.cos(beta) ) + x0 * np.sin(beta) ], 
			[ 0,  0,  0,  1 ] ]

	Rz = [ 	[ np.cos(gamma), -np.sin(gamma),  0,  x0 * (1 - np.cos(gamma) ) + y0 * np.sin(gamma) ], 
			[ np.sin(gamma), np.cos(gamma),  0,  y0 * (1 - np.cos(gamma) ) - x0 * np.sin(gamma) ], 
			[ 0,  0,  1,  0], 
			[ 0,  0,  0,  1] ]

	Rx = np.array( Rx )
	Ry = np.array( Ry )
	Rz = np.array( Rz )

	if type( x ) == np.ndarray:
		P = np.array( [ x, y, z, [1] * x.size ] )

	else:
		P = np.array( [ x, y, z ] )
		P = np.r_[ P, [1] * x.size ]

	##. x, y, z, 1
	out_P = np.dot( np.dot( np.dot( Rx, Ry ), Rz), P )

	return  out_P


### === figs
py = np.array( [ 1., 1., 1., 1. ] )
px = np.array( [ 1., 2., 3., 4. ] )
pz = np.array( [ 1., 2., 3., 4. ] )

angl_0 = 60 * np.pi / 180


##.
P_arr = np.array([ px, py ]).T
P0_arr = np.array([ 1, 1 ]).T


# Pr_arr = point_rotate_2D_func( px, py, P0_arr[0], P0_arr[1], angl_0 )
nx, ny, _dd_ = point_rotate_2D_func( px, py, P0_arr[0], P0_arr[1], angl_0 )

fig = plt.figure()
ax = fig.add_axes([0.12, 0.11, 0.80, 0.80])

ax.plot( px, py, 'ro', markersize = 15, alpha = 0.5,)
ax.plot( nx, ny, 'b*', alpha = 0.5, markersize = 15,)

ax.set_xlim( -1, 5 )
ax.set_ylim( -1, 5 )

plt.savefig('/home/xkchen/rotation_2D_test.png', dpi = 300)
plt.close()

raise

##... 
alpha, beta, gamma = 45 * np.pi / 180, 45 * np.pi / 180, 45 * np.pi / 180

nx, ny, nz, _dd_ = point_rotate_3D_func( px[0], py[0], pz[0], 0, 1, 0, alpha, beta, gamma )

fig = plt.figure()
ax = fig.add_axes([0.12, 0.11, 0.80, 0.80], projection = '3d')

ax.plot( px, py, pz, 'ro', alpha = 0.5, markersize = 20, )
ax.plot( nx, ny, nz, 'g*', alpha = 0.5, markersize = 20, )

ax.set_xlim( -1, 5 )
ax.set_ylim( -1, 5 )
ax.set_zlim( -1, 5 )

ax.view_init( 5, -3 )

plt.savefig('/home/xkchen/rotation_3D_test.png', dpi = 300)
plt.show()

