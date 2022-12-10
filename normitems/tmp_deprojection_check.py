"""
try to model the Surface Brightness profiles of satellites
"""
import sys 
sys.path.append('/home/xkchen/tool/Conda/Tools/normitems')

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

import h5py
import numpy as np
import pandas as pds
import astropy.wcs as awc
import astropy.io.ascii as asc
import astropy.io.fits as fits
import astropy.units as U
import astropy.constants as C

from astropy.table import Table, QTable
from astropy import cosmology as apcy
from scipy import interpolate as interp
from scipy import integrate as integ
from astropy.coordinates import SkyCoord
from pynverse import inversefunc
from scipy import optimize
import scipy.signal as signal

#.
from colossus.cosmology import cosmology as co_cosmos
from colossus.halo import profile_nfw
from colossus.halo import profile_einasto
from colossus.halo import profile_hernquist
from colossus.halo import profile_spline

#.
from Gauss_Legendre_factor import GaussLegendreQuadrature
from Gauss_Legendre_factor import GaussLegendreQuad_arr
from rho_projected import input_cosm_model
from rho_projected import cosmos_param
from rho_projected import rhom_set



### === ### cosmology
rad2asec = U.rad.to(U.arcsec)
Test_model = apcy.Planck15.clone(H0 = 67.74, Om0 = 0.311)
H0 = Test_model.H0.value
h = H0 / 100

Omega_m = Test_model.Om0
Omega_lambda = 1.-Omega_m
Omega_k = 1.- (Omega_lambda + Omega_m)
Omega_b = Test_model.Ob0

##.
params = {'flat': True, 'H0': 67.74, 'Om0': 0.311, 'Ob0': 0.049, 'sigma8': 0.81, 'ns': 0.95}
co_cosmos.addCosmology('myCosmo', params = params )
my_cosmo = co_cosmos.setCosmology( 'myCosmo' )


##.
input_cosm_model( get_model = Test_model )
cosmos_param()


### === ### func.s
def interp_F_func(x_new, x, y):

	tck = interp.splrep( x, y, s = 0 )
	y_new = interp.splev( x_new, tck, der = 0 )

	return y_new

##. deprojection
def I_func( Rx, rx, x, y ):

	tmp_F = interp_F_func( Rx, x, y)
	mf_1 = np.sqrt( Rx**2 - rx**2 )

	return -1 * tmp_F / ( mf_1 * np.pi )

def invers_trans_func( x, y ):
	"""
	x, y : projected array of x, and functions
	"""

	dy = np.gradient( y )
	d_lgx = np.gradient( np.log10( x ) )

	dy_dlgx = dy / d_lgx
	dy_dx = dy_dlgx / ( x * np.log( 10 ) )

	#. at the same distance but for 3D case 
	NN_r = len( x )

	y3d = np.zeros( NN_r, )

	for mm in range( NN_r ):

		# mm_I = integ.quad( I_func, x[ mm ], np.inf, args = ( x[mm], x, dy_dx ),)
		mm_I = integ.quad( I_func, x[ mm ], np.max( x ), args = ( x[mm], x, dy_dx ), )

		y3d[ mm ] = mm_I[ 0 ] + 0.

	return y3d

##. projection
def rho_func( Rx, rx, x, y ):

	tmp_F = interp_F_func( Rx, x, y)
	mf_1 = np.sqrt( Rx**2 - rx**2 )

	return 2 * Rx * tmp_F / mf_1

def trans_func( x, y ):
	"""
	x, y : array need to project
	"""

	NN_r = len( x )

	y2d = np.zeros( NN_r, )

	for mm in range( NN_r ):

		# mm_I = integ.quad( rho_func, x[ mm ], np.inf, args = ( x[mm], x, y ),)
		mm_I = integ.quad( rho_func, x[ mm ], np.max( x ), args = ( x[mm], x, y ), )

		y2d[ mm ] = mm_I[ 0 ] + 0.

	return y2d


### === ### deprojection test
Mh0 = 14.24  # M_sun / h
C_h0 = 5
z_h0 = 0.2

Vm = 200

rho_c, rho_m = rhom_set( z_h0 )

R200m = ( 3 * 10**Mh0 / (4 * np.pi * Vm * rho_m) )**(1/3)


# halo_dens = profile_nfw.NFWProfile( M = 10**Mh0, c = C_h0, z = z_h0, mdef = '200m')
halo_dens = profile_einasto.EinastoProfile( M = 10**Mh0, c = C_h0, z = z_h0, mdef = '200m')


# Nr = 1000
Nr = 500
rx = np.logspace( -2, np.log10( 10 * R200m ), Nr )

pros_3d = halo_dens.density( rx )
pros_2d = halo_dens.surfaceDensity( rx )


#. deprojection and projection
rho_Vm = Vm * rho_m

deproj_3d_pros = invers_trans_func( rx, pros_2d / rho_Vm )
deproj_3d_pros = deproj_3d_pros * rho_Vm

proj_2d_pros = trans_func( rx, pros_3d / rho_Vm )
proj_2d_pros = proj_2d_pros * rho_Vm


##. interpolation with ending points test



##.
fig = plt.figure()
ax = fig.add_axes([0.12, 0.32, 0.80, 0.63])
sub_ax = fig.add_axes([0.12, 0.11, 0.80, 0.21]) 

ax.plot( rx, pros_3d, 'r-', alpha = 0.5, label = '3D')
ax.plot( rx, pros_2d, 'b-', alpha = 0.5, label = '2D')

ax.plot( rx, deproj_3d_pros, 'g--', alpha = 0.5, lw = 2.5, label = 'Me Deprojection')
ax.plot( rx, proj_2d_pros, 'k--', alpha = 0.5, lw = 2.5, label = 'Me Projection')

sub_ax.plot( rx, deproj_3d_pros / pros_3d, ls = '--', lw = 2.5, color = 'g', alpha = 0.5,)

sub_ax.plot( rx, proj_2d_pros / pros_2d, ls = '--', color = 'k', alpha = 0.5,)


# ax.set_xlim( 1e-2, 5e3 )
# ax.set_xlabel('R [kpc / h] or r [kpc / h]')
ax.set_xscale('log')

ax.set_yscale('log')
ax.set_ylim( 1e2, 1e10 )
ax.set_ylabel('$\\Sigma(R)$ or $\\rho(r)$', fontsize = 12,)

ax.legend( loc = 3, frameon = False, fontsize = 12,)
ax.axvline( x = R200m, ls = ':', color = 'k',)

sub_ax.set_xlim( ax.get_xlim() )
sub_ax.set_xlabel('R [kpc / h] or r [kpc / h]', fontsize = 12,)
sub_ax.set_xscale('log')

sub_ax.set_ylim( 0.94, 1.06 )
sub_ax.set_ylabel( 'Ratio', fontsize = 12,)
sub_ax.yaxis.set_minor_locator( ticker.AutoMinorLocator() )

sub_ax.axvline( x = R200m, ls = ':', color = 'k',)
sub_ax.axhline( y = 1, ls = ':', color = 'c', alpha = 0.25,)

ax.tick_params( axis = 'both', which = 'both', direction = 'in', labelsize = 12,)
sub_ax.tick_params( axis = 'both', which = 'both', direction = 'in', labelsize = 12,)
ax.set_xticklabels( labels = [] )

plt.savefig('/home/xkchen/deprojection_test.png', dpi = 300)
plt.close()

