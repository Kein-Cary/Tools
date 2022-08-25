import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

import numpy as np
from scipy import interpolate as interp
from scipy import integrate as integ
#.
from Gauss_Legendre_factor import GaussLegendreQuad_arr


##... stellar mass distribution
def sersic_func( R, n, Ie, Re ):
    """
    for n > 0.36
    """
    bc0 = 4 / (405 * n)
    bc1 = 45 / (25515 * n**2)
    bc2 = 131 / (1148175 * n**3)
    bc3 = 2194697 / ( 30690717750 * n**4 )

    b_n = 2 * n - 1/3 + bc0 + bc1 + bc2 - bc3
    mf0 = -b_n * ( R / Re )**( 1 / n )

    Ir = Ie * np.exp( mf0 + b_n )

    return Ir

def dsf_dR_func( R, n, Ie, Re ):

    bc0 = 4 / (405 * n)
    bc1 = 45 / (25515 * n**2)
    bc2 = 131 / (1148175 * n**3)
    bc3 = 2194697 / ( 30690717750 * n**4 )

    b_n = 2 * n - 1/3 + bc0 + bc1 + bc2 - bc3

    mf0 = sersic_func( R, n, Ie, Re )
    mf1 = ( -b_n / ( n * Re ) ) * ( R / Re )**( 1/n - 1 )

    return mf0 * mf1

def integ_func(R, x, n, Ie, Re ):
    """
    PS: R >= x
    """
    mf0 = dsf_dR_func( R, n, Ie, Re )
    mf1 = np.sqrt( R**2 - x**2 )

    return (-1 / np.pi) * mf0 / mf1

def sersic_3D_func( x, n, Ie, Re ):

    mf0 = integ.quad( integ_func, x, np.inf, args = (x, n, Ie, Re),)

    return mf0


##. Mh or Mstar infer func
def dM_2d_func( R, n, Ie, Re ):
    return R * sersic_func( R, n, Ie, Re )

def infer_Ie_func( Ie, n, Re ):

    mf0 = integ.quad( dM_2d_func, 0, Re, args = (n, Ie, Re),)[0]

    return mf0 * 2 * np.pi

def R_func( R, a, b, c):
    return a + (R + b)**c


### === sersic deprojection~(S. H. Price et.al 2022)
if __name__ == "__main__":

    lg_Mstar = 11.5  ##. M_sun 

    ##. assume inclination angle = 90, 
    ##. c / a = 0.5 is the the intrinsic axis ratio of projected density profile
    ##. n = 4, or n is small values, since those galaxy are satellites, even some of them are star forming

    alon_A = 90
    q_3d = 0.5
    q_obs = 0.5

    n = 4
    Re = 15  ##. kpc

    ##. infer Ie
    Ie_0 = 10

    half_Mass = infer_Ie_func( Ie_0, n, Re )
    eta_M = 0.5 * ( 10**lg_Mstar / half_Mass )

    Ie = Ie_0 * eta_M


    ##. 2D and 3D density profile
    N_ss = 200
    R = np.logspace( -1, 2.5, N_ss )
    rx = np.logspace( -1, 2.5, N_ss )

    rho_2D = sersic_func( R, n, Ie, Re )

    rho_3D = np.zeros( N_ss,)

    for mm in range( N_ss ):

        rho_mm = sersic_3D_func( rx[ mm ], n, Ie, Re )
        rho_3D[ mm ] = rho_mm[0]


    fig = plt.figure()
    ax = fig.add_axes([0.12, 0.11, 0.80, 0.85])

    ax.plot( R, rho_2D, 'r-', label = 'Sersic profile ($\\Sigma_{\\ast}$)')
    ax.plot( rx, rho_3D, 'b--', label = 'Deprojection of Sersic ($\\rho_{\\ast}$)',)
    ax.set_xscale('log')
    ax.set_xlabel('R [kpc]', fontsize = 12,)
    ax.set_xlim( 1e-1, 3e2 )

    ax.set_yscale('log')
    ax.set_ylabel('$\\rho_{\\ast} \, [M_{\\odot} \, / \, kpc^{3}]$ or ' + 
                '$\\Sigma_{\\ast} \, [M_{\\odot} \, / \, kpc^{2}]$', fontsize = 12,)
    ax.set_ylim( 1e3, 3e10 )

    ax.legend( loc = 1, frameon = False, fontsize = 12,)
    ax.tick_params( axis = 'both', which = 'both', direction = 'in', labelsize = 12,)

    plt.savefig('/home/xkchen/rho_star.png', dpi = 300)
    plt.close()


    ##. enclose mass test ?????
    scal_rho = rho_2D / Ie

    tck = interp.splrep( R, scal_rho, s = 0)

    xnew = np.logspace( -1, np.log10( R[-1] ), 600 )
    ynew = interp.splev( xnew, tck, der = 0 )


    order = 7
    cumu_M = np.zeros( N_ss,)

    for dd in range( N_ss ):

        cumu_f = GaussLegendreQuad_arr( R, scal_rho * R, order, 0, R[ dd ] )
        cumu_M[ dd ] = cumu_f[0] * Ie * 2 * np.pi

    plt.figure()
    plt.plot( R, cumu_M, 'r-', alpha = 0.5)
    plt.axvline( x = Re, ls = ':',)

    plt.axhline( 0.5 * 10**lg_Mstar, ls = '--',)
    plt.axhline( 10**lg_Mstar, ls = '-',)

    plt.xlabel('R [kpc]')
    plt.xscale('log')

    plt.ylabel('$\\Sigma_{\\ast} \; [M_{\\odot} / kpc^{2}]$')
    plt.yscale('log')

    plt.savefig('/home/xkchen/cumu_Mass_test.png', dpi = 300)
    plt.close()

