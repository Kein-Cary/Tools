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


### === special string files query

import glob

files = glob.glob('/home/xkchen/mywork/Sat_SB/code/*/*.py')
N_file = len( files )

tag_str = '/home/xkchen/project/tmp_obj_cat'
out_str = []

for ll in range( N_file ):

    dat = open( files[ ll ] )

    lines = dat.readlines()
    N_lines = len( lines )

    for dd in range( N_lines ):

        if tag_str in lines[ dd ]:

            out_str.append( files[ ll ] )

            break

        else:
            continue

print( out_str )

