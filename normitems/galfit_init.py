"""
this file use to build the 'initial feed parameters' of galfit
"""
from io import StringIO

import subprocess as subpro
import os


##.
dps = \
"""
===============================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) %s            # Input data image (FITS file)
B) %s       # Output data image block
C) none                # Sigma image name (made from data if blank or "none") 
D) %s        # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data 
F) none                # Bad pixel mask (FITS image or ASCII coord list)
G) none                # File with parameter constraints (ASCII file) 
H) %d    %d   %d    %d   # Image region to fit (xmin xmax ymin ymax)
I) %d    %d          # Size of the convolution box (x y)
J) %.2f              # Magnitude photometric zeropoint 
K) %.3f  %.3f        # Plate scale (dx dy)    [arcsec per pixel]
O) regular             # Display type (regular, curses, both)
P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

# INITIAL FITTING PARAMETERS
#
#   For object type, the allowed functions are: 
#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat, 
#       ferrer, powsersic, sky, and isophote. 
#  
#   Hidden parameters will only appear when they're specified:
#       C0 (diskyness/boxyness), 
#       Fn (n=integer, Azimuthal Fourier Modes),
#       R0-R10 (PA rotation, for creating spiral structures).
# 
# -----------------------------------------------------------------------------
#   par)    par value(s)    fit toggle(s)    # parameter description 
# -----------------------------------------------------------------------------

# Object number: 1
 0) sersic                 #  object type
 1) %.1f  %.1f  1 1  #  position x, y
 3) 20.0890     1          #  Integrated magnitude	
 4) 5.1160      1          #  R_e (half-light radius)   [pix]
 5) 4.2490      1          #  Sersic index n (de Vaucouleurs n=4) 
 6) 0.0000      0          #     ----- 
 7) 0.0000      0          #     ----- 
 8) 0.0000      0          #     ----- 
 9) 0.7570      1          #  axis ratio (b/a)  
10) 45    1          # -60.3690,  position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 

# Object number: 2
 0) sky                    #  object type
 1) %.3f        1          #  sky background at center of fitting region [ADUs]
 2) 0.0000      0          #  dsky/dx (sky gradient in x)
 3) 0.0000      0          #  dsky/dy (sky gradient in y)
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 

================================================================================
""" 

# % ( fit_file, out_file, psf_file, 
# 		xmin, xmax, ymin, ymax, 
# 		covx, covy, 
# 		ZPts, 
# 		pixsize, pixsize )

# docs = '/home/xkchen/test.feedfit'

# print( dps, file = docs )

# docs.close()

