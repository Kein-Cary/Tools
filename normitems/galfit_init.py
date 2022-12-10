"""
this file use to build the 'initial feed parameters' of galfit
"""
from io import StringIO

import subprocess as subpro
import os


##.
dps_wPSF = \
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

# Sersic function

 0) sersic             # Object type
 1) %d %d      1 1    # position x, y        [pixel]
 3) %.2f      1       # total magnitude    
 4) 4.30       1       #     R_e              [Pixels]
 5) 5.20       1       # Sersic exponent (deVauc=4, expdisk=1)  
 9) 0.30       1       # axis ratio (b/a)   
10) 10.0       1       # position angle (PA)  [Degrees: Up=0, Left=90]
 Z) 0                  #  Skip this model in output image?  (yes=1, no=0)

================================================================================
""" 



##.
dps_woPSF = \
"""
===============================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) %s            # Input data image (FITS file)
B) %s       # Output data image block
C) none                # Sigma image name (made from data if blank or "none") 
D) none                # Input PSF image and (optional) diffusion kernel
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

# Sersic function

 0) sersic             # Object type
 1) %d %d      1 1    # position x, y        [pixel]
 3) %.2f      1       # total magnitude    
 4) 4.30       1       #     R_e              [Pixels]
 5) 5.20       1       # Sersic exponent (deVauc=4, expdisk=1)  
 9) 0.30       1       # axis ratio (b/a)   
10) 10.0       1       # position angle (PA)  [Degrees: Up=0, Left=90]
 Z) 0                  #  Skip this model in output image?  (yes=1, no=0)

================================================================================
""" 


###... others model

# Nuker function
extra_model = \
"""
 0) nuker              # Object type
 1) 250.  475.  1 1    # position x, y        [pixel]
 3) 17.2       1       #    mu(Rb)            [surface brightness mag. at Rb]
 4) 20.5       1       #     Rb               [pixels]
 5) 1.2        1       #    alpha  (sharpness of transition)
 6) 0.5        1       #    beta   (outer powerlaw slope)
 7) 0.7        1       #    gamma  (inner powerlaw slope)
 9) 0.72       1       # axis ratio (b/a)   
10) -25.2      1       # position angle (PA)  [Degrees: Up=0, Left=90]
 Z) 0                  #  Skip this model in output image?  (yes=1, no=0)

# Moffat function
 
 0) moffat             # object type
 1) 372.0  450.0 1 1   # position x, y        [pixel]
 3) 16.5       1       # total magnitude     
 4) 0.5        1       #   FWHM               [Pixels]
 5) 1.5        1       # powerlaw      
 9) 0.3        1       # axis ratio (b/a)   
10) 25         1       # position angle (PA)  [Degrees: Up=0, Left=90]
 Z) 0                  #  Skip this model in output image?  (yes=1, no=0)

# Gaussian function

 0) gaussian           # object type
 1) 402.3  345.9  1 1  # position x, y        [pixel]
 3) 18.5       1       # total magnitude     
 4) 0.5        0       #   FWHM               [pixels]
 9) 0.3        1       # axis ratio (b/a)   
10) 25         1       # position angle (PA)  [Degrees: Up=0, Left=90]
 Z) 0                  # leave in [1] or subtract [0] this comp from data?

"""
