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
extra_model = \
"""
 0) sersic      #  object type
 1) 48.5180  51.2800  1 1  #  position x, y
 3) 20.0890     1          #  Integrated magnitude	
 4) 5.1160      1          #  R_e (effective radius)   [pix]
 5) 4.2500      0          #  index of Sersic 
 6) 0.0000      0          #     ----- 
 7) 0.0000      0          #     ----- 
 8) 0.0000      0          #     ----- 
 9) 0.7570      1          #  axis ratio (b/a)  
10) -60.3690    1          #  position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 


 0) expdisk                #  object type
 1) 48.5180  51.2800  1 1  #  position x, y
 3) 20.0890     1          #  Integrated magnitude	
 4) 5.1160      1          #  R_s (scale length)   [pix]
 5) 0.0000      0          #     ----- 
 6) 0.0000      0          #     ----- 
 7) 0.0000      0          #     ----- 
 8) 0.0000      0          #     ----- 
 9) 0.7570      1          #  axis ratio (b/a)  
10) -60.3690    1          #  position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 


 0) gaussian               #  object type
 1) 48.5180  51.2800  1 1  #  position x, y
 3) 20.0890     1          #  Integrated magnitude	
 4) 5.1160      1          #  FWHM   [pix]
 5) 0.0000      0          #     ----- 
 6) 0.0000      0          #     ----- 
 7) 0.0000      0          #     ----- 
 8) 0.0000      0          #     ----- 
 9) 0.7570      1          #  axis ratio (b/a)  
10) -60.3690    1          #  position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract)


 0) ferrer                 #  object type
 1) 48.5180  51.2800  1 1  #  position x, y
 3) 18.0890     1          #  the surface brightness 	
 4) 5.1160      1          #  R_out  [pix]
 5) 2.0500      1          #  index alpha
 6) 1.250       0          #  index belta 
 7) 0.0000      0          #     ----- 
 8) 0.0000      0          #     ----- 
 9) 0.7570      1          #  axis ratio (b/a)  
10) -60.3690    1          #  position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 



 0) moffat                 #  object type
 1) 48.5180  51.2800  1 1  #  position x, y
 3) 20.0890     1          #  total magnitude 	
 4) 5.1160      1          #  FHWM  [pix]
 5) 2.0500      1          #  index n
 6) 0.0000      0          #     ----- 
 7) 0.0000      0          #     ----- 
 8) 0.0000      0          #     ----- 
 9) 0.7570      1          #  axis ratio (b/a)  
10) -60.3690    1          #  position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 


 0) king                   #  object type
 1) 48.5180  51.2800  1 1  #  position x, y
 3) 18.0890     1          #  the surface brightness 	
 4) 5.1160      1          #  R_c    [pix]
 5) 10.050      1          #  R_t    [pix] 
 6) 1.250       0          #  index alpha 
 7) 0.0000      0          #     ----- 
 8) 0.0000      0          #     ----- 
 9) 0.7570      1          #  axis ratio (b/a)  
10) -60.3690    1          #  position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 


 0) nuker                  #  object type
 1) 48.5180  51.2800  1 1  #  position x, y
 3) 18.0890     1          #  the surface brightness 	
 4) 5.1160      1          #  R_b   [pix]
 5) 1.2490      1          #  power index alpha 
 6) 2.2350      0          #  power index belta 
 7) 0.5750      0          #  power index gamma 
 8) 0.0000      0          #     ----- 
 9) 0.7570      1          #  axis ratio (b/a)  
10) -60.3690    1          #  position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 

"""

