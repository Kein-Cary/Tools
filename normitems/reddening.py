#Junde Chen @ SJTU 2019
import numpy as np

def _isiterable(obj):
    """
    check if object is iterable
    """
    try:
        iter(obj)
        return True
    except TypeError:
        return False

def ccm89(wave,Av,Rv=3.1):
    """
    Cardelli, Clayton & Mathis (1989) extinction function law.

    Parameters:
    ===========
    wave : wavelength in Angstroms
    Av : v-band extinction
    Rv : extinction coefficient, Av/E(B-V), 3.1 for milky way
    
    Returns:
    ========
    extinctions at each input wavelength
    
    """
    x = 10**4/wave #convert Angstrom to um^-1
    if _isiterable(x):
        return np.array([_ccm89_module(i,Rv) for i in x])*Av
    else:
        return _ccm89_module(x,Rv)*Av
    
def _ccm89_module(x,Rv):

    if 0.3<=x<1.1:
        """
        infrared
        """
        a = 0.574*x**1.61
        b = -0.527*x**1.61

    elif 1.1<=x<3.3:
        """
        Optical/NIR
        """
        y = x - 1.82
        a = (1. + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 
             + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7)
        b = (1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4
             - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7)

    elif 3.3<=x<8:
        """
        Ultraviolet
        """
        if x<5.9:
            Fa = 0.
            Fb = 0.
        elif 5.9<=x<8.:
            Fa = -0.04473*(x - 5.9)**2 - 0.009779*(x - 5.9)**3
            Fb = 0.2130*(x - 5.9)**2 + 0.1207*(x - 5.9)**3
            
        a = 1.752 - 0.316*x - 0.104/((x - 4.67)**2 + 0.341) + Fa
        b = -3.090 + 1.825*x + 1.206/((x - 4.62)**2 + 0.263) + Fb
        
    elif 8.<=x<=10.:
        """
        Far-UV
        """
        a = -1.073 - 0.628*(x - 8) + 0.137*(x - 8)**2 - 0.070*(x - 8)**3
        b = 13.670 + 4.257*(x - 8) - 0.420*(x - 8)**2 + 0.374*(x - 8)**3
        
    return a+b/Rv

#################################################################################
    
def calzetti00(wave,Av,Rv):
    """
    Calzetti et al (2000) extinction function law.

    Parameters:
    ===========
    wave : wavelength in Angstroms
    Av : v-band extinction
    Rv : extinction coefficient, Av/E(B-V), 4.05 for starburst galaxies
    
    Returns:
    ========
    extinctions at each input wavelength
    
    """
    x = wave/10**4 #convert Angstrom to um
    if _isiterable(x):
        return np.array([_calzetti00_module(i,Rv) for i in x])*Av/Rv
    else:
        return _calzetti00_module(x,Rv)*Av/Rv

def _calzetti00_module(x,Rv):
    if 0.63<=x:
        return 2.659*(-1.857+1.040/x)+Rv
    elif x<0.63:
        return 2.659*(-2.156+1.509/x-0.198/x/x + 0.011/x/x/x)+Rv
    
#################################################################################

def red_corr_ccm89(wave,flux,oc_ratio,c1,c2,Rv=3.1):
    """
    correct intrinsic reddening using Calzetti00 extinction curve

    Parameters:
    ===========
    wave :  input wavelength in Angstrom
    flux :  the observed flux at input wave
    oc_ratio=r_cri/r_obs : r_obs: flux ratio of line1 and line2, r_cri: intrinsic value of line1/line2
    c1,c2 : line center of 1 and c2
    Rv : extinction coefficient, Av/E(B-V), 4.05 for starburst galaxies

    Returns:
    ===========
    corrected flux
    """
    x = wave
    k1 = ccm89(c1,1.0,Rv)
    k2 = ccm89(c2,1.0,Rv)
    kl = ccm89(x ,1.0,Rv)
    Al = 2.5*kl/(k1-k2)*np.log10(oc_ratio)
    return flux*10**(0.4*Al)

def red_corr_calzetti00(wave,flux,oc_ratio,c1,c2,Rv=4.05):
    """
    correct intrinsic reddening using Calzetti00 extinction curve

    Parameters:
    ===========
    wave :  input wavelength in Angstrom
    flux :  the observed flux at input wave
    oc_ratio=r_cri/r_obs : r_obs: flux ratio of line1 and line2, r_cri: intrinsic value of line1/line2
    c1,c2 : line center of 1 and c2
    Rv : extinction coefficient, Av/E(B-V), 4.05 for starburst galaxies

    Returns:
    ===========
    corrected flux
    """
    x = wave
    k1 = calzetti00(c1,1.0,Rv)
    k2 = calzetti00(c2,1.0,Rv)
    kl = calzetti00(x ,1.0,Rv)
    Al = 2.5*kl/(k1-k2)*np.log10(oc_ratio)
    return flux*10**(0.4*Al)
