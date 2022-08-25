"""
this file used to find the angular diameter in different redshift
"""
import astropy.constants as C
import astropy.units as U
import numpy as np
from astropy import cosmology as apcy
from scipy.integrate import quad as scq


## constant
rad2arcsec = U.rad.to(U.arcsec)
kpc2m = U.kpc.to(U.m)
Mpc2cm = U.Mpc.to(U.cm)

L_speed = C.c.value # m/s
G = C.G.value

M_sun = C.M_sun.value # in unit of kg
Lsun = C.L_sun.value*10**7 # (erg/s/cm^2)
Msun2kg = U.M_sun.to(U.kg)


# get the velocity of light and in unit Km/s
vc = C.c.to(U.km/U.s).value


##. initial parameters
def input_cosm_model( get_model = None ):

    global Test_model

    if get_model is not None:
        Test_model = get_model

    else:
        ### cosmology
        Test_model = apcy.FlatwCDM( H0 = 70.0, Om0 = 0.3, Ob0 = 0.0, w0 = -1.0 )

    return Test_model

def cosmos_param():

    global H0, h, Omega_m, Omega_lambda, Omega_k, DH

    ## cosmology params
    H0 = Test_model.H0.value
    h = H0 / 100
    Omega_m = Test_model.Om0
    Omega_lambda = 1.-Omega_m
    Omega_k = 1.- (Omega_lambda + Omega_m)

    DH = vc / H0

    return


def mark_by_self(z_in,goal_size):

    z_array = z_in

    size_cluster = goal_size

    Alpha = np.zeros(len(z_array), dtype = np.float)

    Da = np.zeros(len(z_array), dtype = np.float)

    Dc = np.zeros(len(z_array), dtype = np.float)

    fEz1 = lambda sysz: 1./np.sqrt(Omega_m*(1+sysz)**3+Omega_lambda+Omega_k*(1+sysz)**2)

    for k in range(len(z_array)):

        if z_array[k] == 0:
            Da[k] = 0.
        else:
            Dc[k] = DH * scq(fEz1,0.,z_array[k])[0]

    if Omega_k == 0.:
        Dm = Dc*1.

    elif Omega_k > 0.:
        Dm = DH*np.sinh(np.asrt(Omega_k)*Dc/DH)*1/np.sqrt(Omega_k)

    else:

        Dm = DH*np.sin(np.asrt(np.abs(Omega_k))*Dc/DH)*1/np.sqrt(np.abs(Omega_k))

    DA = Dm / ( 1 + z_array )
    alpha = ( size_cluster / h ) / DA
    Alpha = alpha
    Da = DA

    return Alpha, Da


def mark_by_plank(z_in,goal_size):
    z = z_in
    size = goal_size
    Alpha = np.zeros(len(z), dtype = np.float)
    Da = np.zeros(len(z), dtype = np.float)
    for k in range(len(z)):
        if z[k] == 0:
            Da[k] = 0.
            Alpha[k] = np.inf
        else:
            DA_conference = Test_model.angular_diameter_distance(z[k]).value
            alpha_conference = (size/h)/DA_conference
            Alpha[k] = alpha_conference
            Da[k] = DA_conference
    return Alpha, Da


def mark_by_self_Noh(z_in,goal_size):
    z_array = z_in
    size_cluster = goal_size
    Alpha = np.zeros(len(z_array), dtype = np.float)
    Da = np.zeros(len(z_array), dtype = np.float)
    Dc = np.zeros(len(z_array), dtype = np.float)
    fEz1 = lambda sysz: 1./np.sqrt(Omega_m*(1+sysz)**3+Omega_lambda+Omega_k*(1+sysz)**2)
    for k in range(len(z_array)):
        if z_array[k] == 0:
            Da[k] = 0.
        else:
            Dc[k] = DH*scq(fEz1,0.,z_array[k])[0]
    if Omega_k == 0.:
        Dm = Dc*1.
    elif Omega_k > 0.:
        Dm = DH*np.sinh(np.asrt(Omega_k)*Dc/DH)*1/np.sqrt(Omega_k)
    else:
        Dm = DH*np.sin(np.asrt(np.abs(Omega_k))*Dc/DH)*1/np.sqrt(np.abs(Omega_k))
    DA = Dm/(1+z_array)
    alpha = size_cluster/DA
    Alpha = alpha
    Da = DA 
    return Alpha, Da


def mark_by_plank_Noh(z_in,goal_size):
    z = z_in
    size = goal_size
    Alpha = np.zeros(len(z), dtype = np.float)
    Da = np.zeros(len(z), dtype = np.float)
    for k in range(len(z)):
        if z[k] == 0:
            Da[k] = 0.
            Alpha[k] = np.inf
        else:
            DA_conference = Test_model.angular_diameter_distance(z[k]).value
            alpha_conference = size/DA_conference
            Alpha[k] = alpha_conference
            Da[k] = DA_conference

    return Alpha, Da


##.
if __name__ == "__main__":

    input_cosm_model( get_model = None )
    cosmos_param()

    print( H0, h, Omega_m, Omega_lambda, Omega_k, DH )

