#  this file use to calculate this type formulation: y = x!; y = x!!
import numpy as np

_all_ = ['n_step', 'n_n_step']

def n_step(x):
    """
    x : the data will use to calculate x!, x can be a data, or an array
    """
    x = x
    try:
        gama = np.zeros(len(x), dtype = np.float)
        for k in range(len(x)):
            gama[k] = 1
            for p in range(1,x[k]+1):
                gama[k] = gama[k]*p
    except TypeError:
        gama = 1
        for k in range(1,x+1):
            gama = gama*k
    return gama

def n_n_step(x):
    """
    x : the data will use to calculate x!!, x can be a data, or an array
    """
    x = x
    try:
        gama = np.zeros(len(x),dtype = np.float)
        for k in range(len(x)):
            gama[k] = 1
            if x[k] % 2 == 0:
                for p in range(1,np.int0(x/2+1)):
                    gama[k] = gama[k]*(p*2)
            if x[k] % 2 == 1:
                for p in range(1,np.int0((x+1)/2)):
                    gama[k] = gama[k]*(2*p+1)  
    except TypeError: 
        if x % 2 == 0:
            gama = 1
            for k in range(1,np.int0(x/2+1)):
                gama = gama*(k*2)
        if x % 2 == 1:
            gama = 1
            for k in range(1,np.int0((x+1)/2)):
                gama = gama*(2*k+1)         
    return gama
