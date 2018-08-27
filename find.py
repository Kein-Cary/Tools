# this file try to find some element in an array
import numpy as np
import scipy.stats as st

_all_ = ['find1d','find1ds','find2d','find2ds']
def find1d(a,b):
    """
    This function is try to find an element in an array
    parameter:
    ----------
    a : the array use to find an element
    b : the element need to find
    return:
    -------
    the index of the element in the given array
    PS:if there are more than one element equal to b,the return value will be the first one
    """
    try :
        da = np.array(a)
        goal = b
        ix = da == goal
        iy = ix.tolist()
        x = iy.index(True)
    except ValueError:
        print('can not find element!')
        x = 'Error!'
    return x
#find1d = (a = True,b = True)
def find1ds(a,b):
    """
    This function will use to find elements those equal to the given value
    parameter:
    ----------
    a : the array use to find an element
    b : the element need to find(b not only one)
    return:
    -------
    data_ip : the index of the element in the given array,in type int 
    for index should be int type
    PS:this function also can be use to find element like NaN,inf,but 
    this part can't be used to find boor type
    """
    try:
        da = np.array(a)
        goal = b
        f = st.find_repeats(da)
        ix = f[0] == goal
        iy = ix.tolist()
        x = iy.index(True)
        y = f[1][x]
        data_ip = np.zeros(y,dtype = np.int0)
        for k in range(y):
            ia = da == goal
            ib = ia.tolist()
            ic = ib.index(True)
            data_ip[k] = ic
            if goal == 0 :
                da[ic] = np.int0(True)
            else:
                da[ic] = np.int0(False)
    except ValueError:
        print('can not find element!')  
        data_ip = 'Error!'         
    return data_ip
##find1ds = (a = True,b = True)
def find2d(a,b):
    """
    This function is used to find a element in a 2D array
    parameter:
    ----------
    a : the array for finding the element
    b : the element need to find
    return :
    --------
    data_ip : the index of the element in the array
    """
    try:
        da = np.array(a)
        goal = b
        A = da.flatten()
        B = da.shape
        ix = A == goal
        iy = ix.tolist()
        x = iy.index(True)
        N = np.floor((x+1)/B[1])
        if (x+1) % B[1] == 0:
            X = N - 1
        else:
            X = N
        Y = x % B[1]
        data_ip = np.array([X,Y])
    except ValueError:
        print('can not find element!')
        data_ip = 'Error!'
    return data_ip
#find2d(a = True,b = True)
def find2ds(a,b):
    """
    This function is used to find element in a 2D array
    parameter:
    ----------
    a : the array in which will to find the element
    b : the elements need to find
    return:
    ---------
    data_ip:give the index all the elements in given array 
    PS :this part can't be used to find boor type
    """
    try:
        da = np.array(a)
        goal = b
        f = st.find_repeats(da)
        ix = f[0] == goal
        iy = ix.tolist()
        x = iy.index(True)
        y = f[1][x]
        data_ip = np.zeros((y,2),dtype = np.int0)
        A = da.flatten()
        B = da.shape
        for k in range(y):
            ia = A == goal
            ib = ia.tolist()
            ip = ib.index(True)
            N = np.floor((ip+1)/B[1])
            if (ip+1) % B[1] == 0:
                X = N - 1
            else:
                X = N
            Y = ip % B[1]
            data_ip[k,0] = X
            data_ip[k,1] = Y
            if goal == 0:
                A[ip] = np.int0(True)
            else:
                A[ip] = np.int0(False)
    except ValueError:
        print('can not find element!')
        data_ip = 'Error!'
    return data_ip
#find2ds(a = True,b = True)