"""
this file use to figure a '+'-like mice symbole at a given piont (in 2D figure)
"""
_all_ = ['mpl', 'plt', 'np', 'mice_symble']
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
def mice_symble(x0,y0,a0,l1,l2,l3,l4, **kwargs):
    """
    x0,y0 : the center of the mice symbole 
    l1,l2,l3,l4 : the length of the four sub-line-symbole of the mice symbole
    represent : right, up, left, bottom (1, 2, 3, 4)
    a0 : 1/2 of the width of the centeral square area
    """
    x0 = x0
    y0 = y0
    a0 = a0 
    lr = l1
    lu = l2
    ll = l3
    lb = l4
    # the two point of right
    xr1 = x0+a0
    yr1 = y0
    xr2 = x0+a0+lr
    yr2 = y0
    # the two point of upper
    xu1 = x0
    yu1 = y0+a0
    xu2 = x0
    yu2 = y0+a0+lu
    # the two point of left
    xl1 = x0-a0
    yl1 = y0
    xl2 = x0-a0-ll
    yl2 = y0
    # the two point of bottom
    xb1 = x0
    yb1 = y0-a0
    xb2 = x0
    yb2 = y0-a0-lb
    array_r = [[xr1,xr2],[yr1,yr2]]
    array_u = [[xu1,xu2],[yu1,yu2]]
    array_l = [[xl1,xl2],[yl1,yl2]]
    array_b = [[xb1,xb2],[yb1,yb2]]
    plt.plot(array_r[0],array_r[1],**kwargs)
    plt.plot(array_u[0],array_u[1],**kwargs)
    plt.plot(array_l[0],array_l[1],**kwargs)
    plt.plot(array_b[0],array_b[1],**kwargs)
    return
    