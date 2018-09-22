"""
# this file try to get the elimation of a number ,for example:
1.223 --> 1.22
1.225 --> 1.23
"""
import numpy as np
def get_round_number_float(value,n):
    """
    value ： the data to calculation
    n : the number of decimal to save the form
    """
    a = value*1
    N = n*1
    A = np.round(a*10**N)/10**N
    return A
#get_round_number_float(value = True,n = True)
def get_round_number_int(value,n):
    """
    value ： the data to calculation
    n : the number of digits from right to left
    This part is also can use to chang a foat number to a round number
    for instance : 234.567 --> 230,if n == 1,234.567 --> 235,if n == 0
    """
    a = value*1
    N = n*1
    A = np.round(a/10**N)*10**N
    return A
#get_round_number_int(value = True,n = True)
