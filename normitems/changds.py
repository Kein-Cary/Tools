import numpy as np
_all_ = ['chidas_int','chidas_float','chidas_array']

def chidas_int(data,length):
    """
    chang a data to a str,while data is in type int
    This may main use in cite the glob of files or documents
    parameter:
    ----------
    data : the data will be change to str
    length: the length of the str 
    
    return:
    -------
    s : the result of change str
    Ps: this use to fill zeros before the data,and if the data if float,see chidas_float
        while array may see chidas_array
    There's no type set for str,and you need to change this by yourself!
    """
    number = data
    sl = length
    ss = str(number)
    len_ss = len(ss)
    delta = sl - len_ss
    if delta < 0.:
       print('Input error,can not change!')
       s = []
    elif delta ==0:
       s = str(number)
    else:
       s = str(number).zfill(sl)
    return s
#chdas_int(data = True,length = True)
def chidas_float(data,length):
    """
    change a float data to str,the result depend on the data,it may be catious about negative data
    parameter:
    ----------
    data : the data will be change to str
    length : the length of the str
    PS : if the length is smaller than len(str(data)),then will be return str(data)
    
    return :
    --------
    s : the result of chang str
    Ps: this use to fill zeros before the data,and if the data if int,see chidas_int
        while array may see chidas_array
    There's no type set for str,and you need to change this by yourself!
    """
    number = data
    sl = length
    s = str(number).zfill(sl)
    return s
#chdas_float(data = True,length = True)
def chidas_array(data,length):
    """
    This file is using to change a np.array element to a 'table' which in element as str
    parameter :
    -----------
    data : the array will be change
    length : the len of str for each element in the goal table
    Ps:you must pay caution to that this function can use only for 1-D array,if len(data.shape)>1,
    please use a repeat way for different dimension
    For signal data,you may need chidas_int or chidas_float
    """
    a = data
    s = {}
    sl = length
    N = np.int0(a.shape[0])
    for k in range(N):
        s[k] = str(a[k]).zfill(sl)
    return s
#chidas_array(data = True,length = True)
def inv_chidas(s):
    """
    change str to data
    parameter:
    ---------
    s : the str need to change
    use like inv_chidas(s),and s need to '0123' and so on,also can be b'123' type
    """
    chs = s
    try:
        a = int(chs)
    except ValueError:
        a = float(chs)
    return a
#inv_chidas(s =True)
            
            
    