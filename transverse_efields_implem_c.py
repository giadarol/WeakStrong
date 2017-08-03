import ctypes
import os
import numpy as np
from scipy.constants import epsilon_0

modulepath=os.path.dirname(os.path.abspath(__file__))
libpath=os.path.join(modulepath, 'cBB6D.so')
BB6D=ctypes.CDLL(libpath)

    
    
def get_Ex_Ey_Gx_Gy_gauss(x, y, sigma_x, sigma_y, min_sigma_diff):

    Ex = ctypes.c_double()
    Ey = ctypes.c_double()
    Gx = ctypes.c_double()
    Gy = ctypes.c_double()

    BB6D.get_Ex_Ey_Gx_Gy_gauss(ctypes.c_double(x), ctypes.c_double(y), 
        ctypes.c_double(sigma_x), ctypes.c_double(sigma_y), 
        ctypes.c_double(min_sigma_diff),
        ctypes.byref(Ex), ctypes.byref(Ey), 
        ctypes.byref(Gx), ctypes.byref(Gy))
                    
    return Ex.value, Ey.value, Gx.value, Gy.value




