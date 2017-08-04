import numpy as np
import os, ctypes

modulepath=os.path.dirname(os.path.abspath(__file__))
libpath=os.path.join(modulepath, 'cBB6D.so')
BB6D=ctypes.CDLL(libpath)

from boost import ParBoost

# Python wrapper for C boost functions


def boost(x, px, y, py, sigma, delta, parboost):
    
    struct_to_pass = parboost.tobuffer()
    
    x_c = ctypes.c_double(x)
    px_c = ctypes.c_double(px)
    y_c = ctypes.c_double(y)
    py_c = ctypes.c_double(py)
    sigma_c = ctypes.c_double(sigma)
    delta_c = ctypes.c_double(delta)

    BB6D.BB6D_boost(struct_to_pass.ctypes.data, 
        ctypes.byref(x_c),
        ctypes.byref(px_c),
        ctypes.byref(y_c),
        ctypes.byref(py_c),
        ctypes.byref(sigma_c),
        ctypes.byref(delta_c))
        
    return (x_c.value,
        px_c.value,
        y_c.value,
        py_c.value,
        sigma_c.value,
        delta_c.value)
        
def inv_boost(x_st, px_st, y_st, py_st, sigma_st, delta_st, parboost):
    
    struct_to_pass = np.zeros(5,dtype=np.float64)
    struct_to_pass[0] = parboost.sphi
    struct_to_pass[1] = parboost.cphi
    struct_to_pass[2] = parboost.tphi
    struct_to_pass[3] = parboost.salpha
    struct_to_pass[4] = parboost.calpha
    
    x_c = ctypes.c_double(x_st)
    px_c = ctypes.c_double(px_st)
    y_c = ctypes.c_double(y_st)
    py_c = ctypes.c_double(py_st)
    sigma_c = ctypes.c_double(sigma_st)
    delta_c = ctypes.c_double(delta_st)
    
    BB6D.BB6D_inv_boost(struct_to_pass.ctypes.data, 
        ctypes.byref(x_c),
        ctypes.byref(px_c),
        ctypes.byref(y_c),
        ctypes.byref(py_c),
        ctypes.byref(sigma_c),
        ctypes.byref(delta_c))
        
    return (x_c.value,
        px_c.value,
        y_c.value,
        py_c.value,
        sigma_c.value,
        delta_c.value)
    
