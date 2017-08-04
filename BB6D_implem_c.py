import os, ctypes
modulepath=os.path.dirname(os.path.abspath(__file__))
libpath=os.path.join(modulepath, 'cBB6D.so')
BB6D=ctypes.CDLL(libpath)

import numpy as np
from BB6D import BB6D_Data, BB6D_init

    
def BB6D_track(x, px, y, py, sigma, delta, q0, p0, bb6ddata):
    
    x_c = ctypes.c_double(x)
    px_c = ctypes.c_double(px)
    y_c = ctypes.c_double(y)
    py_c = ctypes.c_double(py)
    sigma_c = ctypes.c_double(sigma)
    delta_c = ctypes.c_double(delta)

    buf = bb6ddata.tobuffer()

    BB6D.BB6D_track(ctypes.byref(x_c),
            ctypes.byref(px_c),
            ctypes.byref(y_c),
            ctypes.byref(py_c),
            ctypes.byref(sigma_c),
            ctypes.byref(delta_c), 
            ctypes.c_double(q0), 
            ctypes.c_double(p0), 
            buf.ctypes.data)
            
    (x_ret, px_ret, y_ret, py_ret, sigma_ret, delta_ret) = (x_c.value,
            px_c.value,
            y_c.value,
            py_c.value,
            sigma_c.value,
            delta_c.value) 
    
    return x_ret, px_ret, y_ret, py_ret, sigma_ret, delta_ret
