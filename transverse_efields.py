import ctypes
import os
import numpy as np

modulepath=os.path.dirname(os.path.abspath(__file__))
libpath=os.path.join(modulepath, 'cFromsixtracklib.so')
fromsixtracklib=ctypes.CDLL(libpath)

def transverse_field_ellip(x, y, sigma_x, sigma_y, Delta_x, Delta_y):

    Ex_out = ctypes.c_double()
    Ey_out = ctypes.c_double()

    struct_to_pass = np.zeros(4,dtype=np.float64)
    struct_to_pass[0] = sigma_x
    struct_to_pass[1] = sigma_y
    struct_to_pass[2] = Delta_x
    struct_to_pass[3] = Delta_y

    fromsixtracklib.get_transv_field_gauss_ellip(struct_to_pass.ctypes.data, 
                            ctypes.c_double(x), ctypes.c_double(y), 
                            ctypes.byref(Ex_out), ctypes.byref(Ey_out))
                            
    return Ex_out.value, Ey_out.value
vect_transverse_field_ellip = np.vectorize(transverse_field_ellip)

def transverse_field_round(x, y, sigma, Delta_x, Delta_y):

    Ex_out = ctypes.c_double()
    Ey_out = ctypes.c_double()

    struct_to_pass = np.zeros(3,dtype=np.float64)
    struct_to_pass[0] = sigma
    struct_to_pass[1] = Delta_x
    struct_to_pass[2] = Delta_y

    fromsixtracklib.get_transv_field_gauss_round(struct_to_pass.ctypes.data, 
                            ctypes.c_double(x), ctypes.c_double(y), 
                            ctypes.byref(Ex_out), ctypes.byref(Ey_out))
                            
    return Ex_out.value, Ey_out.value
vect_transverse_field_round = np.vectorize(transverse_field_round)


