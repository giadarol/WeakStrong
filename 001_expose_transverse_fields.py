import ctypes
import os
import numpy as np

modulepath=os.path.dirname(os.path.abspath(__file__))
libpath=os.path.join(modulepath, 'cFromsixtracklib.so')
fromsixtracklib=ctypes.CDLL(libpath)

#def transverse_field_ellip(x, y, sigma_x, sigma_y, Delta_x, Delta_y):
    
x = 0.3
y = 0.3
sigma_x=.9
sigma_y=.9000001
Delta_x = 0.
Delta_y = 0.

Ex_out = ctypes.c_double()
Ey_out = ctypes.c_double()

struct_to_pass = np.zeros(4,dtype=np.float64)
struct_to_pass[0] = sigma_x
struct_to_pass[1] = sigma_y
struct_to_pass[2] = Delta_x
struct_to_pass[3] = Delta_y

fromsixtracklib.get_transv_field_gauss_ellip(struct_to_pass.ctypes.data, 
                        ctypes.c_double(x), ctypes.c_double(x), 
                        ctypes.byref(Ex_out), ctypes.byref(Ey_out))
                        
print '%e, %e'%(Ex_out.value, Ey_out.value)

