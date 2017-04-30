import ctypes
import os
import numpy as np
from scipy.constants import epsilon_0

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
    
    
def get_Ex_Ey_Gx_Gy_gauss(x, y, sigma_x, sigma_y, min_sigma_diff):

    flag_round = np.abs(sigma_x-sigma_y)< min_sigma_diff
    
    if flag_round:

        sigma = 0.5*(sigma_x+sigma_y)
        
        #electric fields from sixtracklib
        Ex, Ey = transverse_field_round(x, y, sigma, Delta_x=0., Delta_y=0.)
        Gx = 1/(2.*(x*x+y*y))*(y*Ey-x*Ex+1./(2*np.pi*epsilon_0*sigma*sigma)*x*x*np.exp(-(x*x+y*y)/(2.*sigma*sigma)))
        Gy = 1./(2*(x*x+y*y))*(x*Ex-y*Ey+1./(2*np.pi*epsilon_0*sigma*sigma)*y*y*np.exp(-(x*x+y*y)/(2.*sigma*sigma)))

    else:
        Ex, Ey = transverse_field_ellip(x, y, sigma_x, sigma_y, Delta_x=0., Delta_y=0.)
        Sig_11 = sigma_x*sigma_x
        Sig_33 = sigma_y*sigma_y
        
        Gx =-1./(2*(Sig_11-Sig_33))*(x*Ex+y*Ey+1./(2*np.pi*epsilon_0)*\
                    (sigma_y/sigma_x*np.exp(-x*x/(2*Sig_11)-y*y/(2*Sig_33))-1.))
        Gy =1./(2*(Sig_11-Sig_33))*(x*Ex+y*Ey+1./(2*np.pi*epsilon_0)*\
                    (sigma_x/sigma_y*np.exp(-x*x/(2*Sig_11)-y*y/(2*Sig_33))-1.))
                    
    return Ex, Ey, Gx, Gy




