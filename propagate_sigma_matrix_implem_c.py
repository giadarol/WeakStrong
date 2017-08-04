import numpy as np

import os, ctypes

modulepath=os.path.dirname(os.path.abspath(__file__))
libpath=os.path.join(modulepath, 'cBB6D.so')
BB6D=ctypes.CDLL(libpath)


from propagate_sigma_matrix import Sigmas, boost_sigmas, propagate_full_Sigma_matrix_in_drift

# Python wrapper for the corresponding C function
def propagate_Sigma_matrix(Sigmas_at_0, S, threshold_singular = 1e-16, handle_singularities=True):
    
    struct_to_pass = np.zeros(10,dtype=np.float64)
    struct_to_pass[0] = Sigmas_at_0.Sig_11_0
    struct_to_pass[1] = Sigmas_at_0.Sig_12_0
    struct_to_pass[2] = Sigmas_at_0.Sig_13_0
    struct_to_pass[3] = Sigmas_at_0.Sig_14_0
    struct_to_pass[4] = Sigmas_at_0.Sig_22_0
    struct_to_pass[5] = Sigmas_at_0.Sig_23_0
    struct_to_pass[6] = Sigmas_at_0.Sig_24_0
    struct_to_pass[7] = Sigmas_at_0.Sig_33_0
    struct_to_pass[8] = Sigmas_at_0.Sig_34_0
    struct_to_pass[9] = Sigmas_at_0.Sig_44_0
    
    Sig_11_hat = ctypes.c_double()
    Sig_33_hat = ctypes.c_double()
    costheta = ctypes.c_double()
    sintheta = ctypes.c_double()
    dS_Sig_11_hat = ctypes.c_double()
    dS_Sig_33_hat = ctypes.c_double()
    dS_costheta = ctypes.c_double()
    dS_sintheta = ctypes.c_double()
    
    
    BB6D.BB6D_propagate_Sigma_matrix(struct_to_pass.ctypes.data,
        ctypes.c_double(S), ctypes.c_double(threshold_singular),
        ctypes.c_long(handle_singularities),
        ctypes.byref(Sig_11_hat), ctypes.byref(Sig_33_hat), 
        ctypes.byref(costheta), ctypes.byref(sintheta),
        ctypes.byref(dS_Sig_11_hat), ctypes.byref(dS_Sig_33_hat), 
        ctypes.byref(dS_costheta), ctypes.byref(dS_sintheta))
        
    extra_data={}
    
    return Sig_11_hat.value, Sig_33_hat.value, costheta.value, sintheta.value,\
        dS_Sig_11_hat.value, dS_Sig_33_hat.value, dS_costheta.value, dS_sintheta.value,\
        extra_data
        
propagate_Sigma_matrix_vectorized = np.vectorize(propagate_Sigma_matrix, excluded =['Sigmas_at_0', 'threshold_singular', 'handle_singularities'])
