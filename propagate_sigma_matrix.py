import numpy as np

class Sigmas(object):
    def __init__(self, Sig_11_0, Sig_12_0, Sig_13_0,
                Sig_14_0, Sig_22_0, Sig_23_0, Sig_24_0,
                Sig_33_0, Sig_34_0, Sig_44_0):
        
        self.Sig_11_0 = Sig_11_0
        self.Sig_12_0 = Sig_12_0
        self.Sig_13_0 = Sig_13_0
        self.Sig_14_0 = Sig_14_0
        self.Sig_22_0 = Sig_22_0
        self.Sig_23_0 = Sig_23_0
        self.Sig_24_0 = Sig_24_0
        self.Sig_33_0 = Sig_33_0
        self.Sig_34_0 = Sig_34_0
        self.Sig_44_0 = Sig_44_0

def propagate_Sigma_matrix(Sigmas_at_0, S, threshold_singular = 1e-30):
    
    Sig_11_0 = Sigmas_at_0.Sig_11_0
    Sig_12_0 = Sigmas_at_0.Sig_12_0
    Sig_13_0 = Sigmas_at_0.Sig_13_0
    Sig_14_0 = Sigmas_at_0.Sig_14_0
    Sig_22_0 = Sigmas_at_0.Sig_22_0
    Sig_23_0 = Sigmas_at_0.Sig_23_0
    Sig_24_0 = Sigmas_at_0.Sig_24_0
    Sig_33_0 = Sigmas_at_0.Sig_33_0
    Sig_34_0 = Sigmas_at_0.Sig_34_0
    Sig_44_0 = Sigmas_at_0.Sig_44_0
    
    Sig_11 = Sig_11_0 + 2.*Sig_12_0*S+Sig_22_0*S*S
    Sig_33 = Sig_33_0 + 2.*Sig_34_0*S+Sig_44_0*S*S
    Sig_13 = Sig_13_0 + (Sig_14_0+Sig_23_0)*S+Sig_24_0*S*S

    R = Sig_11-Sig_33
    W = Sig_11+Sig_33
    T = R*R+4*Sig_13*Sig_13
    
    if T<threshold_singular:
        a = Sig_12-Sig_34
        b = Sig_22-Sig_44
        c = Sig_14+Sig_23
        d = Sig_24
    
    else:
        sqrtT = np.sqrt(T)
        signR = np.sign(R)

        cos2theta = signR*R/sqrtT
        costheta = np.sqrt(0.5*(1.+cos2theta))
        sintheta = signR*np.sign(Sig_13)*np.sqrt(0.5*(1.-cos2theta))

        # in sixtrack this line seems to be different different
        #~ sintheta = -np.sign((Sig_11-Sig_33))*np.sqrt(0.5*(1.-cos2theta))

        Sig_11_hat = 0.5*(W+signR*sqrtT)
        Sig_33_hat = 0.5*(W-signR*sqrtT)

        #evaluate derivatives
        dS_R = 2.*(Sig_12_0-Sig_34_0)+2*S*(Sig_22_0-Sig_44_0)
        dS_W = 2.*(Sig_12_0+Sig_34_0)+2*S*(Sig_22_0+Sig_44_0)
        dS_Sig_13 = Sig_14_0 + Sig_23_0 + 2*Sig_24_0*S
        dS_T = 2*R*dS_R+8.*Sig_13*dS_Sig_13

        dS_cos2theta = signR*(dS_R/sqrtT - R/(2*sqrtT*sqrtT*sqrtT)*dS_T)
        dS_costheta = 1/(4*costheta)*dS_cos2theta
        dS_sintheta = -1/(4*sintheta)*dS_cos2theta

        dS_Sig_11_hat = 0.5*(dS_W + signR*0.5/sqrtT*dS_T)
        dS_Sig_33_hat = 0.5*(dS_W - signR*0.5/sqrtT*dS_T)
    
    # This will not be exposed in C
    extra_data = {}
    extra_data['Sig_11'] = Sig_11
    extra_data['Sig_13'] = Sig_13
    extra_data['Sig_33'] = Sig_33
    
    return Sig_11_hat, Sig_33_hat, costheta, sintheta,\
        dS_Sig_11_hat, dS_Sig_33_hat, dS_costheta, dS_sintheta,\
        extra_data


propagate_Sigma_matrix_vectorized = np.vectorize(propagate_Sigma_matrix, excluded =['Sigmas_at_0', 'threshold_singular'])

    
    
