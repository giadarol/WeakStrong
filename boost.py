import numpy as np

# I program as close as possible to C...

class ParBoost(object):
    #it is practically a struct
    def __init__(self, phi, alpha):
        self.sphi = np.sin(phi)
        self.cphi = np.cos(phi)
        self.tphi = np.tan(phi)
        self.salpha = np.sin(alpha)
        self.calpha = np.cos(alpha)

def boost(x, px, y, py, sigma, delta, parboost):
    
    sphi = parboost.sphi
    cphi = parboost.cphi
    tphi = parboost.tphi
    salpha = parboost.salpha
    calpha = parboost.calpha
    
    h = delta + 1. - np.sqrt((1.+delta)**2-px**2-py**2) 

    px_st = px/cphi-h*calpha*tphi/cphi
    py_st = py/cphi-h*salpha*tphi/cphi
    delta_st = delta -px*calpha*tphi-py*salpha*tphi+h*tphi*tphi

    pz_st = np.sqrt((1.+delta_st)**2-px_st**2-py_st**2)
    hx_st = px_st/pz_st
    hy_st = py_st/pz_st
    hsigma_st = 1.-(delta_st+1)/pz_st

    L11 = 1.+hx_st*calpha*sphi
    L12 = hx_st*salpha*sphi
    L13 = calpha*tphi

    L21 = hy_st*calpha*sphi
    L22 = 1.+hy_st*salpha*sphi
    L23 = salpha*tphi

    L31 = hsigma_st*calpha*sphi
    L32 = hsigma_st*salpha*sphi
    L33 = 1./cphi

    x_st = L11*x + L12*y + L13*sigma
    y_st = L21*x + L22*y + L23*sigma
    sigma_st = L31*x + L32*y + L33*sigma
    
    return x_st, px_st, y_st, py_st, sigma_st, delta_st
    
        
    
