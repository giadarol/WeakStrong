
# I program as close as possible to C...

class parboost(object):
    #it is practically a struct
    def __init__(self, phi, alpha):
        self.sphi = np.sin(phi)
        self.cphi = np.cos(phi)
        self.tphi = np.tan(phi)
        self.salpha = np.sin(alpha)
        self.calpha = np.cos(alpha)

def boost(x, px, y, py, sigma, delta, parboost):
    
        
    
