import metaclass
import numpy as np

ob = metaclass.twiss('twiss_s.tfs')

R = ob.SIG11-ob.SIG33
W = ob.SIG11+ob.SIG33
T = R*R+4.*ob.SIG13*ob.SIG13

cos2theta = -np.sign(R)*R/np.sqrt(T)
