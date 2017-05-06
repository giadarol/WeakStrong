import numpy as np
import pylab as pl
import mystyle as ms

pl.close('all')

Sig_11 = 10.
Sig_33 = 5.
Sig_13 = 0.1


R = Sig_11-Sig_33
W = Sig_11+Sig_33
T = R*R+4*Sig_13*Sig_13

sqrtT = np.sqrt(T)
signR = np.sign(R)

cos2theta = -signR*R/sqrtT
costheta = np.sqrt(0.5*(1.+cos2theta))
sintheta = -np.sign((Sig_11-Sig_33)*Sig_13)*np.sqrt(0.5*(1.-cos2theta))

# in sixtrack this line seems to be different different
# sintheta = -np.sign((Sig_11-Sig_33))*np.sqrt(0.5*(1.-cos2theta))

Sig_11_hat = 0.5*(W+signR*sqrtT)
Sig_33_hat = 0.5*(W-signR*sqrtT)


a = np.array([[Sig_11, Sig_13],
              [Sig_13, Sig_33]])
w, v = np.linalg.eig(a)

theta = np.linspace(0, 2*np.pi, 100)



for tt in theta:
    res = np.dot(a, np.array([np.cos(tt), np.sin(tt)]).T)
    pl.plot(res[0], res[1], '.')
pl.axis('equal')

pl.show()
