import numpy as np
import sys, os
import pylab as pl

sys.path.append('../')

import mystyle as ms
import propagate_sigma_matrix as psm


# First test
SIG11 = 10
SIG33 = 10
SIG13 = 0. 
SIG12 = .5
SIG22 = 0.2
SIG14 = 0.5
SIG23 = 0
SIG24 = 0.1
SIG34 = 0.
SIG44 = 0.25

#~ # To be further tested... It makes a step
#~ SIG11 = 10
#~ SIG33 = 10
#~ SIG13 = 0. 
#~ SIG12 = .5
#~ SIG22 = 0.2
#~ SIG14 = 0.5
#~ SIG23 = 0
#~ SIG24 = 0.1
#~ SIG34 = 0.5
#~ SIG44 = 0.25


# Generate S vector for test
S = np.linspace(-2, 2, 201)

# Propagate Sigma matrix
Sigmas_at_0 = psm.Sigmas(SIG11, SIG12, SIG13, SIG14,
                     SIG22, SIG23, SIG24, SIG33,
                     SIG34, SIG44)

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

sqrtT = np.sqrt(T)
signR = np.sign(R)

cos2theta = signR*R/sqrtT
costheta = np.sqrt(0.5*(1.+cos2theta))
sintheta = signR*np.sign(Sig_13)*np.sqrt(0.5*(1.-cos2theta))


dS_R = 2.*(Sig_12_0-Sig_34_0)+2*S*(Sig_22_0-Sig_44_0)
dS_W = 2.*(Sig_12_0+Sig_34_0)+2*S*(Sig_22_0+Sig_44_0)
dS_Sig_13 = Sig_14_0 + Sig_23_0 + 2*Sig_24_0*S
dS_T = 2*R*dS_R+8.*Sig_13*dS_Sig_13
dS_cos2theta = signR*(dS_R/sqrtT - R/(2*sqrtT*sqrtT*sqrtT)*dS_T)
dS_costheta = 1/(4*costheta)*dS_cos2theta
dS_sintheta = -1/(4*sintheta)*dS_cos2theta
dS_Sig_11_hat = 0.5*(dS_W + signR*0.5/sqrtT*dS_T)
dS_Sig_33_hat = 0.5*(dS_W - signR*0.5/sqrtT*dS_T)






# Singular evaluation

a = Sig_12_0-Sig_34_0
b = Sig_22_0-Sig_44_0
c = Sig_14_0+Sig_23_0
d = Sig_24_0

mysign = lambda u: 2*np.float_(u>=0)-1.

T_test = 4*S**2*((a**2+c**2)+S*(a*b+2*c*d))
cos2theta_test = np.abs(2*a+b*S)/(2*np.sqrt(a**2+c**2+S*(a*b+2*c*d)))
#cos2theta_test3 = np.abs(2*a+b*S)/(2*np.sqrt(a**2+c**2))*(1.-S/2*(a*b+2*c*d)/(a**2+c**2))
cos2theta_test2 = mysign(2*a+b*S)/(2*np.sqrt(a**2+c**2))*(2*a+S*(b-a*(a*b+2*c*d)/(a**2+c**2)))

dS_cos2theta_test = 0.5*mysign(2*a+b*S)*(b*(a**2+c**2+S*(a*b+2*c*d))**(-0.5)\
                -0.5*(2*a+b*S)*(a**2+c**2+S*(a*b+2*c*d))**(-1.5)*(a*b+2*c*d))
dS_cos2theta_test2 = mysign(2*a+b*S)/(2*np.sqrt(a**2+c**2))*(b-a*(a*b+2*c*d)/(a**2+c**2))
cos2theta_S0 = np.interp(0, S, cos2theta)


sintheta_test = mysign((2*a*S+b*S**2)*(c*S+d*S**2))*np.sqrt(0.5*(1.-cos2theta_test))
costheta_test = np.sqrt(0.5*(1.+cos2theta_test))

dS_costheta_test = 1/(4*costheta_test)*dS_cos2theta_test
dS_sintheta_test = -1/(4*sintheta_test)*dS_cos2theta_test

dS_Sig_11_hat_test = 0.5*dS_W + mysign(a*S+b*S**2)*mysign(S)*(np.sqrt(a**2+c**2+S*(a*b+2*c*d))+0.5*S*(a*b+2*c*d)/np.sqrt(a**2+c**2+S*(a*b+2*c*d)))



pl.close('all')
pl.figure(1)
sp0 = pl.subplot(1,1,1)
pl.plot(S, cos2theta, '.-', label = 'cos2theta')
pl.plot(S, cos2theta_test, 'r.-', label = 'cos2theta_test')
pl.plot(S, cos2theta_test2, 'g.-', label = 'cos2theta_test2')
pl.plot(S, cos2theta_S0+dS_cos2theta_test*S, 'k', label = 'derivative')
pl.legend(loc='best')

pl.figure(2)
pl.subplot(1,1,1, sharex = sp0)
pl.plot(S, T, '.-', label = 'T')
pl.plot(S, T_test, 'r.-', label = 'T test')
pl.legend(loc='best')
pl.suptitle('a=%.2f, b=%.2f, c=%.2f, d=%.2f'%(a,b,c,d))

pl.figure(3)
pl.subplot(2,1,1, sharex = sp0)
pl.plot(S, sintheta, '.-', label = 'sintheta')
pl.plot(S, sintheta_test, 'r.-', label = 'sintheta test')
pl.legend(loc='best')
pl.subplot(2,1,2, sharex = sp0)
pl.plot(S, costheta, '.-', label = 'costheta')
pl.plot(S, costheta_test, 'r.-', label = 'costheta test')
pl.legend(loc='best')
pl.suptitle('a=%.2f, b=%.2f, c=%.2f, d=%.2f'%(a,b,c,d))

pl.figure(4)
pl.subplot(2,1,1, sharex = sp0)
pl.plot(S, dS_sintheta, '.-', label = 'dS_sintheta')
pl.plot(S, dS_sintheta_test, 'r.-', label = 'dS_sintheta test')
pl.legend(loc='best')
pl.subplot(2,1,2, sharex = sp0)
pl.plot(S, dS_costheta, '.-', label = 'dS_costheta')
pl.plot(S, dS_costheta_test, 'r.-', label = 'dS_costheta test')
pl.legend(loc='best')
pl.suptitle('a=%.2f, b=%.2f, c=%.2f, d=%.2f'%(a,b,c,d))

pl.figure(5)
pl.subplot(2,1,1, sharex = sp0)
pl.plot(S, dS_Sig_11_hat, '.-', label = 'dS_Sig_11_hat')
pl.plot(S, dS_Sig_11_hat_test, 'r.-', label = 'dS_Sig_11_hat_test')
pl.legend(loc='best')
pl.suptitle('a=%.2f, b=%.2f, c=%.2f, d=%.2f'%(a,b,c,d))


pl.legend(loc='best')



pl.show()



