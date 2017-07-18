import numpy as np
import sys, os
import pylab as pl

sys.path.append('../')

import mystyle as ms
import propagate_sigma_matrix as psm


#~ # Case T=0., |c|>0.
#~ SIG11 = 10
#~ SIG33 = 10
#~ SIG13 = 0. 
#~ 
#~ SIG12 = -.5
#~ SIG22 = 0.2
#~ SIG14 = 0.5
#~ SIG23 = 0
#~ SIG24 = 0.1
#~ SIG34 = 0.
#~ SIG44 = 0.2

# Case T=0., c = 0., |a|>0
SIG11 = 10
SIG33 = 10
SIG13 = 0. 

SIG12 = -.5
SIG22 = 0.2
SIG14 = 0.5
SIG23 = -0.5
SIG24 = 0.1
SIG34 = 0.
SIG44 = 0.25



# Propagate by:
DS = -4
SIG11_DS, SIG12_DS, SIG13_DS, SIG14_DS,\
    SIG22_DS, SIG23_DS, SIG24_DS, SIG33_DS,\
    SIG34_DS, SIG44_DS = psm.propagate_full_Sigma_matrix_in_drift(\
                     SIG11, SIG12, SIG13, SIG14,
                     SIG22, SIG23, SIG24, SIG33,
                     SIG34, SIG44, DS)



# Generate S vector for test
S = np.linspace(-5-DS, 5-DS, 21)
#~ S = np.linspace(-5-DS, 5-DS, 201)


# Propagate Sigma matrix
Sigmas_at_0 = psm.Sigmas(SIG11_DS, SIG12_DS, SIG13_DS, SIG14_DS,
                     SIG22_DS, SIG23_DS, SIG24_DS, SIG33_DS,
                     SIG34_DS, SIG44_DS)
                                          
Sig_11_hat, Sig_33_hat, costheta, sintheta, \
    dS_Sig_11_hat, dS_Sig_33_hat, dS_costheta, dS_sintheta,\
    extra_data = psm.propagate_Sigma_matrix_vectorized(Sigmas_at_0, S, handle_singularities=False)

Sig_11_hat_s, Sig_33_hat_s, costheta_s, sintheta_s, \
    dS_Sig_11_hat_s, dS_Sig_33_hat_s, dS_costheta_s, dS_sintheta_s,\
    extra_data_s = psm.propagate_Sigma_matrix_vectorized(Sigmas_at_0, S, handle_singularities=True)
    
#~ # Extract extra data
cos2theta = []
T = []
R = []
for ii in xrange(len(S)):
    cos2theta.append(extra_data_s[ii]['cos2theta'])
    T.append(extra_data_s[ii]['T'])
    R.append(extra_data_s[ii]['R'])
    
T = np.array(T)
R = np.array(R)

# Plot results of the tests
pl.close('all')
fontsz = 14
lw = 3
mks = 10
ms.mystyle(fontsz = fontsz)


pl.close('all')

# some investigations
a = SIG12-SIG34
b = SIG22-SIG44
c = SIG14+SIG23
d = SIG24
ddSS = S + DS

pl.figure(1001)
pl.plot(S, T)
T1 = 4*ddSS**2*(a**2+c**2+ddSS*(a*b+2*c*d))
#~ pl.plot(S, ddSS**2*(2*a+b*ddSS)**2 +4*ddSS**2*(c+d*ddSS)**2)
#~ pl.plot(S, ddSS**2*(4*a**2+4*a*b*ddSS) +4*ddSS**2*(c**2+2*c*d*ddSS))
pl.plot(S, T1)

pl.figure(1002)
pl.plot(S, R)
pl.plot(S, 2*a*ddSS+b*ddSS**2)




pl.figure(1000)
pl.plot(S, cos2theta)
pl.plot(S, np.abs(2*a+b*ddSS)/(np.sqrt((2*a+b*ddSS)**2 +4*(c+d*ddSS)**2)))
print 'It seems that the following expansion is not appropriate!'
#~ pl.plot(S,np.abs(2*a+b*ddSS)/(2*np.sqrt(a**2+c**2+ddSS*(a*b+2*c*d))))


fig2 = pl.figure(2); pl.clf()
fig2.set_facecolor('w')
sp0 = pl.subplot(2,1,1)
pl.plot(S, Sig_11_hat, '.-b', label = 'Sig_11_hat', lw=lw, markersize=mks)
pl.plot(S, Sig_33_hat, '.-r', label = 'Sig_33_hat', lw=lw, markersize=mks)
pl.plot(S, Sig_11_hat_s, '.b', lw=lw, markersize=mks)
pl.plot(S, Sig_33_hat_s, '.r', lw=lw, markersize=mks)
ms.sciy()
pl.grid('on')
pl.legend(loc='best', prop={'size':fontsz})
pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, costheta, '.-b', label='costheta', lw=lw, markersize=mks)
pl.plot(S, sintheta, '.-r', label='sintheta', lw=lw, markersize=mks)
pl.plot(S, costheta_s, '.--b', lw=lw/2, markersize=mks)
pl.plot(S, sintheta_s, '.--r', lw=lw/2, markersize=mks)
pl.legend(loc='best', prop={'size':fontsz})
pl.grid('on')
fig2.subplots_adjust(top=.82)

fig3 = pl.figure(3); pl.clf()
fig3.set_facecolor('w')
pl.subplot(2,1,1, sharex=sp0)
pl.plot(S, dS_Sig_11_hat, '.-b', label = 'dS_Sig_11_hat', lw=lw, markersize=mks)
pl.plot(S, dS_Sig_33_hat, '.-r', label = 'dS_Sig_33_hat', lw=lw, markersize=mks)
pl.plot(S, dS_Sig_11_hat_s, '.--b', lw=lw/2, markersize=mks)
pl.plot(S, dS_Sig_33_hat_s, '.--r', lw=lw/2, markersize=mks)
ms.sciy()
pl.grid('on')
pl.legend(loc='best', prop={'size':fontsz})
pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, dS_costheta, '.-b', label='dS_costheta', lw=lw, markersize=mks)
pl.plot(S, dS_sintheta, '.-r', label='dS_sintheta', lw=lw, markersize=mks)
pl.plot(S, dS_costheta_s, '.--b', lw=lw/2, markersize=mks)
pl.plot(S, dS_sintheta_s, '.--r', lw=lw/2, markersize=mks)
pl.legend(loc='best', prop={'size':fontsz})
pl.grid('on')
fig2.subplots_adjust(top=.82)


pl.show()



