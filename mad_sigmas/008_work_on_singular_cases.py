import numpy as np
import sys, os
import pylab as pl

sys.path.append('../')

import mystyle as ms
import propagate_sigma_matrix as psm


# Case T=0., |c|>0.
SIG11 = 10
SIG33 = 10
SIG13 = 0. 

SIG12 = .5
SIG22 = 0.2
SIG14 = 0.5
SIG23 = 0
SIG24 = 0.1
SIG34 = 0.
SIG44 = 0.2


# Generate S vector for test
S = np.linspace(-5, 5, 21)

# Propagate Sigma matrix
Sigmas_at_0 = psm.Sigmas(SIG11, SIG12, SIG13, SIG14,
                     SIG22, SIG23, SIG24, SIG33,
                     SIG34, SIG44)
                                          
Sig_11_hat, Sig_33_hat, costheta, sintheta, \
    dS_Sig_11_hat, dS_Sig_33_hat, dS_costheta, dS_sintheta,\
    extra_data = psm.propagate_Sigma_matrix_vectorized(Sigmas_at_0, S)
    
# Extract extra data
Sig_11 = extra_data['Sig_11']
Sig_33 = extra_data['Sig_33']
Sig_13 = extra_data['Sig_13']


# Plot results of the tests
pl.close('all')
fontsz = 14
lw = 3
mks = 10
ms.mystyle(fontsz = fontsz)

pl.figure(1); pl.clf()
sp0 = pl.subplot(2,1,1)
pl.plot(S, Sig_11)
pl.plot(S, Sig_33, 'r')


ms.sciy()

pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, Sig_13)

ms.sciy()
pl.suptitle('Check optics propagation against MAD-X')



#dS_Sig_11_hat, dS_Sig_33_hat, dS_costheta, dS_sintheta,\


pl.close('all')

fig2 = pl.figure(2); pl.clf()
fig2.set_facecolor('w')
pl.subplot(2,1,1, sharex=sp0)
pl.plot(S, Sig_11_hat, '.-b', label = 'Sig_11_hat', lw=lw, markersize=mks)
pl.plot(S, Sig_33_hat, '.-r', label = 'Sig_33_hat', lw=lw, markersize=mks)
ms.sciy()
pl.grid('on')
pl.legend(loc='best', prop={'size':fontsz})
pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, costheta, '.-b', label='costheta', lw=lw, markersize=mks)
pl.plot(S, sintheta, '.-r', label='sintheta', lw=lw, markersize=mks)
pl.legend(loc='best', prop={'size':fontsz})
pl.grid('on')
fig2.subplots_adjust(top=.82)

fig3 = pl.figure(3); pl.clf()
fig3.set_facecolor('w')
pl.subplot(2,1,1, sharex=sp0)
pl.plot(S, dS_Sig_11_hat, '.-b', label = 'dS_Sig_11_hat', lw=lw, markersize=mks)
pl.plot(S, dS_Sig_33_hat, '.-r', label = 'dS_Sig_33_hat', lw=lw, markersize=mks)
ms.sciy()
pl.grid('on')
pl.legend(loc='best', prop={'size':fontsz})
pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, dS_costheta, '.-b', label='dS_costheta', lw=lw, markersize=mks)
pl.plot(S, dS_sintheta, '.-r', label='dS_sintheta', lw=lw, markersize=mks)
pl.legend(loc='best', prop={'size':fontsz})
pl.grid('on')
fig2.subplots_adjust(top=.82)


pl.show()



