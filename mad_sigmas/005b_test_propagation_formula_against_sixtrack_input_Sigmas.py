import numpy as np
import sys, os
import pylab as pl

sys.path.append('../')

import mystyle as ms
import propagate_sigma_matrix as psm


# Case T=0., c = 0., |a|>0
SIG11 = 10.
SIG33 = 10.
SIG13 = 0. 
SIG12 = -.5
SIG22 = 0.2
SIG14 = 0.5
SIG23 = -0.5
SIG24 = -0.1
SIG34 = 0.15
SIG44 = 0.25

# Case T=0., c = 0., |a|>0, d = 0 decoupled case
SIG11 = 10
SIG33 = 10
SIG13 = 0. 
SIG12 = -.5
SIG22 = 0.2
SIG14 = 0.5
SIG23 = -0.5
SIG24 = 0.
SIG34 = 0.
SIG44 = 0.25

# Case T=0., c = 0., a = 0.
SIG11 = 10.
SIG33 = 10.
SIG13 = 0. 
SIG12 = -.5
SIG22 = 0.2
SIG14 = 0.5
SIG23 = -0.5
SIG24 = 0.1
SIG34 = -.5
SIG44 = 0.25

#decoupled case
SIG11 = 10.
SIG33 = 15.
SIG13 = 0. 
SIG12 = 0.
SIG22 = 0.2
SIG14 = 0.
SIG23 = 0.
SIG24 = 0.
SIG34 = -.5
SIG44 = 0.25

# Case |T|>0., |S31|>0
SIG11 = 10.
SIG33 = 12.
SIG13 = 1.
SIG12 = -.5
SIG22 = 0.2
SIG14 = -0.5
SIG23 = 0.9
SIG24 = 0.1
SIG34 = -0.7
SIG44 = 0.23

#~ # Case |T|>0., S31=0
#~ SIG11 = 10
#~ SIG33 = 12
#~ SIG13 = 0. 
#~ SIG12 = -.5
#~ SIG22 = 0.2
#~ SIG14 = -0.5
#~ SIG23 = 0.2
#~ SIG24 = 0.1
#~ SIG34 = 0.
#~ SIG44 = 0.2

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
    extra_data = psm.propagate_Sigma_matrix_vectorized(Sigmas_at_0, S)

# Extract extra data
Sig_11 = extra_data['Sig_11']
Sig_33 = extra_data['Sig_33']
Sig_13 = extra_data['Sig_13']



#check against sixtrack
import sigmas_sixtrak as sigs
bcu = np.array([\
Sigmas_at_0.Sig_11_0,
Sigmas_at_0.Sig_33_0,
Sigmas_at_0.Sig_13_0,
Sigmas_at_0.Sig_12_0,
Sigmas_at_0.Sig_14_0,
Sigmas_at_0.Sig_22_0,
Sigmas_at_0.Sig_23_0,
Sigmas_at_0.Sig_24_0,
Sigmas_at_0.Sig_34_0,
Sigmas_at_0.Sig_44_0])



propagate_sigma_sixtrack_vectorized = np.vectorize(sigs.propagate_sigma, excluded=['bcu'])

sinth_sixt, costh_sixt, sx_sixt, sy_sixt, costhp_sixt, sinthp_sixt,\
 ds_sx_sixt, ds_sy_sixt, Sig_11_sixt, Sig_33_sixt, Sig_13_sixt = propagate_sigma_sixtrack_vectorized(sp=S,bcu=bcu)

# Plot results of the tests
pl.close('all')
fontsz = 14
lw = 3
ms.mystyle(fontsz = fontsz)




fig2 = pl.figure(2)
fig2.set_facecolor('w')
sp0 = pl.subplot(2,1,1)
pl.plot(S, Sig_11_hat, 'b', label = 'Sig_11_hat', lw=lw)
pl.plot(S, sx_sixt, 'c--', label = 'Sig_11_hat', lw=lw)
pl.plot(S, Sig_33_hat, 'r', label = 'Sig_33_hat', lw=lw)
pl.plot(S, sy_sixt, 'm--', label = 'Sig_33_hat', lw=lw)
ms.sciy()
pl.legend(loc='best', prop={'size':fontsz})
pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, costheta, 'b', label='costheta', lw=lw)
pl.plot(S, costh_sixt, 'c--', label='costheta_sixt', lw=lw)
pl.plot(S, sintheta, 'r', label='sintheta', lw=lw)
pl.plot(S, sinth_sixt, 'm--', label='sintheta_sixt', lw=lw)
pl.legend(loc='best', prop={'size':fontsz})

fig2.subplots_adjust(top=.82)

fig20 = pl.figure(20)
fig20.set_facecolor('w')
sp0 = pl.subplot(2,1,1)
pl.plot(S, Sig_11, 'b', label = 'Sig_11', lw=lw)
pl.plot(S, Sig_11_sixt, 'c--', label = 'Sig_11_sixt', lw=lw)
pl.plot(S, Sig_33, 'r', label = 'Sig_33', lw=lw)
pl.plot(S, Sig_33_sixt, 'm--', label = 'Sig_33_sixt', lw=lw)
ms.sciy()
pl.legend(loc='best', prop={'size':fontsz})
pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, Sig_13, 'b', label = 'Sig_13', lw=lw)
pl.plot(S, Sig_13_sixt, 'c--', label = 'Sig_13_sixt', lw=lw)
ms.sciy()


fig30 = pl.figure(30)
fig30.set_facecolor('w')
sp0 = pl.subplot(2,1,1)
pl.plot(S, dS_Sig_11_hat, 'b', label = 'dS_Sig_11_hat', lw=lw)
pl.plot(S, ds_sx_sixt, 'c--', label = 'dS_Sig_11_hat_sixt', lw=lw)
pl.plot(S, dS_Sig_33_hat, 'r', label = 'Sig_33_hat', lw=lw)
pl.plot(S, ds_sy_sixt, 'm--', label = 'dS_Sig_33_hat_sixt', lw=lw)
ms.sciy()
pl.legend(loc='best', prop={'size':fontsz})
pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, dS_costheta, 'b', label='dS_costheta', lw=lw)
pl.plot(S, costhp_sixt, 'c--', label='dS_costheta_sixt', lw=lw)
pl.plot(S, dS_sintheta, 'r', label='dS_sintheta', lw=lw)
pl.plot(S, sinthp_sixt, 'm--', label='dS_sintheta_sixt', lw=lw)
pl.legend(loc='upper left', prop={'size':fontsz})


fig3 = pl.figure(3)
fig3.set_facecolor('w')
theta = np.arctan2(sintheta, costheta)
theta_sixt = np.arctan2(sinth_sixt, costh_sixt)
pl.plot(S, theta*180/np.pi, lw=lw, label='Theta')
pl.plot(S, theta_sixt*180/np.pi, 'c--', lw=lw, label='Theta_sixt')
pl.legend(loc='upper left', prop={'size':fontsz})






pl.show()
