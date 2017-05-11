import numpy as np
import sys, os
import pylab as pl

sys.path.append('../')

import mystyle as ms
import propagate_sigma_matrix as psm


L_line = 100
Ds = 2.
skew_at = 10.
L_skew = 1.
k_skew = .02
betx_enter = 1500
bety_enter = 3000
alfx_enter = 100
alfy_enter = -50


#Build a line with strong couplig using MAD-X
inserted_skew = False
with open('madauto.madx', 'w') as fid:
    fid .write('''
q1: quadrupole, l=%.2f, k1=%.2e, tilt=0.5;
s: sequence, l=%e;
'''%(L_skew, k_skew, L_line))
    for i_s, ss in enumerate(np.arange(0., L_line, Ds)):
        fid .write('m%d: marker, at=%.2f;\n'%(i_s, ss))
        if not inserted_skew:
            if ss>skew_at:
                fid .write('q1, at=%.2f;\n'%(ss+L_skew+0.1))
                inserted_skew = True
    fid .write('endsequence;\n')

    fid .write('''
beam,ex=5e-10,ey=5e-10;
use,sequence=s;
twiss,betx=%2f,bety=%.2f,alfx=%.2f, alfy=%.2f, ripken,file=twiss_s.tfs;
'''%(betx_enter, bety_enter, alfx_enter, alfy_enter))

os.system('madx madauto.madx')
import metaclass
ob = metaclass.twiss('twiss_s.tfs')

# Choose starting point
s_ref_close_to = 60.
i_ref = np.argmin(np.abs(ob.S-s_ref_close_to))
s_ref = ob.S[i_ref]

# Generate S vector for test
S = np.linspace(skew_at+L_skew+2., max(ob.S), 1000)-s_ref
#S = np.array([0., 0.1])

# Propagate Sigma matrix
Sigmas_at_0 = psm.Sigmas(ob.SIG11[i_ref], ob.SIG12[i_ref], ob.SIG13[i_ref], ob.SIG14[i_ref],
                     ob.SIG22[i_ref], ob.SIG23[i_ref], ob.SIG24[i_ref], ob.SIG33[i_ref],
                     ob.SIG34[i_ref], ob.SIG44[i_ref])
                                          
Sig_11_hat, Sig_33_hat, costheta, sintheta, \
    dS_Sig_11_hat, dS_Sig_33_hat, dS_costheta, dS_sintheta,\
    extra_data = psm.propagate_Sigma_matrix(Sigmas_at_0, S)

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
