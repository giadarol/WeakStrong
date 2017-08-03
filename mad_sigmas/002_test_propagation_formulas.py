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

#~ L_line = 100
#~ Ds = 2.
#~ skew_at = 10.
#~ L_skew = 1.
#~ k_skew = -.02
#~ betx_enter = 1500
#~ bety_enter = 3000
#~ alfx_enter = 100
#~ alfy_enter = -50

#~ L_line = 100
#~ Ds = 2.
#~ skew_at = 10.
#~ L_skew = 1.
#~ k_skew = .02
#~ betx_enter = 1500
#~ bety_enter = 3000
#~ alfx_enter = -100
#~ alfy_enter = 50

#~ # To have a positive Sigma_13
#~ L_line = 100
#~ Ds = 2.
#~ skew_at = 10.
#~ L_skew = 1.
#~ k_skew = -.02
#~ betx_enter = 3000
#~ bety_enter = 1500
#~ alfx_enter = -100
#~ alfy_enter = 50

# To have pass through a roung beam
#~ L_line = 100
#~ Ds = 2.
#~ skew_at = 1.
#~ L_skew = 1.
#~ k_skew = -.001
#~ betx_enter = 400
#~ bety_enter = 200
#~ alfx_enter = 5
#~ alfy_enter = 0


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

# Propagate Sigma matrix
Sigmas_at_0 = psm.Sigmas(ob.SIG11[i_ref], ob.SIG12[i_ref], ob.SIG13[i_ref], ob.SIG14[i_ref],
                     ob.SIG22[i_ref], ob.SIG23[i_ref], ob.SIG24[i_ref], ob.SIG33[i_ref],
                     ob.SIG34[i_ref], ob.SIG44[i_ref])
                                          
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
ms.mystyle(fontsz = fontsz)

pl.figure(1); pl.clf()
sp0 = pl.subplot(2,1,1)
pl.plot(S, Sig_11)
pl.plot(ob.S-s_ref, ob.SIG11, '.c')
pl.plot(S, Sig_33, 'r')
pl.plot(ob.S-s_ref, ob.SIG33, '.m')
ms.sciy()

pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, Sig_13)
pl.plot(ob.S-s_ref, ob.SIG13, '.r')
ms.sciy()
pl.suptitle('Check optics propagation against MAD-X')


#Check with linear algebra
lam_1 = []
lam_2 = []
cos_1 = []
sin_1 = []
cos_2 = []
sin_2 = []
for i_s, ss in enumerate(S):
    a = np.array([[Sig_11[i_s], Sig_13[i_s]],
                 [Sig_13[i_s], Sig_33[i_s]]])
    w, v = np.linalg.eig(a)
    lam_1.append(w[0])
    lam_2.append(w[1])
    v1x = v[0,0]
    v1y = v[1,0]
    if v1x<0: 
        v1x *= -1
        v1y *= -1
    v2x = v[0,1]
    v2y = v[1,1]
    if v2x<0: 
        v2x *= -1
        v2y *= -1
    cos_1.append(v1x/np.sqrt(v1x**2+v1y**2))
    cos_2.append(v2x/np.sqrt(v2x**2+v2y**2))
    sin_1.append(v1y/np.sqrt(v1x**2+v1y**2))    
    sin_2.append(v2y/np.sqrt(v2x**2+v2y**2))
    
fig2 = pl.figure(2); pl.clf()
fig2.set_facecolor('w')
pl.subplot(2,1,1, sharex=sp0)
pl.plot(S, Sig_11_hat, 'b', label = 'Sig_11_hat', lw=lw)
pl.plot(S, lam_1, 'c--', label= 'eigen 1', lw=lw)
pl.plot(S, Sig_33_hat, 'r', label = 'Sig_33_hat', lw=lw)
pl.plot(S, lam_2, 'm--', label= 'eigen 2', lw=lw)
ms.sciy()
pl.legend(loc='best', prop={'size':fontsz})
pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, costheta, 'b', label='costheta', lw=lw)
pl.plot(S, sintheta, 'r', label='sintheta', lw=lw)
pl.plot(S, sin_1, 'm--', label='sin (diag)', lw=lw)
pl.plot(S, cos_1, 'c--', label='cos (diag)', lw=lw)
pl.legend(loc='best', prop={'size':fontsz})
pl.suptitle('Check rotation against matrix diagonalization\n'+\
            'At s=0: Sig11=%.1e,Sig22=%.1e,Sig33=%.1e,Sig44=%.1e\n'%(Sigmas_at_0.Sig_11_0, Sigmas_at_0.Sig_22_0, Sigmas_at_0.Sig_33_0, Sigmas_at_0.Sig_44_0)+\
            'Sig12=%.1e,Sig13=%.1e,Sig14=%.1e,\nSig23=%.1e,Sig24=%.1e,Sig34=%.1e'%(Sigmas_at_0.Sig_12_0, Sigmas_at_0.Sig_13_0, Sigmas_at_0.Sig_14_0, Sigmas_at_0.Sig_23_0, Sigmas_at_0.Sig_24_0, Sigmas_at_0.Sig_34_0))
theta = np.arctan2(sintheta, costheta)
fig2.subplots_adjust(top=.82)

pl.figure(3)
pl.plot(S, theta*180/np.pi)

pl.figure(4)
pl.subplot(2,1,1, sharex=sp0)
pl.plot(S, dS_costheta, 'm-')
pl.plot(S[:-1], np.diff(costheta)/np.diff(S), 'r--')
pl.plot(S, dS_sintheta, 'b-')
pl.plot(S[:-1], np.diff(sintheta)/np.diff(S), 'c--')
pl.ylim(-.02, .02)
ms.sciy()

pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, dS_Sig_11_hat, 'm-')
pl.plot(S[:-1], np.diff(Sig_11_hat)/np.diff(S), 'r--')
pl.plot(S, dS_Sig_33_hat, 'b-')
pl.plot(S[:-1], np.diff(Sig_33_hat)/np.diff(S), 'c--')
pl.suptitle('Check derivatives against finite differeces')
pl.ylim(-1e-6, 1e-6)
ms.sciy()

pl.show()



