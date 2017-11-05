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
#~ k_skew = .02
#~ betx_enter = 1500
#~ bety_enter = 3000
#~ alfx_enter = -100
#~ alfy_enter = 50
#~ 
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
#~ 
#~ # To have pass through a round beam
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
                pos_skew = ss+L_skew+0.1
                fid .write('q1, at=%.2f;\n'%(pos_skew))
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
S = np.linspace(pos_skew, max(ob.S), 1000)-s_ref

# Propagate Sigma matrix
Sig_11, Sig_12, Sig_13, Sig_14, Sig_22, Sig_23, Sig_24,\
    Sig_33, Sig_34, Sig_44 = psm.propagate_full_Sigma_matrix_in_drift(ob.SIG11[i_ref], 
                     ob.SIG12[i_ref], ob.SIG13[i_ref], ob.SIG14[i_ref],
                     ob.SIG22[i_ref], ob.SIG23[i_ref], ob.SIG24[i_ref], ob.SIG33[i_ref],
                     ob.SIG34[i_ref], ob.SIG44[i_ref], S)
                                          


# Propagate Sigma matrix and diagonalize
Sigmas_at_0 = psm.Sigmas(ob.SIG11[i_ref], ob.SIG12[i_ref], ob.SIG13[i_ref], ob.SIG14[i_ref],
                     ob.SIG22[i_ref], ob.SIG23[i_ref], ob.SIG24[i_ref], ob.SIG33[i_ref],
                     ob.SIG34[i_ref], ob.SIG44[i_ref])
Sig_11_hat, Sig_33_hat, costheta, sintheta, \
    dS_Sig_11_hat, dS_Sig_33_hat, dS_costheta, dS_sintheta,\
    extra_data = psm.propagate_Sigma_matrix_vectorized(Sigmas_at_0, S, handle_singularities=True, threshold_singular = 1e-30)


# Plot results of the tests
pl.close('all')
fontsz = 14
lw = 3
mks = 8
ms.mystyle_arial(fontsz = fontsz, dist_tick_lab=5)

fig1 = pl.figure(1); pl.clf()
fig1.set_facecolor('w')
sp0 = pl.subplot(2,1,1)
pl.plot(S, Sig_11, 'b', label='Sig_11', lw=lw)
pl.plot(ob.S-s_ref, ob.SIG11, '.c', markersize=mks, lw=lw)#, label='Sig_11 (MAD)'
pl.plot(S, Sig_33, 'r', label='Sig_33', lw=lw)
pl.plot(ob.S-s_ref, ob.SIG33, '.m', markersize=mks)#, label='Sig_33 (MAD)')
pl.axvline(x=pos_skew-s_ref, lw=lw, color='grey', linestyle='--')
pl.legend(loc='best', prop = {'size':fontsz})
pl.grid('on')
ms.sciy()

sp1 = pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, Sig_13, label='Sig_13', lw=lw)
pl.plot(ob.S-s_ref, ob.SIG13, '.c', markersize=mks)
pl.axvline(x=pos_skew-s_ref, lw=lw, color='grey', linestyle='--')
pl.xlabel('s [m]')
pl.legend(loc='best', prop = {'size':fontsz})
pl.grid('on')
ms.sciy()
pl.suptitle('Check optics propagation against MAD-X')



fig101 = pl.figure(101); pl.clf()
fig101.set_facecolor('w')
sp0 = pl.subplot(2,1,1)
pl.plot(S, Sig_12, 'b', label='Sig_12', lw=lw, markersize=mks)
pl.plot(ob.S-s_ref, ob.SIG12, '.b', lw=lw, markersize=mks)
pl.plot(S, Sig_14, 'r', label='Sig_14', lw=lw, markersize=mks)
pl.plot(ob.S-s_ref, ob.SIG14, '.m', lw=lw, markersize=mks)
pl.plot(S, Sig_23, 'g', label='Sig_23', lw=lw, markersize=mks)
pl.plot(ob.S-s_ref, ob.SIG23, '.g', lw=lw, markersize=mks)
pl.plot(S, Sig_34, 'c', label='Sig_34', lw=lw, markersize=mks)
pl.plot(ob.S-s_ref, ob.SIG34, '.c', lw=lw, markersize=mks)
pl.axvline(x=pos_skew-s_ref, lw=lw, color='grey', linestyle='--')
pl.grid('on')

pl.legend(loc='center right', prop = {'size':fontsz})
ms.sciy()

pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, Sig_22, label='Sig_22', lw=lw, markersize=mks)
pl.plot(ob.S-s_ref, ob.SIG22, '.b', lw=lw, markersize=mks)
pl.plot(S, Sig_24,'r', label='Sig_24', lw=lw, markersize=mks)
pl.plot(ob.S-s_ref, ob.SIG24, '.r', lw=lw, markersize=mks)
pl.plot(S, Sig_44,'g', label='Sig_44', lw=lw, markersize=mks)
pl.plot(ob.S-s_ref, ob.SIG44, '.g', lw=lw, markersize=mks)
pl.axvline(x=pos_skew-s_ref, lw=lw, color='grey', linestyle='--')
pl.xlabel('s [m]')
pl.legend(loc='center right', prop = {'size':fontsz})
ms.sciy()
pl.grid('on')
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
pl.grid('on')
ms.sciy()
pl.legend(loc='best', prop={'size':fontsz})
pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, costheta, 'b', label='costheta', lw=lw)
pl.plot(S, sintheta, 'r', label='sintheta', lw=lw)
pl.plot(S, sin_1, 'm--', label='sin (diag)', lw=lw)
pl.plot(S, cos_1, 'c--', label='cos (diag)', lw=lw)
pl.xlabel('s [m]')
pl.grid('on')
pl.legend(loc='best', prop={'size':fontsz})
pl.suptitle('Check rotation against matrix diagonalization\n'+\
            'At s=0: Sig11=%.1e,Sig22=%.1e,Sig33=%.1e,Sig44=%.1e\n'%(Sigmas_at_0.Sig_11_0, Sigmas_at_0.Sig_22_0, Sigmas_at_0.Sig_33_0, Sigmas_at_0.Sig_44_0)+\
            'Sig12=%.1e,Sig13=%.1e,Sig14=%.1e,\nSig23=%.1e,Sig24=%.1e,Sig34=%.1e'%(Sigmas_at_0.Sig_12_0, Sigmas_at_0.Sig_13_0, Sigmas_at_0.Sig_14_0, Sigmas_at_0.Sig_23_0, Sigmas_at_0.Sig_24_0, Sigmas_at_0.Sig_34_0))
theta = np.arctan2(sintheta, costheta)
fig2.subplots_adjust(top=.82)

fig3 = pl.figure(3)
fig3.set_facecolor('w')
pl.plot(S, theta*180/np.pi, lw=lw)
pl.xlabel('s [m]')
pl.ylabel('Theta [deg]')
pl.grid('on')


fig4 = pl.figure(4)
fig4.set_facecolor('w')
pl.subplot(2,1,1, sharex=sp0)
pl.plot(S, dS_costheta, 'm-', lw=lw, label = 'dS_costheta')
pl.plot(S[:-1], np.diff(costheta)/np.diff(S), 'r--', lw=lw)
pl.plot(S, dS_sintheta, 'b-', lw=lw, label = 'dS_sintheta')
pl.plot(S[:-1], np.diff(sintheta)/np.diff(S), 'c--', lw=lw)
pl.ylim(-.02, .02)
pl.grid('on')
pl.legend(loc='best', prop={'size':fontsz})
ms.sciy()

pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, dS_Sig_11_hat, 'm-', lw=lw, label = 'dS_Sig_11_hat')
pl.plot(S[:-1], np.diff(Sig_11_hat)/np.diff(S), 'r--', lw=lw)
pl.plot(S, dS_Sig_33_hat, 'b-', lw=lw, label = 'dS_Sig_33_hat')
pl.plot(S[:-1], np.diff(Sig_33_hat)/np.diff(S), 'c--', lw=lw)
pl.suptitle('Check derivatives against finite differeces')
pl.xlabel('s [m]')
pl.ylim(-1e-6, 1e-6)
pl.grid('on')
pl.legend(loc='best', prop={'size':fontsz})
ms.sciy()


for fig in [fig1, fig101, fig2, fig3, fig4]:
    fig.subplots_adjust(bottom=.12)



pl.show()



