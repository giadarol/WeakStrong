import numpy as np
import sys, os
import pylab as pl

sys.path.append('../')

import mystyle as ms
import propagate_sigma_matrix as psm


L_line = 150
Ds = 2.
skew_at = 10.
L_skew = 1.
k_skew = .02
betx_enter = 1500
bety_enter = 3000
alfx_enter = 100
alfy_enter = -50

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
    extra_data = psm.propagate_Sigma_matrix(Sigmas_at_0, S)
    
# Extract extra data
Sig_11 = extra_data['Sig_11']
Sig_33 = extra_data['Sig_33']
Sig_13 = extra_data['Sig_13']

theta = np.arctan2(sintheta, costheta)


# Plot results of the tests
pl.close('all')
fontsz = 14
lw = 3
ms.mystyle(fontsz = fontsz)

S_test_vect = np.arange(np.min(S), np.max(S), 5.)

folder_pngs = 'ellip_pngs'
if not os.path.exists(folder_pngs):
    os.makedirs(folder_pngs)

fig2 = pl.figure(2, figsize=(8*1.8, 6))
fig2.set_facecolor('w')

for i_png, S_test in enumerate(S_test_vect):
    
    fig2.clf()

    sp0 = pl.subplot(2,2,1)
    pl.plot(S, Sig_11_hat, 'b', label = 'Sig_11_hat', lw=lw)
    pl.plot(S, Sig_33_hat, 'r', label = 'Sig_33_hat', lw=lw)
    pl.ylabel('[m]')
    pl.xlabel('s [m]')
    pl.axvline(S_test, linestyle='--', lw=lw, color='k', alpha=.5)
    pl.grid('on')

    ms.sciy()
    pl.legend(loc='best', prop={'size':fontsz})
    pl.subplot(2,2,3, sharex=sp0)
    #~ pl.plot(S, costheta, 'b', label='costheta', lw=lw)
    #~ pl.plot(S, sintheta, 'r', label='sintheta', lw=lw)
    pl.plot(S, theta*180/np.pi, lw=lw, color='g')
    pl.xlabel('[m]')
    pl.ylabel('Theta [deg]')
    pl.axvline(S_test, linestyle='--', lw=lw, color='k', alpha=.5)
    pl.grid('on')
    pl.legend(loc='best', prop={'size':fontsz})


    pl.subplot2grid(shape=(2,2), loc=(0,1), rowspan=2)

    Sig_11_ellip = np.interp(S_test, S, Sig_11)
    Sig_33_ellip = np.interp(S_test, S, Sig_33)
    Sig_13_ellip = np.interp(S_test, S, Sig_13)
    Sig_11_hat_ellip = np.interp(S_test, S, Sig_11_hat)
    Sig_33_hat_ellip = np.interp(S_test, S, Sig_33_hat)
    sintheta_ellip = np.interp(S_test, S, sintheta)
    costheta_ellip = np.interp(S_test, S, costheta)
    
    scale = 1/(0.5*(np.sqrt(Sig_11_ellip)+np.sqrt(Sig_33_ellip)))

    Sigma = np.array([[Sig_11_ellip, Sig_13_ellip],
                  [Sig_13_ellip, Sig_33_ellip]])
    w, v = np.linalg.eig(Sigma)
    # This matrix transform the unitary circle into an ellipse 
    # in which the quadratic form is constant 
    # (same eigenvectors as Sigma, eigenvalues are the sqrt) 
    a = np.dot(v, np.dot(np.diag(np.sqrt(w)), v.T))

    tt_ellip = np.linspace(0, 2*np.pi, 100)

    x_ellip = []
    y_ellip = []
    for tt in tt_ellip:
        res = np.dot(a, np.array([[np.cos(tt), np.sin(tt)]]).T)
        x_ellip.append(res[0])
        y_ellip.append(res[1])
    x_ellip = np.squeeze(x_ellip)
    y_ellip = np.squeeze(y_ellip)
        
    pl.plot(scale*x_ellip, scale*y_ellip, lw=lw, color='k', alpha=.4)
    pl.plot([0, scale*np.sqrt(Sig_11_hat_ellip)*costheta_ellip], [0, scale*np.sqrt(Sig_11_hat_ellip)*sintheta_ellip],'b', lw=lw, )
    pl.plot([0, -scale*np.sqrt(Sig_33_hat_ellip)*sintheta_ellip], [0, scale*np.sqrt(Sig_33_hat_ellip)*costheta_ellip],'r', lw=lw)
    pl.xlabel('x normalized to average sigma')
    pl.ylabel('y normalized to average sigma')

    pl.grid('on')
    ms.sciy(); ms.scix()
    pl.axis('equal')
    pl.xlim(-2, 2)
    pl.ylim(-2, 2)
    pl.gca().set_aspect('equal',adjustable='box')

    
    fig2.savefig(folder_pngs+'/ellip_%02d.png'%i_png, dpi=200)


pl.show()



