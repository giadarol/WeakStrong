import numpy as np
import os
import pylab as pl
import mystyle as ms


L_line = 100
Ds = 2.
skew_at = 10.
L_skew = 1.

inserted_skew = False

with open('madauto.madx', 'w') as fid:
    fid .write('''
q1: quadrupole, l=%.2f, k1=0.02, tilt=0.5;
s: sequence, l=%e;
'''%(L_skew, L_line))
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
twiss,betx=1500,bety=3000,alfx=100, alfy=-50, ripken,file=twiss_s.tfs;
''')

os.system('madx madauto.madx')

import metaclass
ob = metaclass.twiss('twiss_s.tfs')

s_ref_close_to = 60.
i_ref = np.argmin(np.abs(ob.S-s_ref_close_to))

s_ref = ob.S[i_ref]
S = np.linspace(skew_at+L_skew+2., max(ob.S), 1000)-s_ref


Sig_11_0 = ob.SIG11[i_ref]
Sig_12_0 = ob.SIG12[i_ref]
Sig_13_0 = ob.SIG13[i_ref]
Sig_14_0 = ob.SIG14[i_ref]
Sig_22_0 = ob.SIG22[i_ref]
Sig_23_0 = ob.SIG23[i_ref]
Sig_24_0 = ob.SIG24[i_ref]
Sig_33_0 = ob.SIG33[i_ref]
Sig_34_0 = ob.SIG34[i_ref]
Sig_44_0 = ob.SIG44[i_ref]

Sig_11 = Sig_11_0 + 2.*Sig_12_0*S+Sig_22_0*S*S
Sig_33 = Sig_33_0 + 2.*Sig_34_0*S+Sig_44_0*S*S
Sig_13 = Sig_13_0 + (Sig_14_0+Sig_23_0)*S+Sig_24_0*S*S

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

#evaluate derivatives
dS_R = 2.*(Sig_12_0-Sig_34_0)+2*S*(Sig_22_0-Sig_44_0)
dS_W = 2.*(Sig_12_0+Sig_34_0)+2*S*(Sig_22_0+Sig_44_0)
dS_Sig_13 = Sig_14_0 + Sig_23_0 + 2*Sig_24_0*S
dS_T = 2*R*dS_R+8.*Sig_13*dS_Sig_13


dS_cos2theta = -signR*(dS_R/sqrtT - R/(2*sqrtT*sqrtT*sqrtT)*dS_T)
dS_costheta = 1/(4*costheta)*dS_cos2theta
dS_sintheta = -1/(4*sintheta)*dS_cos2theta

dS_Sig_11_hat = 0.5*(dS_W + signR*0.5/sqrtT*dS_T)
dS_Sig_33_hat = 0.5*(dS_W - signR*0.5/sqrtT*dS_T)



pl.close('all')
pl.figure(1)
sp0 = pl.subplot(3,1,1)
pl.plot(S, Sig_11)
pl.plot(ob.S-s_ref, ob.SIG11, '.r')
ms.sciy()

pl.subplot(3,1,2, sharex=sp0)
pl.plot(S, Sig_33)
pl.plot(ob.S-s_ref, ob.SIG33, '.r')
ms.sciy()

pl.subplot(3,1,3, sharex=sp0)
pl.plot(S, Sig_13)
pl.plot(ob.S-s_ref, ob.SIG13, '.r')
ms.sciy()
pl.suptitle('Check optics propagation against MAD-X')


#Check with linear algebra
lam_1 = []
lam_2 = []
cos_1 = []
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
    v2x = v[0,1]
    v2y = v[1,1]
    cos_1.append(v1x/np.sqrt(v1x**2+v1y**2))
    cos_2.append(v2x/np.sqrt(v2x**2+v2y**2))
    
    sin_2.append(v2y/np.sqrt(v2x**2+v2y**2))
    
pl.figure(2)
pl.subplot(2,1,1, sharex=sp0)
pl.plot(S, lam_1)
pl.plot(S, lam_2)
pl.plot(S, Sig_11_hat)
pl.plot(S, Sig_33_hat)
pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, costheta, 'b')
pl.plot(S, sintheta, 'r')
#pl.plot(S, cos_1)
pl.plot(S, cos_2, 'c--')
pl.plot(S, sin_2, 'm--')
pl.suptitle('Check rotation against matrix diagonalization')

pl.figure(3)
pl.subplot(2,1,1, sharex=sp0)
pl.plot(S, dS_costheta, 'm-')
pl.plot(S[:-1], np.diff(costheta)/np.diff(S), 'r--')
pl.plot(S, dS_sintheta, 'b-')
pl.plot(S[:-1], np.diff(sintheta)/np.diff(S), 'c--')
pl.ylim(-.02, .02)

pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, dS_Sig_11_hat, 'm-')
pl.plot(S[:-1], np.diff(Sig_11_hat)/np.diff(S), 'r--')
pl.plot(S, dS_Sig_33_hat, 'b-')
pl.plot(S[:-1], np.diff(Sig_33_hat)/np.diff(S), 'c--')
pl.suptitle('Check derivatives against finite differeces')
pl.ylim(-1e-6, 1e-6)



pl.show()
