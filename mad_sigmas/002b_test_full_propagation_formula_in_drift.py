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
Sig_11, Sig_12, Sig_13, Sig_14, Sig_22, Sig_23, Sig_24,\
    Sig_33, Sig_34, Sig_44 = psm.propagate_full_Sigma_matrix_in_drift(ob.SIG11[i_ref], 
                     ob.SIG12[i_ref], ob.SIG13[i_ref], ob.SIG14[i_ref],
                     ob.SIG22[i_ref], ob.SIG23[i_ref], ob.SIG24[i_ref], ob.SIG33[i_ref],
                     ob.SIG34[i_ref], ob.SIG44[i_ref], S)
                                          




# Plot results of the tests
pl.close('all')
fontsz = 14
lw = 3
ms.mystyle(fontsz = fontsz)

pl.figure(1); pl.clf()
sp0 = pl.subplot(2,1,1)
pl.plot(S, Sig_11, 'b', label='Sig_11')
pl.plot(ob.S-s_ref, ob.SIG11, '.c')
pl.plot(S, Sig_33, 'r', label='Sig_33')
pl.plot(ob.S-s_ref, ob.SIG33, '.m')
pl.legend(loc='best', prop = {'size':fontsz})
ms.sciy()

pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, Sig_13, label='Sig_13')
pl.plot(ob.S-s_ref, ob.SIG13, '.r')
pl.legend(loc='best', prop = {'size':fontsz})
ms.sciy()
pl.suptitle('Check optics propagation against MAD-X')

pl.figure(2); pl.clf()
sp0 = pl.subplot(2,1,1)
pl.plot(S, Sig_12, 'b', label='Sig_12')
pl.plot(ob.S-s_ref, ob.SIG12, '.b')
pl.plot(S, Sig_14, 'r', label='Sig_14')
pl.plot(ob.S-s_ref, ob.SIG14, '.m')
pl.plot(S, Sig_23, 'g', label='Sig_23')
pl.plot(ob.S-s_ref, ob.SIG23, '.g')
pl.plot(S, Sig_34, 'c', label='Sig_34')
pl.plot(ob.S-s_ref, ob.SIG34, '.c')

pl.legend(loc='best', prop = {'size':fontsz})
ms.sciy()

pl.subplot(2,1,2, sharex=sp0)
pl.plot(S, Sig_22, label='Sig_22')
pl.plot(ob.S-s_ref, ob.SIG22, '.b')
pl.plot(S, Sig_24,'r', label='Sig_24')
pl.plot(ob.S-s_ref, ob.SIG24, '.r')
pl.plot(S, Sig_44,'g', label='Sig_44')
pl.plot(ob.S-s_ref, ob.SIG44, '.g')
pl.legend(loc='best', prop = {'size':fontsz})
ms.sciy()
pl.suptitle('Check optics propagation against MAD-X')





pl.show()



