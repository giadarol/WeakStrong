import numpy as np
import os
import pylab as pl


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

pl.close('all')
pl.figure(1)
sp1 = pl.subplot(1,1,1)
pl.plot(ob.S, ob.BETX, 'b', label='beta_x')
pl.plot(ob.S, ob.BETY, 'r', label='beta_y')

pl.figure(2)
pl.subplot(2,2,1, sharex=sp1)
pl.plot(ob.S, ob.SIG11, 'b')
pl.plot(ob.S, ob.SIG33, 'r')
pl.plot(ob.S, ob.SIG13, 'g')
pl.grid('on')

pl.subplot(2,2,2, sharex=sp1)
pl.plot(ob.S, ob.SIG22, 'b')
pl.plot(ob.S, ob.SIG44, 'r')
pl.plot(ob.S, ob.SIG24, 'g')

pl.subplot(2,2,3, sharex=sp1)
pl.plot(ob.S, ob.SIG12, 'b')
pl.plot(ob.S, ob.SIG14, 'g')
pl.plot(ob.S, ob.SIG23, 'm')
pl.plot(ob.S, ob.SIG34, 'r')



R = ob.SIG11-ob.SIG33
W = ob.SIG11+ob.SIG33
T = R*R+4.*ob.SIG13*ob.SIG13

cos2theta = -np.sign(R)*R/np.sqrt(T)
theta_chk = np.arccos(cos2theta)/2.

sintheta = np.sign((ob.SIG11-ob.SIG33)*ob.SIG13)*np.sqrt(0.5*(1.-cos2theta))
costheta = np.sign((ob.SIG11-ob.SIG33)*ob.SIG13)*np.sqrt(0.5*(1.+cos2theta))

SIG11_hat = 0.5*(W+np.sign(R)*np.sqrt(T))
SIG33_hat = 0.5*(W-np.sign(R)*np.sqrt(T))

# not needed in implementation
theta = np.arctan2(sintheta, costheta)
theta[theta<0]+=np.pi #it is irrelevant


pl.figure(3)
pl.subplot(3,1,1, sharex=sp1)
pl.plot(ob.S, R)
pl.subplot(3,1,2, sharex=sp1)
pl.plot(ob.S, W)
pl.subplot(3,1,3, sharex=sp1)
pl.plot(ob.S, T)


pl.figure(4)
pl.subplot(3,1,1, sharex=sp1)
pl.plot(ob.S, cos2theta)
pl.subplot(3,1,2, sharex=sp1)
pl.plot(ob.S, theta/np.pi*180.)
pl.plot(ob.S, theta_chk/np.pi*180.)
pl.subplot(3,1,3, sharex=sp1)
pl.plot(ob.S, SIG11_hat, '.-')
pl.plot(ob.S, SIG33_hat, '.-')



pl.show()
