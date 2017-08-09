import transverse_efields as tef
import from_PyECLOUD.BassErsk as PyECBE

import numpy as np
import pylab as pl
import mystyle as ms


#~ sigma_x=.9
#~ sigma_y=.7

sigma_x=.5
sigma_y=.9

r = np.linspace(-20, 20., 1000)
theta = 20*np.pi/180
x = r * np.cos(theta)
y = r * np.sin(theta)

nv = np.vectorize


Ex_PyEC, Ey_PyEC = nv(PyECBE.BassErsk)(x,y, sigma_x, sigma_y)
Ex_sl, Ey_sl = nv(tef.transverse_field_ellip)(x, y, sigma_x, sigma_y, Delta_x=0., Delta_y=0.)


pl.close('all')
fontsz = 14
lw = 3
ms.mystyle_arial(fontsz=fontsz, dist_tick_lab=5)

fig1 = pl.figure(1)
fig1.set_facecolor('w')
ax1 = pl.subplot(2,1,1)
pl.plot(r, Ex_PyEC, 'b', lw = lw, label='PyECLOUD')
pl.plot(r, Ex_sl, '--r', lw = lw, label='Library')
ax2 = pl.subplot(2,1,2,sharex = ax1)
pl.plot(r, Ey_PyEC, 'b', lw = lw)
pl.plot(r, Ey_sl, '--r', lw = lw)

for sp in [ax1, ax2]:
    sp.grid('on')
ax1.legend(loc='best', prop={'size':fontsz})
ax2.set_xlabel('r [m]')
ax1.set_ylabel('Ex [V/m]')
ax2.set_ylabel('Ey [V/m]')

pl.suptitle('sigmax = %.1e m sigma_y=%.1e m, theta_dir=%.1f deg\nFields for 1.0 C/m'%(sigma_x, sigma_y, theta*180/np.pi))

pl.show()

#~ Ex, Ey = tef.transverse_field_ellip(x, y, sigma_x, sigma_y, Delta_x, Delta_y)
                        #~ 
#~ print '%e, %e'%(Ex, Ey)

