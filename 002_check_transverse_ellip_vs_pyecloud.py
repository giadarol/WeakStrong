import transverse_efields as tef
import from_PyECLOUD.BassErsk as PyECBE

import numpy as np
import pylab as pl


#~ sigma_x=.9
#~ sigma_y=.7

sigma_x=.5
sigma_y=.9

r = np.linspace(-20, 20., 1000)
theta = 20*np.pi/180
x = r * np.cos(theta)
y = r * np.sin(theta)


Ex_PyEC, Ey_PyEC = np.vectorize(PyECBE.BassErsk)(x,y, sigma_x, sigma_y)
Ex_sl, Ey_sl = tef.vect_transverse_field_ellip(x, y, sigma_x, sigma_y, Delta_x=0., Delta_y=0.)


pl.close('all')
pl.figure(1)
ax1 = pl.subplot(2,1,1)
pl.plot(r, Ex_PyEC)
pl.plot(r, Ex_sl, '.r')
pl.subplot(2,1,2,sharex = ax1)
pl.plot(r, Ey_PyEC)
pl.plot(r, Ey_sl, '.r')
pl.show()

#~ Ex, Ey = tef.transverse_field_ellip(x, y, sigma_x, sigma_y, Delta_x, Delta_y)
                        #~ 
#~ print '%e, %e'%(Ex, Ey)

