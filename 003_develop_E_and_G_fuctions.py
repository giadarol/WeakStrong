import transverse_efields as tef
import beambeam_force_sixtrack as bbfs

import numpy as np
import pylab as pl

nv = np.vectorize


sigma_x=.9
sigma_y=.9
theta = 20*np.pi/180
r_max = 20.

flag_round = np.abs(sigma_x-sigma_y)<1e-10


r = np.linspace(-r_max, r_max, 1000)

x = r * np.cos(theta)
y = r * np.sin(theta)

if flag_round:
    print "Round beam"
    Ex_sl, Ey_sl = tef.vect_transverse_field_round(x, y, sigma_x, Delta_x=0., Delta_y=0.)
else:
    print "Elliptic beam"
    Ex_sl, Ey_sl = tef.vect_transverse_field_ellip(x, y, sigma_x+1e-4, sigma_y, Delta_x=0., Delta_y=0.)


bbfx, bbfy, bbgx, bbgy = nv(bbfs.bbf)(x,y,sigma_x**2,sigma_y**2)
from scipy.constants import epsilon_0
fact_sixtrack = 1./(4*np.pi*epsilon_0)
bbfx*=fact_sixtrack
bbfy*=fact_sixtrack
bbgx*=fact_sixtrack
bbgy*=fact_sixtrack


pl.close('all')
pl.figure(1)
ax1 = pl.subplot(2,1,1)
pl.plot(r, Ex_sl)
pl.plot(r, bbfx)
pl.ylabel('Ex')
ax2 = pl.subplot(2,1,2, sharex=ax1)
pl.plot(r, Ey_sl)
pl.plot(r, bbfy)
pl.ylabel('Ey')


#~ pl.figure(2)
#~ pl.plot(r, np.abs(Ex_sl-bbfx)/np.abs(bbfx))

pl.show()
