import transverse_efields as tef
import beambeam_force_sixtrack as bbfs

import numpy as np
import pylab as pl

nv = np.vectorize


sigma_x=.9
sigma_y=.9


r = np.linspace(-20, 20., 1000)
theta = 0.000*np.pi/180
x = r * np.cos(theta)
y = r * np.sin(theta)

Ex_sl, Ey_sl = tef.vect_transverse_field_ellip(x, y, sigma_x+1e-4, sigma_y, Delta_x=0., Delta_y=0.)
bbfx, bbfy, bbgx, bbgy = nv(bbfs.bbf)(x,y,sigma_x,sigma_y)

from scipy.constants import epsilon_0
fact_sixtrack = 1./(4*np.pi*epsilon_0)

bbfx = fact_sixtrack * bbfx

pl.close('all')
pl.figure(1)
pl.plot(r, Ex_sl)
pl.plot(r, bbfx)

pl.figure(2)
pl.plot(r, Ex_sl/bbfx)

pl.show()
