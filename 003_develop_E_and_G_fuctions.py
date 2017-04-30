import transverse_efields as tef
import beambeam_force_sixtrack as bbfs

from scipy.constants import epsilon_0

import numpy as np
import pylab as pl

nv = np.vectorize


sigma_x=.5
sigma_y=.9
theta = 23*np.pi/180
r_max = 20.

flag_round = np.abs(sigma_x-sigma_y)<1e-10


r = np.linspace(-r_max, r_max, 1000)

x = r * np.cos(theta)
y = r * np.sin(theta)

if flag_round:
    print "Round beam"
    
    sigma = 0.5*(sigma_x+sigma_y)
    
    #electric fields from sixtracklib
    Ex, Ey = tef.vect_transverse_field_round(x, y, sigma, Delta_x=0., Delta_y=0.)
    Gx = 1/(2.*(x*x+y*y))*(y*Ey-x*Ex+1./(2*np.pi*epsilon_0*sigma*sigma)*x*x*np.exp(-(x**2+y**2)/(2.*sigma*sigma)))
    Gy = 1./(2*(x*x+y*y))*(x*Ex-y*Ey+1./(2*np.pi*epsilon_0*sigma*sigma)*y*y*np.exp(-(x**2+y**2)/(2.*sigma*sigma)))
    
    
else:
    print "Elliptic beam"
    Ex, Ey = tef.vect_transverse_field_ellip(x, y, sigma_x, sigma_y, Delta_x=0., Delta_y=0.)
    Sig_11 = sigma_x*sigma_x
    Sig_33 = sigma_y*sigma_y
    
    Gx =-1./(2*(Sig_11-Sig_33))*(x*Ex+y*Ey+1./(2*np.pi*epsilon_0)*\
                (sigma_y/sigma_x*np.exp(-x*x/(2*Sig_11)-y*y/(2*Sig_33))-1.))
    Gy =1./(2*(Sig_11-Sig_33))*(x*Ex+y*Ey+1./(2*np.pi*epsilon_0)*\
                (sigma_x/sigma_y*np.exp(-x*x/(2*Sig_11)-y*y/(2*Sig_33))-1.))


# In sixtrack it is used in this way
if sigma_x>sigma_y:
    bbfx, bbfy, bbgx, bbgy = nv(bbfs.bbf)(x,y,sigma_x**2,sigma_y**2)
else:
    bbfy, bbfx, bbgy, bbgx = nv(bbfs.bbf)(y,x,sigma_y**2,sigma_x**2)


fact_sixtrack = 1./(4*np.pi*epsilon_0)
bbfx*=fact_sixtrack
bbfy*=fact_sixtrack
bbgx*=fact_sixtrack
bbgy*=fact_sixtrack


pl.close('all')
pl.figure(1)
ax1 = pl.subplot(2,1,1)
pl.plot(r, Ex)
pl.plot(r, bbfx)
pl.ylabel('Ex')
ax2 = pl.subplot(2,1,2, sharex=ax1)
pl.plot(r, Ey)
pl.plot(r, bbfy)
pl.ylabel('Ey')


pl.figure(2)
ax1 = pl.subplot(2,1,1)
pl.plot(r, Gx)
pl.plot(r, bbgx, label='sixtrack')
pl.ylabel('Gx')
pl.grid('on')
ax2 = pl.subplot(2,1,2, sharex=ax1)
pl.plot(r, Gy)
pl.plot(r, bbgy, label='sixtrack')
pl.ylabel('Gy')
pl.grid('on')
pl.legend()
#~ pl.figure(3)
#~ pl.plot(r, np.abs(Ex_sl-bbfx)/np.abs(bbfx))

#~ pl.figure(4)
#~ pl.subplot(111, sharex=ax1)
#~ pl.plot(r, Gx/(bbgx))
#~ pl.ylim(0,2)

pl.show()
