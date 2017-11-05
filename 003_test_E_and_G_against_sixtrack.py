import transverse_efields_implem_c as tef
import beambeam_force_sixtrack as bbfs
import mystyle as ms

from scipy.constants import epsilon_0

import numpy as np
import pylab as pl

nv = np.vectorize

min_sigma_diff = 1e-10

# Tall beam
sigma_x=.5
sigma_y=.9
theta = 23*np.pi/180
r_max = 20.

#~ # Flat beam
#~ sigma_x=.9
#~ sigma_y=.5
#~ theta = 23*np.pi/180
#~ r_max = 20.

# Round beam
sigma_x=.7
sigma_y=.7
theta = 23*np.pi/180
r_max = 20.

r = np.linspace(-r_max, r_max, 1000)
x = r * np.cos(theta)
y = r * np.sin(theta)

Ex, Ey, Gx, Gy = nv(tef.get_Ex_Ey_Gx_Gy_gauss)(x, y, sigma_x, sigma_y, min_sigma_diff)

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


fontsz = 14
lw = 3
pl.close('all')
ms.mystyle_arial(fontsz=fontsz, dist_tick_lab=5)

fig1 = pl.figure(1)
fig1.set_facecolor('w')

ax1 = pl.subplot(2,1,1)
pl.plot(r, Ex, 'b', label = 'Library', lw=lw)
pl.plot(r, bbfx, 'r--',  label='Sixtrack', lw=lw)
pl.ylabel('Ex')
ax2 = pl.subplot(2,1,2, sharex=ax1)
pl.plot(r, Ey, 'b', label = 'Library', lw=lw)
pl.plot(r, bbfy, 'r--',  label='Sixtrack', lw=lw)
pl.ylabel('Ey')
for ax in [ax1, ax2]: ax.grid('on')
ax1.legend(loc='best', prop={'size':fontsz})


fig2 = pl.figure(2)
fig2.set_facecolor('w')
ax21 = pl.subplot(2,1,1)
pl.plot(r, Gx, 'b', label = 'Library', lw=lw)
pl.plot(r, bbgx, 'r--',  label='Sixtrack', lw=lw)
pl.ylabel('Gx')
pl.grid('on')
ax22 = pl.subplot(2,1,2, sharex=ax1)
pl.plot(r, Gy, 'b', lw=lw)
pl.plot(r, bbgy, 'r--', label='sixtrack', lw=lw)
pl.ylabel('Gy')
pl.grid('on')

ax21.legend(loc='best', prop={'size':fontsz})

for fig in [fig1, fig2]:
    fig.suptitle('sigmax = %.1e m sigma_y=%.1e m, theta_dir=%.1f deg\nFor 1.0 C/m'%(sigma_x, sigma_y, theta*180/np.pi))
#~ pl.figure(3)
#~ pl.plot(r, np.abs(Ex_sl-bbfx)/np.abs(bbfx))

#~ pl.figure(4)
#~ pl.subplot(111, sharex=ax1)
#~ pl.plot(r, Gx/(bbgx))
#~ pl.ylim(0,2)

pl.show()
