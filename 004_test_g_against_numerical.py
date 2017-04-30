import transverse_efields as tef
import beambeam_force_sixtrack as bbfs

from scipy.constants import epsilon_0

import numpy as np
import pylab as pl

nv = np.vectorize

min_sigma_diff = 1e-10

#~ # Tall beam
#~ sigma_x=.5
#~ sigma_y=.9
#~ theta = 20*np.pi/180
#~ r_max = 20.

#~ # Flat beam
#~ sigma_x=.9
#~ sigma_y=.5
#~ theta = 20*np.pi/180
#~ r_max = 20.

# Round beam
sigma_x=.5
sigma_y=.5
theta = 20*np.pi/180
r_max = 20.



r = np.linspace(-r_max, r_max, 10000)
x = r * np.cos(theta)
y = r * np.sin(theta)
r_centers = 0.5*(r[1:]+r[:-1])

D_sigma_x = sigma_x*1e-2
D_sigma_y = sigma_y*1e-2

Ex, Ey, Gx, Gy = nv(tef.get_Ex_Ey_Gx_Gy_gauss)(x, y, sigma_x, sigma_y, min_sigma_diff)
phi = -np.cumsum(0.5*(Ex[:-1]+Ex[1:])*np.diff(x)+0.5*(Ey[:-1]+Ey[1:])*np.diff(y))

Ex_plusx, Ey_plusx, _, _ = nv(tef.get_Ex_Ey_Gx_Gy_gauss)(x, y, sigma_x + D_sigma_x, sigma_y, min_sigma_diff)
phi_plusx = -np.cumsum(0.5*(Ex_plusx[:-1]+Ex_plusx[1:])*np.diff(x)+0.5*(Ey_plusx[:-1]+Ey_plusx[1:])*np.diff(y))
Ex_minusx, Ey_minusx, _, _ = nv(tef.get_Ex_Ey_Gx_Gy_gauss)(x, y, sigma_x - D_sigma_x, sigma_y, min_sigma_diff)
phi_minusx = -np.cumsum(0.5*(Ex_minusx[:-1]+Ex_minusx[1:])*np.diff(x)+0.5*(Ey_minusx[:-1]+Ey_minusx[1:])*np.diff(y))
Gx_num = -(phi_plusx-phi_minusx)/((sigma_x+D_sigma_x)**2 - (sigma_x-D_sigma_x)**2) #The derivative is with respect to the capital Sigma!!!

Ex_plusy, Ey_plusy, _, _ = nv(tef.get_Ex_Ey_Gx_Gy_gauss)(x, y, sigma_x, sigma_y + D_sigma_y, min_sigma_diff)
phi_plusy = -np.cumsum(0.5*(Ex_plusy[:-1]+Ex_plusy[1:])*np.diff(x)+0.5*(Ey_plusy[:-1]+Ey_plusy[1:])*np.diff(y))
Ex_minusy, Ey_minusy, _, _ = nv(tef.get_Ex_Ey_Gx_Gy_gauss)(x, y, sigma_x, sigma_y - D_sigma_y, min_sigma_diff)
phi_minusy = -np.cumsum(0.5*(Ex_minusy[:-1]+Ex_minusy[1:])*np.diff(x)+0.5*(Ey_minusy[:-1]+Ey_minusy[1:])*np.diff(y))
Gy_num = -(phi_plusy-phi_minusy)/((sigma_y+D_sigma_y)**2 - (sigma_y-D_sigma_y)**2)


pl.close('all')
pl.figure(1)
pl.plot(r_centers, phi)
pl.plot(r_centers, phi_plusx)
pl.plot(r_centers, phi_minusx)

pl.figure(2)
pl.subplot(2,1,1)
pl.plot(r_centers, phi_plusx-phi)
pl.plot(r_centers, phi_minusx-phi)
pl.subplot(2,1,2)
pl.plot(r_centers, phi_plusy-phi)
pl.plot(r_centers, phi_minusy-phi)



pl.figure(3)
pl.subplot(2,1,1)
pl.plot(r_centers, Gx_num)
pl.plot(r_centers, 0.5*(Gx[:-1]+Gx[1:]), '--')
pl.subplot(2,1,2)
pl.plot(r_centers, Gy_num)
pl.plot(r_centers, 0.5*(Gy[:-1]+Gy[1:]), '--')



#~ pl.figure(10)
#~ pl.plot(r, Ex)
#~ pl.plot(r, Ey)


pl.show()
