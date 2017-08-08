import numpy as np
import pylab as pl

import slicing_sixtrack
import slicing as slc
import boost as bo

import mystyle as ms


#~ sigmaz = .075
#~ N_slices = 7
#~ N_part_tot =1e11

sigmaz = 100*.075
N_slices = 50
N_part_tot =1e11

#~ phi = .22
#~ alpha = .7

phi = 0.8
alpha = 0.7

cphi = np.cos(phi)
sphi = np.sin(phi)
calpha = np.cos(alpha)
salpha = np.sin(alpha)


star = np.array(np.zeros((3,N_slices)), order='F')
slicing_sixtrack.stsld(star,cphi2=cphi,sphi2=sphi,sigzs=sigmaz,calpha=calpha,salpha=salpha)
z_centr_strk = star[2,:]
x_centr_strk = star[0,:]
y_centr_strk = star[1,:]

z_centroids, z_cuts, N_part_per_slice = slc.constant_charge_slicing_gaussian(N_part_tot, sigmaz, N_slices)

# sort according to z, head at the first position in the arrays
ind_sorted = np.argsort(z_centroids)[::-1]
z_centroids = np.take(z_centroids, ind_sorted)
N_part_per_slice = np.take(N_part_per_slice, ind_sorted)

# boost slice coordinates
parboost = bo.ParBoost(phi, alpha)
x_st, px_st, y_st, py_st, sigma_st, delta_st = bo.boost(x=0*z_centroids, px=0*z_centroids, 
                        y=0*z_centroids, py=0*z_centroids, sigma=z_centroids, delta=0*z_centroids, parboost=parboost)
                        



pl.close('all')
fontsz = 14
lw = 2
ms.mystyle_arial(fontsz=fontsz, dist_tick_lab=5)


fig1 = pl.figure(1)
fig1.set_facecolor('w')


ax1 = pl.subplot(3,1,1)
pl.plot(sigma_st, '-ob', lw=lw, label='Library')
pl.plot(z_centr_strk, 'vr', label='Sixtrack')
ax2 = pl.subplot(3,1,2, sharex=ax1)
pl.plot(x_st, '-ob', lw=lw)
pl.plot(x_centr_strk, 'vr')
ax3 = pl.subplot(3,1,3, sharex=ax1)
pl.plot(y_st, '-ob', lw=lw)
pl.plot(y_centr_strk, 'vr')

for ax in [ax1,ax2,ax3]:
    ax.grid('on')
    
ax3.set_xlabel('# slice')
ax1.set_ylabel('z* [m]')
ax2.set_ylabel('x* [m]')
ax3.set_ylabel('y* [m]')

ax1.legend(loc='best', prop={'size':fontsz})


fig100 = pl.figure(100)
fig100.set_facecolor('w')
ax21 = pl.subplot(3,1,1, sharex=ax1)
pl.plot((sigma_st - z_centr_strk)/sigma_st, '.-', lw=lw)
ms.sciy()
ax22 = pl.subplot(3,1,2, sharex=ax1)
pl.plot((x_st - x_centr_strk)/x_st, '.-', lw=lw)
ms.sciy()
ax23 = pl.subplot(3,1,3, sharex=ax1)
pl.plot((y_st - y_centr_strk)/x_st, '.-', lw=lw)
ms.sciy()

for ax in [ax21,ax22,ax23]:
    ax.grid('on')

    
ax23.set_xlabel('# slice')
ax21.set_ylabel('Error on z* [m]')
ax22.set_ylabel('Error on x* [m]')
ax23.set_ylabel('Error on y* [m]')

fig100.subplots_adjust(top=.95, hspace=.35)


pl.show()
