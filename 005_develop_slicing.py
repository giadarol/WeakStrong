import numpy as np
import pylab as pl
from scipy.integrate import cumtrapz

from scipy.special import erf, erfinv

import mystyle as ms

sigmaz = .3

z = np.linspace(-5*sigmaz, 5*sigmaz, 10000)

# line desity
lam = 1/(sigmaz*np.sqrt(2*np.pi))*np.exp(-z*z/(2*sigmaz*sigmaz))

# accumulated charge
Q_ana = 0.5+0.5*erf(z/np.sqrt(2)/sigmaz)
Q_num=cumtrapz(lam,z)

pl.close('all')
fontsz = 14
lw = 3
ms.mystyle_arial(fontsz=fontsz, dist_tick_lab=5)

fig1 = pl.figure(1)
fig1.set_facecolor('w')
ax1 = pl.subplot(1,1,1)
pl.plot(z, lam, lw=lw)

pl.figure(2)
pl.plot(z, Q_ana)
pl.plot(z[1:], Q_num, '--')

# check identification of charge cut
Q1 = 0.8
z1 = np.sqrt(2)*sigmaz*erfinv(2*Q1-1.)
pl.axvline(x=z1)
pl.axhline(y=Q1)

# Constant charge slicing
N_slices = 5
Qi = (np.arange(N_slices)/float(N_slices))[1:]
z_cuts = np.sqrt(2)*sigmaz*erfinv(2*Qi-1.)

z_cuts_with_infs = np.array([-np.inf] + list(z_cuts) + [np.inf])
for i_cut, z_cut in enumerate(z_cuts_with_infs):
    if i_cut ==0:
        continue
    mask_z = np.logical_and(z<=z_cuts_with_infs[i_cut], z>z_cuts_with_infs[i_cut-1])
    if np.mod(i_cut, 2)>0:col = 'b'
    else:col = 'grey'
    ax1.fill_between(z[mask_z], lam[mask_z], alpha=.5, color=col)
    #~ ax1.axvline(x=z_cut, color='k')
ax1.grid('on')


z_centroids = []
first_centroid = -sigmaz/np.sqrt(2*np.pi)*np.exp(-z_cuts[0]**2/(2*sigmaz*sigmaz))*float(N_slices)
z_centroids.append(first_centroid)
for ii in xrange(N_slices-2):
    this_centroid = -sigmaz/np.sqrt(2*np.pi)*(np.exp(-z_cuts[ii+1]**2/(2*sigmaz*sigmaz))-
                                              np.exp(-z_cuts[ii]**2/(2*sigmaz*sigmaz)))*float(N_slices) #the multiplication times n slices comes from the fact that we have to divide by the slice charge, i.e. 1./N
    z_centroids.append(this_centroid)                                         

last_centroid = sigmaz/np.sqrt(2*np.pi)*np.exp(-z_cuts[-1]**2/(2*sigmaz*sigmaz))*float(N_slices)
z_centroids.append(last_centroid)


for z_cen in z_centroids:
    ax1.axvline(x=z_cen, color='k', lw=lw, linestyle='--')
ax1.set_ylim(bottom=0)
ax1.set_ylabel('Longitudinal profile')
ax1.set_xlabel('z [m]')
pl.show()
