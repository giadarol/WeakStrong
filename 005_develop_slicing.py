import numpy as np
import pylab as pl
from scipy.integrate import cumtrapz

from scipy.special import erf, erfinv

sigmaz = 2.

z = np.linspace(-5*sigmaz, 5*sigmaz, 10000)

# line desity
lam = 1/(sigmaz*np.sqrt(2*np.pi))*np.exp(-z*z/(2*sigmaz*sigmaz))

# accumulated charge
Q_ana = 0.5+0.5*erf(z/np.sqrt(2)/sigmaz)
Q_num=cumtrapz(lam,z)

pl.close('all')
pl.figure(1)
ax1 = pl.subplot(1,1,1)
pl.plot(z, lam)

pl.figure(2)
pl.plot(z, Q_ana)
pl.plot(z[1:], Q_num, '--')

# check identification of charge cut
Q1 = 0.8
z1 = np.sqrt(2)*sigmaz*erfinv(2*Q1-1.)
pl.axvline(x=z1)
pl.axhline(y=Q1)

# Constant charge slicing
N_slices = 7
Qi = (np.arange(N_slices)/float(N_slices))[1:]
z_cuts = np.sqrt(2)*sigmaz*erfinv(2*Qi-1.)

for z_cut in z_cuts:
    ax1.axvline(x=z_cut, color='k')

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
    ax1.axvline(x=z_cen, color='r')


pl.show()
