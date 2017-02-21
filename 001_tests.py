import python_wrapper as pw

import numpy as np
import matplotlib.pyplot as pl

# Prepare a function that can act on sequences of points
def vectorized_weak_strong(x, y, sigmax, sigmay):

    assert(len(x)==len(y))
    
    part = type('', (), {})() #create empty object
    px = 0*x
    py = 0*y
    for ii in xrange(len(x)):
        part.x = x[ii]
        part.y = y[ii]
        part.px = 0.
        part.py = 0.
        pw.weak_strong_single_particle(part, sigmax, sigmay)
        px[ii] = part.px
        py[ii] = part.py
    return px, py
    
def vectorized_efield_gauss_round(x, y, sigma, Delta_x, Delta_y):

    assert(len(x)==len(y))
    
    Ex = 0.*x
    Ey = 0.*x
    
    for ii in xrange(len(x)):
        Ex[ii], Ey[ii] = pw.test_efield_gauss_round(x=x[ii], y=y[ii], sigma=sigma, Delta_x=Delta_x, Delta_y=Delta_y)
        
    return Ex, Ey
        
    
# First plot
sigmax = 0.30000001
sigmay = 0.30000002

pl.close('all')
pl.figure(1)
xx = np.linspace(-20*sigmax, 20*sigmax, 100)
yy = 0*xx
px, py = vectorized_weak_strong(xx, yy, sigmax, sigmay)
pl.plot(xx, px, '.-')
pl.plot(xx, py, '.-')
pl.suptitle('On the x axis')

pl.figure(2)
yy = np.linspace(-20*sigmay, 20*sigmay, 100)
xx = 0*yy
px, py = vectorized_weak_strong(xx, yy, sigmax, sigmay)
pl.plot(yy, px, '.-')
pl.plot(yy, py, '.-')
pl.suptitle('On the y axis')

x_vec = np.linspace(-20*sigmax, 20*sigmax, 1000)
y_vec = np.linspace(-20*sigmay, 20*sigmay, 1001)

Ex_mat = np.zeros((len(x_vec), len(y_vec)))
Ey_mat = 0*Ex_mat

for i_x in xrange(len(x_vec)):
   Ex_mat[i_x, :], Ey_mat[i_x, :] = vectorized_weak_strong(x_vec[i_x]+y_vec*0, y_vec, sigmax, sigmay)
    
pl.figure(10)
sp1 = pl.subplot2grid((2,2),(0,0))
pl.pcolormesh(x_vec, y_vec, Ex_mat.T)
pl.axis('equal')
sp2 = pl.subplot2grid((2,2),(0,1), sharex=sp1, sharey=sp1)
pl.pcolormesh(x_vec, y_vec, Ey_mat.T)
pl.axis('equal')
sp3 = pl.subplot2grid((2,2),(1,0), colspan=2, sharex=sp1, sharey=sp1)
pl.pcolormesh(x_vec, y_vec, np.sqrt(Ex_mat**2+Ey_mat**2).T)
pl.axis('equal')


sigma = 0.3
Delta_x = 0.
Delta_y = 0.
x = np.linspace(-20*sigma, 20*sigma, 100)
y = 0*x
Ex, Ey = vectorized_efield_gauss_round(x, y, sigma, Delta_x, Delta_y)

pl.figure(101)
pl.plot(x, Ex, '.-');
pl.plot(x, Ey, '.-');
pl.suptitle('Round test')
pl.show()
