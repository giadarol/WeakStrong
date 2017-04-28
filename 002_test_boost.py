import boost_sixtrack as bs
import numpy as np

phi = 0.5#150e-6
alpha = 0.

x = 1e-3
px = 50e-6
y = 2e-3
py = 20e-6
sigma = 1e-2
delta = 0.5e-3

sphi = np.sin(phi)
cphi = np.cos(phi)
tphi = np.tan(phi)
salpha = np.sin(alpha)
calpha = np.cos(alpha)

track = np.array([x, px, y, py, sigma, delta])

bs.boost(sphi,cphi,tphi,salpha,calpha,track)


h = delta + 1. - np.sqrt((1.+delta)**2-px**2-py**2) 

px_st = px/cphi-h*calpha*tphi/cphi
py_st = py/cphi-h*salpha*tphi/cphi
delta_st = delta -px*calpha*tphi-py*salpha*tphi+h*tphi*tphi

pz_st = np.sqrt((1.+delta_st)**2-px_st**2-py_st**2)
hx_st = px_st/pz_st
hy_st = py_st/pz_st
hsigma_st = 1.-(delta_st+1)/pz_st

L11 = 1.+hx_st*calpha*sphi
L12 = hx_st*salpha*sphi
L13 = calpha*tphi

L21 = hy_st*calpha*sphi
L22 = 1.+hy_st*salpha*sphi
L23 = salpha*tphi

L31 = hsigma_st*calpha*sphi
L32 = hsigma_st*salpha*sphi
L33 = 1./cphi

x_st = L11*x + L12*y + L13*sigma
y_st = L21*x + L22*y + L23*sigma
sigma_st = L31*x + L32*y + L33*sigma

print "Check x", x_st, track[0], 
print "Check px", px_st, track[1]
print "Check y", y_st, track[2]
print "Check py", py_st, track[3]
print "Check sigma", sigma_st, track[4]
print "Check delta", delta_st, track[5]
