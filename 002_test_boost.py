import boost_sixtrack as bs
import boost as bo
import numpy as np

phi = 0.8#150e-6
alpha = 23*np.pi/180.

x = 1e-3
px = 50e-6*10000
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

# test direct boost
bs.boost(sphi,cphi,tphi,salpha,calpha,track)
parboost = bo.ParBoost(phi, alpha)
x_st, px_st, y_st, py_st, sigma_st, delta_st = bo.boost(x, px, y, py, sigma, delta, parboost)

print '\n\n ***** Direct boost'
print "Check x", x_st, track[0], np.abs(x_st-track[0]) 
print "Check px", px_st, track[1], np.abs(px_st-track[1]) 
print "Check y", y_st, track[2], np.abs(y_st - track[2]) 
print "Check py", py_st, track[3], np.abs(py_st - track[3]) 
print "Check sigma", sigma_st, track[4], np.abs(sigma_st - track[4]) 
print "Check delta", delta_st, track[5], np.abs(delta_st - track[5]) 


# test inverse boost (sixtrack vs initial)
bs.boosti(sphi,cphi,tphi,salpha,calpha,track)

print '\n\n ***** Inverse boost (sixtrack vs initial)'
print "Check x", x, track[0], np.abs(x-track[0]) 
print "Check px", px, track[1], np.abs(px-track[1]) 
print "Check y", y, track[2], np.abs(y - track[2]) 
print "Check py", py, track[3], np.abs(py - track[3]) 
print "Check sigma", sigma, track[4], np.abs(sigma - track[4]) 
print "Check delta", delta, track[5], np.abs(delta - track[5]) 


pz_st = np.sqrt((1.+delta_st)**2-px_st**2-py_st**2)
hx_st = px_st/pz_st
hy_st = py_st/pz_st
hsigma_st = 1.-(delta_st+1)/pz_st

Det_L = 1./cphi + (hx_st*calpha + hy_st*salpha-hsigma_st*sphi)*tphi

Linv_11 = (1./cphi + salpha*tphi*(hy_st-hsigma_st*salpha*sphi))/Det_L
Linv_12 = (salpha*tphi*(hsigma_st*calpha*sphi-hx_st))/Det_L
Linv_13 = -tphi*(calpha - hx_st*salpha*salpha*sphi + hy_st*calpha*salpha*sphi)/Det_L

Linv_21 = (calpha*tphi*(-hy_st + hsigma_st*salpha*tphi*sphi))/Det_L
Linv_22 = (1./cphi + calpha*tphi*(hx_st-hsigma_st*calpha*calpha*sphi))/Det_L
Linv_23 = -tphi*(salpha - hy_st*calpha*calpha*sphi + hx_st*calpha*salpha*sphi)/Det_L

Linv_31 = -hsigma_st*calpha*sphi/Det_L
Linv_32 = -hsigma_st*salpha*sphi/Det_L
Linv_33 = (1. + hx_st*calpha*sphi + hy_st*salpha*sphi)/Det_L

x_i = Linv_11*x_st + Linv_12*y_st + Linv_13*sigma_st
y_i = Linv_21*x_st + Linv_22*y_st + Linv_23*sigma_st
sigma_i = Linv_31*x_st + Linv_32*y_st + Linv_33*sigma_st

h = (delta_st+1.-pz_st)*cphi*cphi

px_i = px_st*cphi+h*calpha*tphi
py_i = py_st*cphi+h*salpha*tphi

delta_i = delta_st + px_i*calpha*tphi + py*salpha*tphi - h*tphi*tphi


# test inverse boost (python vs initial)
print '\n\n ***** Inverse boost (python vs initial)'
print "Check x", x, x_i, np.abs(x-x_i) 
print "Check px", px, px_i, np.abs(px-px_i) 
print "Check y", y, y_i, np.abs(y-y_i) 
print "Check py", py, py_i, np.abs(py-py_i) 
print "Check sigmas", sigma, sigma_i, np.abs(sigma-sigma_i) 
print "Check delta", delta, delta_i, np.abs(delta-delta_i) 

# test inverse boost (python vs sixtracks)
print '\n\n ***** Inverse boost (python vs sixtrak)'
print "Check x", x_i, track[0], np.abs(x_i-track[0]) 
print "Check px", px_i, track[1], np.abs(px_i-track[1]) 
print "Check y", y_i, track[2], np.abs(y_i-track[2]) 
print "Check py", py_i, track[3], np.abs(py_i-track[3]) 
print "Check sigmas", sigma_i, track[4], np.abs(sigma_i-track[4]) 
print "Check delta", delta_i, track[5], np.abs(delta_i-track[5]) 
