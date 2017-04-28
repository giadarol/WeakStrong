import boost_sixtrack as bs
import boost as bo
import numpy as np

phi = 0.5#150e-6
alpha = 23*np.pi/180.

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
print '\n\n ***** Inverse boost'
print "Check x", x, track[0], np.abs(x-track[0]) 
print "Check px", px, track[1], np.abs(px-track[1]) 
print "Check y", y, track[2], np.abs(y - track[2]) 
print "Check py", py, track[3], np.abs(py - track[3]) 
print "Check sigma", sigma, track[4], np.abs(sigma - track[4]) 
print "Check delta", delta, track[5], np.abs(delta - track[5]) 
