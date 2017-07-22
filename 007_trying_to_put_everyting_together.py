import numpy as np
import pylab as pl

import propagate_sigma_matrix as psm
import boost
import slicing

# Description of the 6D interaction

#crossing plane
alpha = 0.7

#crossing angle
phi = 0.8

#Intensity strong beam
N_part_tot = 1.1e11

#bunch length strong beam (assumed gaussian)
sigmaz = 0.075

# N slices
N_slices = 7

# strong beam shape at the IP
(Sig_11_0, Sig_12_0, Sig_13_0, 
Sig_14_0, Sig_22_0, Sig_23_0, 
Sig_24_0, Sig_33_0, Sig_34_0, Sig_44_0) = (
8.4282060230000004e-06,  1.8590458800000001e-07,  -3.5512334410000001e-06,
 -3.8254462239999997e-08, 4.101510281e-09, -7.5517657920000006e-08,
 -8.1134615060000002e-10, 1.031446898e-05, 1.177863077e-07, 1.3458251810000001e-09)
 
#Coordinates weak particle that I want to treat
x = 1e-3
px = 50e-3
y = 2e-3
py = 27e-3
sigma = 3.
delta = 2e-4
 

###############################
#    Initialization stage     #
###############################

# Prepare data for Lorentz transformation
parboost = boost.ParBoost(phi=phi, alpha=alpha)

# Prepare data with strong beam shape
Sigmas_0 = psm.Sigmas(Sig_11_0, Sig_12_0, Sig_13_0, 
                        Sig_14_0, Sig_22_0, Sig_23_0, 
                        Sig_24_0, Sig_33_0, Sig_34_0, Sig_44_0)
                        
# Boost strong beam shape
Sigmas_0_star = psm.boost_sigmas(Sigmas_0, parboost.cphi)

# Generate info about slices
z_centroids, _, N_part_per_slice = slicing.constant_charge_slicing_gaussian(N_part_tot, sigmaz, N_slices)

# By boosting the strong z and all zeros, I get the transverse coordinates of the strong beam in the ref system of the weak
# (still I need to fully clarify to myself why this happens)
x_slices_star, px_slices_star, y_slices_star, py_slices_star, sigma_slices_star, delta_slices_star = boost.boost(x=0*z_centroids, px=0*z_centroids, 
                        y=0*z_centroids, py=0*z_centroids, sigma=z_centroids, delta=0*z_centroids, parboost=parboost)
                        
                        

###############################
#      Computation stage      #
###############################

# Boost coordinates of the weak beam
x_star, px_star, y_star, py_star, sigma_star, delta_star = boost.boost(x, px, y, py, sigma, delta, parboost)

for i_slice in xrange(N_slices):
    sigma_slice_star = sigma_slices_star[i_slice]
    x_slice_star = x_slices_star[i_slice]
    y_slice_star = y_slices_star[i_slice]
    
    # Identify the Collision Point (CP)
    S = 0.5*(sigma_star - sigma_slice_star)
    
    # Get info on strong beam shape at the CP
    Sig_11_hat, Sig_33_hat, costheta, sintheta,\
        dS_Sig_11_hat, dS_Sig_33_hat, dS_costheta, dS_sintheta,\
        extra_data = psm.propagate_Sigma_matrix(Sigmas_0_star, S)
        
    # Evaluate transverse coordinates of the weake baem w.r.t. the strong beam centroid
    x_bar_star = x_star + px_star*S - x_slice_star
    y_bar_star = y_star + py_star*S - y_slice_star
    
    # Move to the uncoupled frame
    x_bar_hat_star = x_bar_star*costheta +y_bar_star*sintheta
    y_bar_hat_star = -x_bar_star*sintheta +y_bar_star*costheta
    
    
        
    
 
 
 
 


