import numpy as np
import pylab as pl

import propagate_sigma_matrix as psm
import boost
import slicing
import transverse_efields as tef

from scipy.constants import e as qe
from scipy.constants import c as c_light

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

# Single particle properties
q0 = qe
qsl = qe
p0 = 6.5e12 *qe/c_light

# Minimum difference to fall on round
min_sigma_diff = 1e-16



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

p0c = p0*c_light

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
                        

# Record coordinates before interaction (for comparison against sixtrack)   
coord_init = np.array([x, px, y, py, sigma, delta])                     

###############################
#      Computation stage      #
###############################

# Boost coordinates of the weak beam
x_star, px_star, y_star, py_star, sigma_star, delta_star = boost.boost(x, px, y, py, sigma, delta, parboost)
#~ print x_star, px_star, y_star, py_star, sigma_star, delta_star
for i_slice in xrange(N_slices):
    sigma_slice_star = sigma_slices_star[i_slice]
    x_slice_star = x_slices_star[i_slice]
    y_slice_star = y_slices_star[i_slice]
    
    # Compute force scaling factor
    Ksl = N_part_per_slice[i_slice]*qsl*q0/p0c
    
    # Identify the Collision Point (CP)
    S = 0.5*(sigma_star - sigma_slice_star)
    
    # Get strong beam shape at the CP
    Sig_11_hat_star, Sig_33_hat_star, costheta, sintheta,\
        dS_Sig_11_hat_star, dS_Sig_33_hat_star, dS_costheta, dS_sintheta,\
        extra_data = psm.propagate_Sigma_matrix(Sigmas_0_star, S)
        
    # Evaluate transverse coordinates of the weake baem w.r.t. the strong beam centroid
    x_bar_star = x_star + px_star*S - x_slice_star
    y_bar_star = y_star + py_star*S - y_slice_star
    
    # Move to the uncoupled reference frame
    x_bar_hat_star = x_bar_star*costheta +y_bar_star*sintheta
    y_bar_hat_star = -x_bar_star*sintheta +y_bar_star*costheta
    
    # Compute derivatives of the transformation
    dS_x_bar_hat_star = x_bar_star*dS_costheta +y_bar_star*dS_sintheta
    dS_y_bar_hat_star = -x_bar_star*dS_sintheta +y_bar_star*dS_costheta
    
    # Compute normalized field
    Ex, Ey, Gx, Gy = tef.get_Ex_Ey_Gx_Gy_gauss(x=x_bar_star, y=y_bar_star, 
                        sigma_x=np.sqrt(Sig_11_hat_star), sigma_y=np.sqrt(Sig_33_hat_star),
                        min_sigma_diff = min_sigma_diff)
                        
    # Compute kicks
    Fx_hat_star = Ksl*Ex
    Fy_hat_star = Ksl*Ey
    Gx_hat_star = Ksl*Gx
    Gy_hat_star = Ksl*Gy
    
    # Move kisks to coupled reference frame
    Fx_star = Fx_hat_star*costheta - Fy_hat_star*sintheta
    Fy_star = Fx_hat_star*sintheta + Fy_hat_star*costheta
    
    # Compute longitudinal kick
    Fz_star = 0.5*(Fx_hat_star*dS_x_bar_hat_star  + Fy_hat_star*dS_y_bar_hat_star+\
                   Gx_hat_star*dS_Sig_11_hat_star + Gy_hat_star*dS_Sig_33_hat_star)
                   
    # Apply the kicks (Hirata's synchro-beam)
    delta_star = delta_star + Fz_star+0.5*(\
                Fx_star*(px_star+0.5*Fx_star)+\
                Fy_star*(py_star+0.5*Fy_star))
    x_star = x_star - S*Fx_star
    px_star = px_star + Fx_star
    y_star = y_star - S*Fy_star
    py_star = py_star + Fy_star
    
# Inverse boost on the coordinates of the weak beam
x, px, y, py, sigma, delta = boost.inv_boost(x_star, px_star, y_star, py_star, sigma_star, delta_star, parboost)

coord_fin = np.array([x, px, y, py, sigma, delta])  

###############################################
## Use sixtrack to make the same interaction ##
###############################################

from scipy.constants import epsilon_0
from scipy.constants import e as qe

npa = 1

track = coord_init.copy()

paramarr = np.float_(np.zeros(18))
paramarr[0] = phi
paramarr[1] = N_slices
paramarr[2] = alpha
paramarr[3] = qe*qe/(4*np.pi*epsilon_0)*N_part_tot ##f
paramarr[17] = phi ##phi2

param = np.float_(np.array([paramarr], order='F'))


sigzs = sigmaz

bcu = np.float_(np.array([[\
Sigmas_0.Sig_11_0,
Sigmas_0.Sig_33_0,
Sigmas_0.Sig_13_0,
Sigmas_0.Sig_12_0,
Sigmas_0.Sig_14_0,
Sigmas_0.Sig_22_0,
Sigmas_0.Sig_23_0,
Sigmas_0.Sig_24_0,
Sigmas_0.Sig_34_0,
Sigmas_0.Sig_44_0,
0.,# not used
0.,]], order='F'))

ibb = 1
ne = 1
ibtyp = 0
ibbc=1
mbea = 1
beam_expflag = 0
pieni = 1e-30


import full_interaction_sixtrack as fis
fis.beamint(np=npa,track=track,param=param,sigzs=sigzs,bcu=bcu,ibb=ibb,
            ne=ne,ibtyp=ibtyp,ibbc=ibbc,mbea=mbea,beam_expflag=beam_expflag,pieni=pieni, \
            npart=1, nele=1, nbb=1)
    
    

    
    
        
    
 
 
 
 


