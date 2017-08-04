import BB6D
from scipy.constants import e as qe
from scipy.constants import c as c_light

import numpy as np

#crossing plane
alpha = 0.7

#crossing angle
phi = 0.8

#Intensity strong beam
N_part_tot = 1.1e15

#bunch length strong beam (assumed gaussian)
sigmaz = 0.075*100

# N slices
N_slices = 50

# Single particle properties
q_part = qe

# Minimum difference to fall on round
min_sigma_diff = 1e-16

threshold_singular = 1e-16

 
# strong beam shape at the IP (decoupled round beam)
(Sig_11_0, Sig_12_0, Sig_13_0, 
Sig_14_0, Sig_22_0, Sig_23_0, 
Sig_24_0, Sig_33_0, Sig_34_0, Sig_44_0) = (
20e-06,  0.,  0.,
0., 0., 0.,
0., 20e-6, 0., 0.)




bb6ddata = BB6D.BB6D_init(q_part, N_part_tot, sigmaz, N_slices, min_sigma_diff, threshold_singular, 
                phi, alpha, 
                Sig_11_0, Sig_12_0, Sig_13_0, 
                Sig_14_0, Sig_22_0, Sig_23_0, 
                Sig_24_0, Sig_33_0, Sig_34_0, Sig_44_0)

self = bb6ddata

def int_to_float64arr(val):
    temp = np.zeros(1, (np.float64, {'i64':('i8',0)}))
    temp['i64'][0] = val
    return temp

buffer_list = []
# Buffers corresponding to BB6D struct
buffer_list.append(np.array([self.q_part], dtype=np.float64))
buffer_list.append(self.parboost.tobuffer())
buffer_list.append(self.Sigmas_0_star.tobuffer())
buffer_list.append(np.array([self.min_sigma_diff], dtype=np.float64))
buffer_list.append(np.array([self.threshold_singular], dtype=np.float64))
buffer_list.append(int_to_float64arr(self.N_slices))
buffer_list.append(int_to_float64arr(3))# offset to N_part_per_slice
buffer_list.append(int_to_float64arr(2+N_slices))# offset to x_slices_star
buffer_list.append(int_to_float64arr(1+2*N_slices))# offset to y_slices_star
buffer_list.append(int_to_float64arr(0+3*N_slices))# offset to sigma_slices_star

# Buffers corresponding to arrays
buffer_list.append(np.array(self.N_part_per_slice, dtype=np.float64))
buffer_list.append(np.array(self.x_slices_star, dtype=np.float64))
buffer_list.append(np.array(self.y_slices_star, dtype=np.float64))
buffer_list.append(np.array(self.sigma_slices_star, dtype=np.float64))

buf = np.concatenate(buffer_list)

### Try to call the function

import os, ctypes
modulepath=os.path.dirname(os.path.abspath(__file__))
libpath=os.path.join(modulepath, 'cBB6D.so')
BB6D=ctypes.CDLL(libpath)
                
#Coordinates of and other properties weak particle that I want to treat
x = 1e-3
px = 50e-3
y = 2e-3
py = 27e-3
sigma = 3.
delta = 2e-4
q0 = qe
p0 = 6.5e12 *qe/c_light

coord_init = np.array([x, px, y, py, sigma, delta])                     

x_c = ctypes.c_double(x)
px_c = ctypes.c_double(px)
y_c = ctypes.c_double(y)
py_c = ctypes.c_double(py)
sigma_c = ctypes.c_double(sigma)
delta_c = ctypes.c_double(delta)

BB6D.BB6D_track(ctypes.byref(x_c),
        ctypes.byref(px_c),
        ctypes.byref(y_c),
        ctypes.byref(py_c),
        ctypes.byref(sigma_c),
        ctypes.byref(delta_c), 
        ctypes.c_double(q0), 
        ctypes.c_double(p0), 
        buf.ctypes.data)
        
(x, px, y, py, sigma, delta) = (x_c.value,
        px_c.value,
        y_c.value,
        py_c.value,
        sigma_c.value,
        delta_c.value) 

coord_fin = np.array([x, px, y, py, sigma, delta]) 

        
###############################################
## Use sixtrack to make the same interaction ##
###############################################

from scipy.constants import epsilon_0
from scipy.constants import e as qe

sixtrack_slicing = False

npa = 1

track = coord_init.copy()

paramarr = np.float_(np.zeros(18))
paramarr[0] = phi
paramarr[1] = N_slices
paramarr[2] = alpha
paramarr[3] = -qe*qe/(4*np.pi*epsilon_0)*N_part_tot/p0/c_light ##f, sixtrack has a minus when it applies the sign
paramarr[17] = phi ##phi2

param = np.float_(np.array([paramarr], order='F'))


sigzs = sigmaz

bcu = np.float_(np.array([[\
Sig_11_0,
Sig_33_0,
Sig_13_0,
Sig_12_0,
Sig_14_0,
Sig_22_0,
Sig_23_0,
Sig_24_0,
Sig_34_0,
Sig_44_0,
0.,# not used
0.,]], order='F'))

ibb = 1
ne = 1
ibtyp = 0
ibbc=1
mbea = 1
beam_expflag = 0
pieni = 1e-30

star_input = np.array([\
    bb6ddata.x_slices_star,
    bb6ddata.y_slices_star,
    bb6ddata.sigma_slices_star], order='F')

if sixtrack_slicing:
    star_input = star_input*0.+9999.

import full_interaction_sixtrack as fis
fis.beamint(np=npa,track=track,param=param,sigzs=sigzs,bcu=bcu,ibb=ibb,
            ne=ne,ibtyp=ibtyp,ibbc=ibbc,mbea=mbea,beam_expflag=beam_expflag,pieni=pieni, \
            npart=1, nele=1, nbb=1, star_input=star_input)

#Compare the two
print '\nCompare kicks against sixtrack:'
names_list = 'x px y py sigma delta'.split()
for name, err, err_sixtr in zip(names_list, coord_fin-coord_init, track-coord_init):
    print 'D_'+name, err, err_sixtr



