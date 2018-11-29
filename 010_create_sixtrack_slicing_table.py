import numpy as np
import pylab as pl

import slicing_sixtrack
import slicing as slc
import boost as bo

import mystyle as ms


N_slices_max = 50

sigmaz = 1.
N_part_tot = 1.


phi = 0.
alpha = 0.
cphi = 1.
sphi = 0.
calpha = 1.
salpha = 0.

z_centroids = np.zeros((N_slices_max+1, N_slices_max), dtype=np.float64)

for N_slices in range(0, N_slices_max+1):

    if N_slices == 0 :
        continue

    star = np.array(np.zeros((3,N_slices)), order='F')
    slicing_sixtrack.stsld(star,cphi2=cphi,sphi2=sphi,sigzs=sigmaz,calpha=calpha,salpha=salpha)
    z_centr_strk = star[2,:]
    x_centr_strk = star[0,:]
    y_centr_strk = star[1,:]

    z_centroids[N_slices, :N_slices] = z_centr_strk

import StringIO
sid = StringIO.StringIO()

np.savetxt(sid, z_centroids, delimiter=',', newline='],\n[')

ss = sid.getvalue()
sid.close()

with open('sixtrack_slicing_table.py', 'w') as fid:
    fid.write('\n'.join([
        'import numpy as np'
        ' ',
        'table = np.array([[',
        ss[:-4],
        ']])',
        ]))

