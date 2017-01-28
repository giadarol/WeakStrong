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

# First plot
sigmax = 1.
sigmay = 2.

xx = np.linspace(-20*sigmax, 20*sigmax, 100)
yy = 0*xx

px, py = vectorized_weak_strong(xx, yy, sigmax, sigmay)

pl.close('all')
pl.figure(1)
pl.plot(xx, px, '.-')
pl.show()
