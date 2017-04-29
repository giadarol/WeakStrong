import transverse_efields as tef

x = 0.3
y = 0.4
sigma_x=.9
sigma_y=.7
Delta_x = 0.
Delta_y = 0.

Ex, Ey = tef.transverse_field_ellip(x, y, sigma_x, sigma_y, Delta_x, Delta_y)
                        
print '%e, %e'%(Ex, Ey)

