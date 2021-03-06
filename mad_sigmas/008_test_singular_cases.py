import numpy as np
import sys, os
import pylab as pl

sys.path.append('../')

save_figs = True
mode = 'check_singularities'#'check_singularities'/'vs_sixtrack'

import mystyle as ms
import propagate_sigma_matrix_implem_c as psm # C implem
#~ import propagate_sigma_matrix as psm # Python implem

# Case T=0., |c|>0.
SIG11 = 10.
SIG33 = 10.
SIG13 = 0. 
SIG12 = -.5
SIG22 = 0.2
SIG14 = 0.5
SIG23 = 0.7
SIG24 = 0.1
SIG34 = -0.9
SIG44 = 0.2

# Case T=0., c = 0., |a|>0
SIG11 = 10.
SIG33 = 10.
SIG13 = 0. 
SIG12 = -.5
SIG22 = 0.2
SIG14 = 0.5
SIG23 = -0.5
SIG24 = -0.1
SIG34 = 0.15
SIG44 = 0.25

# Case T=0., c = 0., |a|>0, d = 0 decoupled case
SIG11 = 10
SIG33 = 10
SIG13 = 0. 
SIG12 = -.5
SIG22 = 0.2
SIG14 = 0.5
SIG23 = -0.5
SIG24 = 0.
SIG34 = 0.
SIG44 = 0.25

# Case T=0., c = 0., a = 0.
SIG11 = 10.
SIG33 = 10.
SIG13 = 0. 
SIG12 = -.5
SIG22 = 0.2
SIG14 = 0.5
SIG23 = -0.5
SIG24 = 0.1
SIG34 = -.5
SIG44 = 0.25

#decoupled case
SIG11 = 10.
SIG33 = 15.
SIG13 = 0. 
SIG12 = 0.
SIG22 = 0.2
SIG14 = 0.
SIG23 = 0.
SIG24 = 0.
SIG34 = -.5
SIG44 = 0.25

# Case |T|>0., |S31|>0
SIG11 = 10.
SIG33 = 12.
SIG13 = 1.
SIG12 = -.5
SIG22 = 0.2
SIG14 = -0.5
SIG23 = 0.9
SIG24 = 0.1
SIG34 = -0.7
SIG44 = 0.23

# Case |T|>0., S31=0
SIG11 = 10
SIG33 = 12
SIG13 = 0. 
SIG12 = -.5
SIG22 = 0.2
SIG14 = -0.5
SIG23 = 0.2
SIG24 = 0.1
SIG34 = 0.
SIG44 = 0.2



# Propagate by:
DS = -4
SIG11_DS, SIG12_DS, SIG13_DS, SIG14_DS,\
    SIG22_DS, SIG23_DS, SIG24_DS, SIG33_DS,\
    SIG34_DS, SIG44_DS = psm.propagate_full_Sigma_matrix_in_drift(\
                     SIG11, SIG12, SIG13, SIG14,
                     SIG22, SIG23, SIG24, SIG33,
                     SIG34, SIG44, DS)



# Generate S vector for test
S = np.linspace(-5-DS, 5-DS, 21)
#~ S = np.linspace(-5-DS, 5-DS, 201)


# Propagate Sigma matrix
Sigmas_at_0 = psm.Sigmas(SIG11_DS, SIG12_DS, SIG13_DS, SIG14_DS,
                     SIG22_DS, SIG23_DS, SIG24_DS, SIG33_DS,
                     SIG34_DS, SIG44_DS)



    
if mode == 'check_singularities':
    Sig_11_hat, Sig_33_hat, costheta, sintheta, \
        dS_Sig_11_hat, dS_Sig_33_hat, dS_costheta, dS_sintheta,\
        extra_data = psm.propagate_Sigma_matrix_vectorized(Sigmas_at_0, S, handle_singularities=False)
        
    Sig_11_hat_s, Sig_33_hat_s, costheta_s, sintheta_s, \
    dS_Sig_11_hat_s, dS_Sig_33_hat_s, dS_costheta_s, dS_sintheta_s,\
    extra_data_s = psm.propagate_Sigma_matrix_vectorized(Sigmas_at_0, S, handle_singularities=True)
        
elif mode == 'vs_sixtrack':
    import sigmas_sixtrak as sigs
    bcu = np.array([\
        Sigmas_at_0.Sig_11_0,
        Sigmas_at_0.Sig_33_0,
        Sigmas_at_0.Sig_13_0,
        Sigmas_at_0.Sig_12_0,
        Sigmas_at_0.Sig_14_0,
        Sigmas_at_0.Sig_22_0,
        Sigmas_at_0.Sig_23_0,
        Sigmas_at_0.Sig_24_0,
        Sigmas_at_0.Sig_34_0,
        Sigmas_at_0.Sig_44_0])    
    propagate_sigma_sixtrack_vectorized = np.vectorize(sigs.propagate_sigma, excluded=['bcu'])
    sintheta_s, costheta_s, Sig_11_hat_s, Sig_33_hat_s, dS_costheta_s, dS_sintheta_s,\
    dS_Sig_11_hat_s, dS_Sig_33_hat_s, Sig_11_sixt, Sig_33_sixt, Sig_13_sixt = propagate_sigma_sixtrack_vectorized(sp=S,bcu=bcu)
   
    Sig_11_hat, Sig_33_hat, costheta, sintheta, \
        dS_Sig_11_hat, dS_Sig_33_hat, dS_costheta, dS_sintheta,\
        extra_data = psm.propagate_Sigma_matrix_vectorized(Sigmas_at_0, S, handle_singularities=True)
        
else:
    raise ValueError('Mode not valid!!!')

# Plot results of the tests
pl.close('all')
fontsz = 14
lw = 3
mks = 10
ms.mystyle(fontsz = fontsz)


pl.close('all')




fig2 = pl.figure(2, figsize = (8*1.8, 6)); pl.clf()
fig2.set_facecolor('w')
sp0 = pl.subplot(2,2,1)
pl.plot(S, Sig_11_hat, '-b', label = 'Sig_11_hat', lw=lw, markersize=mks)
pl.plot(S, Sig_33_hat, '-r', label = 'Sig_33_hat', lw=lw, markersize=mks)
pl.plot(S, Sig_11_hat_s, '.--b', lw=lw/2, markersize=mks)
pl.plot(S, Sig_33_hat_s, '.--r', lw=lw/2, markersize=mks)
ms.sciy()
pl.grid('on')
pl.legend(loc='best', prop={'size':fontsz})
pl.subplot(2,2,3, sharex=sp0)
pl.plot(S, costheta, '-b', label='costheta', lw=lw, markersize=mks)
pl.plot(S, sintheta, '-r', label='sintheta', lw=lw, markersize=mks)
pl.plot(S, costheta_s, '.--b', lw=lw/2, markersize=mks)
pl.plot(S, sintheta_s, '.--r', lw=lw/2, markersize=mks)
pl.legend(loc='best', prop={'size':fontsz})
pl.ylim(-1.1, 1.1)
pl.xlabel('s [m]')
pl.grid('on')
fig2.subplots_adjust(top=.82, left=.07, right=.92)


pl.subplot(2,2,2, sharex=sp0)
pl.plot(S, dS_Sig_11_hat, '-b', label = 'dS_Sig_11_hat', lw=lw, markersize=mks)
pl.plot(S, dS_Sig_33_hat, '-r', label = 'dS_Sig_33_hat', lw=lw, markersize=mks)
pl.plot(S, dS_Sig_11_hat_s, '.--b', lw=lw/2, markersize=mks)
pl.plot(S, dS_Sig_33_hat_s, '.--r', lw=lw/2, markersize=mks)
ms.sciy()
pl.grid('on')
pl.legend(loc='best', prop={'size':fontsz})
pl.subplot(2,2,4, sharex=sp0)
pl.plot(S, dS_costheta, '-b', label='dS_costheta', lw=lw, markersize=mks)
pl.plot(S, dS_sintheta, '-r', label='dS_sintheta', lw=lw, markersize=mks)
pl.plot(S, dS_costheta_s, '.--b', lw=lw/2, markersize=mks)
pl.plot(S, dS_sintheta_s, '.--r', lw=lw/2, markersize=mks)
pl.legend(loc='best', prop={'size':fontsz})
pl.xlabel('s [m]')
pl.grid('on')

# Prepare title
R = SIG11-SIG33; W = SIG11+SIG33;T = R*R+4*SIG13*SIG13
a = SIG12-SIG34;b = SIG22-SIG44;c = SIG14+SIG23;d = SIG24
tit_str = 'Mode: %s At s=%.1e:\nSIG13=%.1e T=%.1e, a=%.1e, b=%.1e, c=%.1e, d=%.1e'%(mode,-DS, SIG13, T,a,b,c,d)
tit_str = tit_str.replace('e+00', '')
fname = tit_str.replace(' ', '_').replace('+', '').replace('-', 'm').replace('=', '').replace('\n', '__').replace(',','').replace(':','')
for fig in [fig2]: 
    fig.suptitle(tit_str)
    fig.subplots_adjust(top=.85)

if save_figs:
    fig2.savefig(fname+'.png', dpi=200)


# some investigations
enable_extra_plots = False

if enable_extra_plots:
#~ # Extract extra data
    cos2theta = []
    T = []
    R = []
    Sig11_calc = []
    Sig13_calc = []
    for ii in xrange(len(S)):
        cos2theta.append(extra_data_s[ii]['cos2theta'])
        T.append(extra_data_s[ii]['T'])
        R.append(extra_data_s[ii]['R'])
        Sig11_calc.append(extra_data_s[ii]['Sig_11'])
        Sig13_calc.append(extra_data_s[ii]['Sig_13'])
        
    T = np.array(T)
    R = np.array(R)

    a = SIG12-SIG34
    b = SIG22-SIG44
    c = SIG14+SIG23
    d = SIG24
    ddSS = S + DS

    pl.figure(1001)
    pl.plot(S, T)
    T1 = 4*ddSS**2*(a**2+c**2+ddSS*(a*b+2*c*d))
    pl.plot(S, T1)

    pl.figure(1002)
    pl.plot(S, R)
    pl.plot(S, 2*a*ddSS+b*ddSS**2)

    pl.figure(1000)
    pl.plot(S, cos2theta)
    pl.plot(S, np.abs(2*a+b*ddSS)/(np.sqrt((2*a+b*ddSS)**2 +4*(c+d*ddSS)**2)))
    
    pl.figure(1003)
    pl.plot(S, Sig13_calc, label='Sig13_calc')
    pl.legend(loc='best')

pl.show()


