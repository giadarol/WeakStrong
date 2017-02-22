cimport python_wrapper

cpdef weak_strong_single_particle(part, sigmax, sigmay, D_px_over_Ex=1., D_py_over_Ey=1., Delta_x=0., Delta_y=0.):
    cdef particle p;
    cdef weak_strong_4d_config conf
    #python particle to c particle
    p.x = part.x
    p.y = part.y
    p.px = part.px
    p.py = part.py
    
    conf.sigmax = sigmax
    conf.sigmay = sigmay
    conf.D_px_over_Ex = D_px_over_Ex
    conf.D_py_over_Ey = D_py_over_Ey  
    conf.Delta_x = Delta_x
    conf.Delta_y = Delta_y 
    
    #act on C particle
    weak_strong_4d(&p, &conf);
    
    #C particle to python particle
    part.x = p.x
    part.y = p.y
    part.px = p.px
    part.py = p.py
    
cpdef test_efield_gauss_round(x, y, sigma, Delta_x, Delta_y):
    
    cdef transv_field_gauss_round_data data;
    cdef double Ex, Ey;
    
    data.sigma = sigma
    data.Delta_x = Delta_x
    data.Delta_y = Delta_y
    
    get_transv_field_gauss_round(&data, x, y, &Ex, &Ey);
    
    return Ex, Ey
