cimport python_wrapper

cpdef weak_strong_single_particle(part, sigmax, sigmay, D_px_over_Ex=1., D_py_over_Ey=1.):
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
    
    #act on C particle
    weak_strong_4d(&p, &conf);
    
    #C particle to python particle
    part.x = p.x
    part.y = p.y
    part.px = p.px
    part.py = p.py
