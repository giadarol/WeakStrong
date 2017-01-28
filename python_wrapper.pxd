cdef extern from "weak_strong_4d_c.h":
    
    ctypedef struct particle:
        double x;
        double px;
        double y;
        double py;
        
    ctypedef struct weak_strong_4d_config:
        double sigmax;
        double sigmay;
        double D_px_over_Ex;
        double D_py_over_Ey;

    
    cdef void weak_strong_4d(particle*, weak_strong_4d_config*);
