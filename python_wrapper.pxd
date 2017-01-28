cdef extern from "weak_strong_4d_c.h":
    
    ctypedef struct particle:
        double x;
        double px;
        double y;
        double py;
    
    cdef void weak_strong_4d(particle*);
