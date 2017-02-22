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
    
    
cdef extern from "transverse_field_gauss_round.h":
    
    ctypedef struct transv_field_gauss_round_data:
        double sigma;
        double Delta_x;
        double Delta_y;

    cdef void get_transv_field_gauss_round(transv_field_gauss_round_data* data, double x, double y, double* Ex, double* Ey);
