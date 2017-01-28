cimport python_wrapper

cpdef run():
    cdef particle p;
    p.x = 10.;
    p.y = 20.;
    weak_strong_4d(&p);
