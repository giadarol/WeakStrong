#include "weak_strong_4d_c.h"


void weak_strong_4d(particle* p, weak_strong_4d_config* conf)
{

    transv_field_gauss_ellip_data data;
    
    data.sigma_x = conf->sigmax,
    data.sigma_y = conf->sigmay;
    data.Delta_x = conf->Delta_x;
    data.Delta_y = conf->Delta_y;
    
    double xx = p->x;
    double yy = p->y;
    double Ex, Ey;
    
    get_transv_field_gauss_ellip(&data, xx, yy, &Ex, &Ey);
    
    
    
    p->px = p->px + Ex*(conf->D_px_over_Ex);
    p->py = p->py + Ey*(conf->D_py_over_Ey);

}
