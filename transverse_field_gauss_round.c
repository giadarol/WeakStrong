#include "transverse_field_gauss_round.h"

void get_transv_norm_efield(transv_norm_efield_data* data, double x, double y, double* Ex, double* Ey){
    
    double r2, temp;
    
    r2 = (x-data->Delta_x)*(x-data->Delta_x)+(y-data->Delta_y)*(y-data->Delta_y);
    if (r2<1e-20) temp = sqrt(r2)/(2.*pi*epsilon_0*data->sigma);
    else          temp = (1-exp(-0.5*r2/(data->sigma*data->sigma)))/(2.*pi*epsilon_0*r2);
    
    (*Ex) = temp*( x-data->Delta_x );
    (*Ey) = temp*( y-data->Delta_y );
}
