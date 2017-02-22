#ifndef _TRANS_FIELD_GUASS_ELLIP_
#define _TRANS_FIELD_GUASS_ELLIP_

#include "constants.h"
#include <math.h>
#include "cmpx.h"

typedef struct{
    double sigma_x;
    double sigma_y;
    double Delta_x;
    double Delta_y;
}transv_field_gauss_ellip_data;

// To use MIT Faddeeva library uncomment following
//#include <complex.h>
//#include "Faddeeva.h" 

// To use CERNLIB Faddeeva 
#include "faddeeva_cern.h"

#define EPSILON_0 (8.854187817620e-12)
#define SQRT_PI (1.772453850905515881919427556568)

void get_transv_field_gauss_ellip(transv_field_gauss_ellip_data* data, double x, double y, double* Ex_out, double* Ey_out);

#endif
