#ifndef _TRANS_FIELD_GUASS_ROUND_
#define _TRANS_FIELD_GUASS_ROUND_

#include "constants.h"
#include <math.h>

typedef struct{
    double sigma;
    double Delta_x;
    double Delta_y;
}transv_field_gauss_round_data;

void get_transv_field_gauss_round(transv_field_gauss_round_data* data, double x, double y, double* Ex, double* Ey);

#endif
