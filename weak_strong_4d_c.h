#ifndef __WEAKSTRONG4D
#define __WEAKSTRONG4D

#include <stdio.h>

#include <math.h>
#include <stdlib.h>

#include "particle.h"

#include "transverse_field_gauss_ellip.h"

typedef struct{
    double Delta_x;
    double Delta_y;
    double sigmax;
    double sigmay;
    double D_px_over_Ex;
    double D_py_over_Ey;
}weak_strong_4d_config;

void weak_strong_4d(particle*, weak_strong_4d_config*);



#endif
