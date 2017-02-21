#ifndef __WEAKSTRONG4D
#define __WEAKSTRONG4D

#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include "cmpx.h"


#include "Faddeeva.h" 
// re-define differently if use another implementation of the w function

#include "particle.h"

typedef struct{
    double sigmax;
    double sigmay;
    double D_px_over_Ex;
    double D_py_over_Ey;
}weak_strong_4d_config;

void weak_strong_4d(particle*, weak_strong_4d_config*);



#endif
