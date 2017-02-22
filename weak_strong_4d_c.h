#ifndef __WEAKSTRONG4D
#define __WEAKSTRONG4D

#include <stdio.h>

#include <math.h>
#include <stdlib.h>
#include "cmpx.h"

// To use MIT Faddeeva library uncomment following
//#include <complex.h>
//#include "Faddeeva.h" 

// To use CERNLIB Faddeeva 
#include "faddeeva_cern.h"

// End Faddeeva choice 

#include "particle.h"

typedef struct{
    double sigmax;
    double sigmay;
    double D_px_over_Ex;
    double D_py_over_Ey;
}weak_strong_4d_config;

void weak_strong_4d(particle*, weak_strong_4d_config*);



#endif
