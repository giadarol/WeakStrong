#ifndef _CMPX
#define _CMPX

#include <math.h>

typedef struct {
    double r;
    double i;
    } cmpx;

cmpx makecmpx(double r, double i);
cmpx cadd(cmpx a, cmpx b);
cmpx csub(cmpx a, cmpx b);
cmpx cmul(cmpx a, cmpx b);
cmpx cdiv(cmpx a, cmpx b);
cmpx expc(cmpx a);
cmpx cscale(cmpx a, double f);


#endif
