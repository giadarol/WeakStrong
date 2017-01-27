#include "weak_strong_4d_c.h"

#include <complex.h>


void weak_strong_4d(int a)
{

double complex z, res;
z = 3.+2.*I;

res = Faddeeva_w(z, 0.);

printf("Hello World\n");
printf("W(%.3e+j%.3e=%.3e+j%.3e\n", creal(z), cimag(z), creal(res), cimag(res));
printf("It should be: 0.0927108 + 0.128317i\n");



}
