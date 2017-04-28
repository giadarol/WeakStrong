#include "cmpx.h"
#include <stdio.h>

int main()
{
cmpx a, b, c;
//~ a.r = -3.;
//~ a.i = 5.;
a = makecmpx(-3.,5.);

//~ b.r = 7.;
//~ b.i = -6.;
b = makecmpx(7.,-6);

c = cadd(a,b);
printf("a+b = %f + I %f\n", c.r, c.i);
c = csub(a,b);
printf("a-b = %f + I %f\n", c.r, c.i);
c = cmul(a,b);
printf("a*b = %f + I %f\n", c.r, c.i);
c = cdiv(a,b);
printf("a/b = %f + I %f\n", c.r, c.i);
c = expc(a);
printf("e^a = %f + I %f\n", c.r, c.i);
c = scale(a, 5.);
printf("5*a = %f + I %f\n", c.r, c.i);
}
 
