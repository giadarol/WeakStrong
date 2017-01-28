#include "weak_strong_4d_c.h"


int main()
{
printf("In main!\n");
particle p;
weak_strong_4d_config conf;

p.x = 1.;
p.y = 2.;

conf.sigmax = 0.2;
conf.sigmay = 0.1;
conf.D_px_over_Ex = 1.;
conf.D_py_over_Ey = 1.;

weak_strong_4d(&p, &conf);
}

