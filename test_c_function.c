#include "weak_strong_4d_c.h"
#include "transverse_field_gauss_round.h"

int main()
{
printf("In main!\n");
particle p;
weak_strong_4d_config conf;

p.x = 1.;
p.y = 20.;
p.px = 0.;
p.py = 0.;

conf.sigmax = 0.2;
conf.sigmay = 0.25;
conf.D_px_over_Ex = 1.;
conf.D_py_over_Ey = 1.;

weak_strong_4d(&p, &conf);

printf("Ex=%.2e Ey=%.2e\n", p.px, p.py);


transv_norm_efield_data data;
double Ex, Ey;

data.sigma = 0.3;
data.Delta_x = 0.;
data.Delta_y = 0;

get_transv_norm_efield(&data, 0.4*data.sigma, 0.7*data.sigma, &Ex, &Ey);
    



}


