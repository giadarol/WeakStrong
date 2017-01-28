#include "weak_strong_4d_c.h"



#define EPSILON_0 (8.854187817620e-12)
#define SQRT_PI (1.772453850905515881919427556568)



void weak_strong_4d(particle* p, weak_strong_4d_config* conf)
{

    double sigmax = conf->sigmax;
    double sigmay = conf->sigmay;
    
    double x = p->x;
    double y = p->y;
    
    
    // I always go to the first quadrant and then apply the signs a posteriori
    // not necessary I think... to be simplified:wq

    double abx = abs(x);
    double aby = abs(y);
    
    double S, factBE, Ex, Ey;
    double complex etaBE, zetaBE, val;
    
    if (sigmax>sigmay){
        
        S=sqrt(2.*(sigmax*sigmax-sigmay*sigmay));
        factBE=1./(2.*EPSILON_0*SQRT_PI*S);
        etaBE=sigmay/sigmax*x+1j*sigmax/sigmay*y;
        zetaBE=x+I*y;
        
        val=factBE*(wfun(zetaBE/S)-cexp( -x*x/(2*sigmax*sigmax)-y*y/(2*sigmay*sigmay))*wfun(etaBE/S) );
           
        Ex=abs(cimag(val));
        Ey=abs(creal(val));
    }
    else{

        S=sqrt(2.*(sigmay*sigmay-sigmax*sigmax));
        factBE=1./(2.*EPSILON_0*SQRT_PI*S);
        etaBE=sigmax/sigmay*y+1j*sigmay/sigmax*x;
        zetaBE=y+I*x;
        
        val=factBE*(wfun(zetaBE/S)-cexp( -y*y/(2*sigmay*sigmay)-x*x/(2*sigmax*sigmax))*wfun(etaBE/S) );
           
        Ey=cimag(val);
        Ex=creal(val);
    }

    if(x<0) Ex=-Ex;
    if(y<0) Ey=-Ey;
    
    p->px = p->px + Ex*(conf->D_px_over_Ex);
    p->py = p->py + Ey*(conf->D_py_over_Ey);

/*
// Some basic testing:
double complex z, res;
z = 3.+2.*I;
res = Faddeeva_w(z, 0.);
printf("Hello World\n");
printf("W(%.3e+j%.3e=%.3e+j%.3e\n", creal(z), cimag(z), creal(res), cimag(res));
printf("It should be: 0.0927108 + 0.128317i\n");
*/


}
