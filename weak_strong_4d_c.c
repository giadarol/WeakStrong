#include "weak_strong_4d_c.h"

#include <complex.h>

#define EPSILON_0 (8.854187817620e-12)


void weak_strong_4d(particle* p, weak_strong_4d_config* conf)
{

    double sigmax = conf->sigmax;
    double sigmay = conf->sigmay;
    
    double x = p->x;
    double y = p->y;
    
    
    
    double abx = abs(x);
    double aby = abs(y);
    
    
    
    
    /*
    eps0=8.854187817620e-12;
    
    
    if sigmax>sigmay:
    
        S=sqrt(2*(sigmax*sigmax-sigmay*sigmay));
        factBE=1/(2*eps0*sqrt(pi)*S);
        etaBE=sigmay/sigmax*x+1j*sigmax/sigmay*y;
        zetaBE=x+1j*y;
        
        val=factBE*(wfun(zetaBE/S)-exp( -x*x/(2*sigmax*sigmax)-y*y/(2*sigmay*sigmay))*wfun(etaBE/S) );
           
        Ex=abs(val.imag)*sign(xin);
        Ey=abs(val.real)*sign(yin);
    
    else:
    
        S=sqrt(2*(sigmay*sigmay-sigmax*sigmax));
        factBE=1/(2*eps0*sqrt(pi)*S);
        etaBE=sigmax/sigmay*y+1j*sigmay/sigmax*x;
        yetaBE=y+1j*x;
        
        val=factBE*(wfun(yetaBE/S)-exp( -y*y/(2*sigmay*sigmay)-x*x/(2*sigmax*sigmax))*wfun(etaBE/S) );
           
        Ey=abs(val.imag)*sign(yin);
        Ex=abs(val.real)*sign(xin);
*/

// Some basic testing:
 
double complex z, res;
z = 3.+2.*I;
res = Faddeeva_w(z, 0.);
printf("Hello World\n");
printf("W(%.3e+j%.3e=%.3e+j%.3e\n", creal(z), cimag(z), creal(res), cimag(res));
printf("It should be: 0.0927108 + 0.128317i\n");



}
