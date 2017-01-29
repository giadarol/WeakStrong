#include "weak_strong_4d_c.h"


#define EPSILON_0 (8.854187817620e-12)
#define SQRT_PI (1.772453850905515881919427556568)



void weak_strong_4d(particle* p, weak_strong_4d_config* conf)
{

    double sigmax = conf->sigmax;
    double sigmay = conf->sigmay;
    
    double xx = p->x;
    double yy = p->y;
    
    
    // I always go to the first quadrant and then apply the signs a posteriori
    // numerically more stable (see http://inspirehep.net/record/316705/files/slac-pub-5582.pdf)

    double abx = fabs(xx);
    double aby = fabs(yy);
    
    //printf("x = %.2e y = %.2e abx = %.2e aby = %.2e", xx, yy, abx, aby);
    
    
    double S, factBE, Ex, Ey;
    double complex etaBE, zetaBE, val;
    
    
    
    if (sigmax>sigmay){
        
        S=sqrt(2.*(sigmax*sigmax-sigmay*sigmay));
        factBE=1./(2.*EPSILON_0*SQRT_PI*S);
        etaBE=sigmay/sigmax*abx+I*sigmax/sigmay*aby;
        zetaBE=abx+I*aby;
        val=factBE*(wfun(zetaBE/S)-cexp( -abx*abx/(2*sigmax*sigmax)-aby*aby/(2*sigmay*sigmay))*wfun(etaBE/S) );
           
        Ex=cimag(val);
        Ey=creal(val);
    }
    else if (sigmax<sigmay){

        S=sqrt(2.*(sigmay*sigmay-sigmax*sigmax));
        factBE=1./(2.*EPSILON_0*SQRT_PI*S);
        etaBE=sigmax/sigmay*aby+I*sigmay/sigmax*abx;
        zetaBE=aby+I*abx;
        
        val=factBE*(wfun(zetaBE/S)-cexp( -aby*aby/(2*sigmay*sigmay)-abx*abx/(2*sigmax*sigmax))*wfun(etaBE/S) );
           
        Ey=cimag(val);
        Ex=creal(val);
    }
    else{
        printf("Round beam not implemented!\n");
        exit(1);
    }

    if(xx<0) Ex=-Ex;
    if(yy<0) Ey=-Ey;
    
    //printf("Ex = %.2e Ey = %.2e \n", Ex, Ey);
    
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
