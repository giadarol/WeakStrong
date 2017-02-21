#include "weak_strong_4d_c.h"


#define EPSILON_0 (8.854187817620e-12)
#define SQRT_PI (1.772453850905515881919427556568)

cmpx wfun(cmpx zz){
    double complex temp;
    cmpx res;
    temp = Faddeeva_w((zz.r + I*zz.i), 0.);
    res = makecmpx(creal(temp), cimag(temp));
    return res;
}   



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
    cmpx etaBE, zetaBE, val, first_term, second_term;
    
    
    
    if (sigmax>sigmay){
        
        S=sqrt(2.*(sigmax*sigmax-sigmay*sigmay));
        factBE = 1./(2.*EPSILON_0*SQRT_PI*S);
        
        etaBE = makecmpx(sigmay/sigmax*abx, sigmax/sigmay*aby);
        zetaBE = makecmpx(abx, aby);
        
        first_term = wfun(cscale(zetaBE, 1./S));
        second_term = cscale(wfun(cscale(etaBE,1./S)), exp( -abx*abx/(2*sigmax*sigmax)-aby*aby/(2*sigmay*sigmay)));
        
        val = csub(first_term,second_term);
           
        Ex=factBE*val.i;
        Ey=factBE*val.r;
    }
    else if (sigmax<sigmay){

        S=sqrt(2.*(sigmay*sigmay-sigmax*sigmax));
        factBE = 1./(2.*EPSILON_0*SQRT_PI*S);
        
        etaBE = makecmpx(sigmax/sigmay*aby, sigmay/sigmax*abx);
        zetaBE = makecmpx(aby, abx);
        
        first_term = wfun(cscale(zetaBE, 1./S));
        second_term = cscale(wfun(cscale(etaBE,1./S)), exp( -aby*aby/(2*sigmay*sigmay)-abx*abx/(2*sigmax*sigmax)));
        
        val = csub(first_term,second_term);
           
        Ey=factBE*val.i;
        Ex=factBE*val.r;
        
        
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
