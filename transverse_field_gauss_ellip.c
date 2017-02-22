#include "transverse_field_gauss_ellip.h"


// To use MIT Faddeeva library uncomment following
// cmpx wfun(cmpx zz){
//     double complex temp;
//     cmpx res;
//     temp = Faddeeva_w((zz.r + I*zz.i), 0.);
//     res = makecmpx(creal(temp), cimag(temp));
//     return res;
// }   


// To use CERNLIB Faddeeva 
 cmpx wfun(cmpx zz){
    cmpx res;
    cerrf(zz.r, zz.i , &(res.r), &(res.i));
    return res;
}

void get_transv_field_gauss_ellip(transv_field_gauss_ellip_data* data, double x, double y, double* Ex_out, double* Ey_out){
    
    double sigmax = data->sigma_x;
    double sigmay = data->sigma_y;
    

    
    // I always go to the first quadrant and then apply the signs a posteriori
    // numerically more stable (see http://inspirehep.net/record/316705/files/slac-pub-5582.pdf)

    double abx = fabs(x - data->Delta_x);
    double aby = fabs(y - data->Delta_y);
    
    //printf("x = %.2e y = %.2e abx = %.2e aby = %.2e", xx, yy, abx, aby);
    
    
    double S, factBE, Ex, Ey;
    cmpx etaBE, zetaBE, val;
    
    
    
    if (sigmax>sigmay){
        
        S=sqrt(2.*(sigmax*sigmax-sigmay*sigmay));
        factBE = 1./(2.*EPSILON_0*SQRT_PI*S);
        
        etaBE = makecmpx(sigmay/sigmax*abx, sigmax/sigmay*aby);
        zetaBE = makecmpx(abx, aby);
        
        val = csub( wfun(cscale(zetaBE, 1./S)),
                    cscale(wfun(cscale(etaBE,1./S)), exp( -abx*abx/(2*sigmax*sigmax)-aby*aby/(2*sigmay*sigmay))) );
         
           
        Ex=factBE*val.i;
        Ey=factBE*val.r;
    }
    else if (sigmax<sigmay){

        S=sqrt(2.*(sigmay*sigmay-sigmax*sigmax));
        factBE = 1./(2.*EPSILON_0*SQRT_PI*S);
        
        etaBE = makecmpx(sigmax/sigmay*aby, sigmay/sigmax*abx);
        zetaBE = makecmpx(aby, abx);
        
        val = csub( wfun(cscale(zetaBE, 1./S)),
                    cscale(wfun(cscale(etaBE,1./S)), exp( -aby*aby/(2*sigmay*sigmay)-abx*abx/(2*sigmax*sigmax))) );
                   
        Ey=factBE*val.i;
        Ex=factBE*val.r;
        
        
    }
    else{
        printf("Round beam not implemented!\n");
        exit(1);
    }

    if((x - data->Delta_x)<0) Ex=-Ex;
    if((y - data->Delta_y)<0) Ey=-Ey;
    
    (*Ex_out) = Ex;
    (*Ey_out) = Ey;
    
    


}
