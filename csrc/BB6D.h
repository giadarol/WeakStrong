#ifndef _BB6D_
#define _BB6D_

#include "BB6D_boost.h"
#include "BB6D_propagateSigma.h"
#include "BB6D_transverse_fields.h"
#include "constants.h"

typedef struct{
    double q_part;
    BB6D_boost_data parboost;
    BB6D_Sigmas Sigmas_0_star;
    double min_sigma_diff;
    double threshold_singular;
    long int N_slices;
    double* N_part_per_slice;
    double* x_slices_star;
    double* y_slices_star;
    double* sigma_slices_star;
}BB6D_data;

void BB6D_track(double* x, double* px, double* y, double* py, double* sigma, 
                double* delta, double q0, double p0, BB6D_data *bb6ddata){

    int N_slices = (int)(bb6ddata->N_slices);
    int i_slice;
    
    double x_star = *x;
    double px_star = *px;
    double y_star = *y;
    double py_star = *py;              
    double sigma_star = *sigma;
    double delta_star = *delta ;  
    
    
    // Boost coordinates of the weak beam
    BB6D_boost(&(bb6ddata->parboost), &x_star, &px_star, &y_star, &py_star, 
                &sigma_star, &delta_star);
    
    // Synchro beam
    for (i_slice=0; i_slice<N_slices; i_slice++)
    {
        double sigma_slice_star = bb6ddata->sigma_slices_star[i_slice];
        double x_slice_star = bb6ddata->x_slices_star[i_slice];
        double y_slice_star = bb6ddata->y_slices_star[i_slice];
        
        //Compute force scaling factor
        double Ksl = bb6ddata->N_part_per_slice[i_slice]*bb6ddata->q_part*q0/(p0*C_LIGHT);

        //Identify the Collision Point (CP)
        double S = 0.5*(sigma_star - sigma_slice_star);
        
        // Propagate sigma matrix
        double Sig_11_hat_star, Sig_33_hat_star, costheta, sintheta;
        double dS_Sig_11_hat_star, dS_Sig_33_hat_star, dS_costheta, dS_sintheta;
        
        // Get strong beam shape at the CP
        BB6D_propagate_Sigma_matrix(&(bb6ddata->Sigmas_0_star),
            S, bb6ddata->threshold_singular, 1,
            &Sig_11_hat_star, &Sig_33_hat_star, 
            &costheta, &sintheta,
            &dS_Sig_11_hat_star, &dS_Sig_33_hat_star, 
            &dS_costheta, &dS_sintheta);
            
        // Evaluate transverse coordinates of the weake baem w.r.t. the strong beam centroid
        double x_bar_star = x_star + px_star*S - x_slice_star;
        double y_bar_star = y_star + py_star*S - y_slice_star;
        
        // Move to the uncoupled reference frame
        double x_bar_hat_star = x_bar_star*costheta +y_bar_star*sintheta;
        double y_bar_hat_star = -x_bar_star*sintheta +y_bar_star*costheta;
        
        // Compute derivatives of the transformation
        double dS_x_bar_hat_star = x_bar_star*dS_costheta +y_bar_star*dS_sintheta;
        double dS_y_bar_hat_star = -x_bar_star*dS_sintheta +y_bar_star*dS_costheta;
        
        // Get transverse fieds
        double Ex, Ey, Gx, Gy;
        get_Ex_Ey_Gx_Gy_gauss(x_bar_hat_star, y_bar_hat_star, 
            sqrt(Sig_11_hat_star), sqrt(Sig_33_hat_star), bb6ddata->min_sigma_diff,
            &Ex, &Ey, &Gx, &Gy);
            
        // Compute kicks
        double Fx_hat_star = Ksl*Ex;
        double Fy_hat_star = Ksl*Ey;
        double Gx_hat_star = Ksl*Gx;
        double Gy_hat_star = Ksl*Gy;
        
        // Move kisks to coupled reference frame
        double Fx_star = Fx_hat_star*costheta - Fy_hat_star*sintheta;
        double Fy_star = Fx_hat_star*sintheta + Fy_hat_star*costheta;
        
        // Compute longitudinal kick
        double Fz_star = 0.5*(Fx_hat_star*dS_x_bar_hat_star  + Fy_hat_star*dS_y_bar_hat_star+
                       Gx_hat_star*dS_Sig_11_hat_star + Gy_hat_star*dS_Sig_33_hat_star);
                       
        // Apply the kicks (Hirata's synchro-beam)
        delta_star = delta_star + Fz_star+0.5*(
                    Fx_star*(px_star+0.5*Fx_star)+
                    Fy_star*(py_star+0.5*Fy_star));
        x_star = x_star - S*Fx_star;
        px_star = px_star + Fx_star;
        y_star = y_star - S*Fy_star;
        py_star = py_star + Fy_star;
        

    }
    
    // Inverse boost on the coordinates of the weak beam
    BB6D_inv_boost(&(bb6ddata->parboost), &x_star, &px_star, &y_star, &py_star, 
                &sigma_star, &delta_star);
                
    *x = x_star;
    *px = px_star;
    *y = y_star;
    *py = py_star;
    *sigma = sigma_star;
    *delta = delta_star;
                    
}

#endif
