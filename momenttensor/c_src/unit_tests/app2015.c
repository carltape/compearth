#include <stdio.h>
#include <stdlib.h>
#include "compearth.h"

int main()
{
    // These are the coordinates you start with in the regular space
    double u = 3.0*M_PI/8.0;
    double v =-1.0/9.0;
    double h = 3.0/4.0;
    double kappa = 4.0*M_PI/5.0;
    double sigma =-M_PI/2.0;
    double M0 = 1.0/sqrt(2.0);
    // This is some junk for obtaining beta from u 
    double tol = 1.e-12;
    int maxit = 20, ierr;
    const double pi180 = 180.0/M_PI;///180.0; 
    double M_use[6], M_use_ref[6], M_ned[6], M_temp[6], lam[3], U[9],
           angle, beta, delta, gamma, theta, xnorm, w;
    int i;
    // Map from regular space to stretched space
    w = 3.0*M_PI/8.0 - u;
    compearth_rect2lune(1, &v, 1, &w, &gamma, &delta);
    gamma = gamma*M_PI/180.0;
    beta = (90.0 - delta)*M_PI/180.0;
    printf("--------------------tape appendix 2015------------------\n");
    printf("(beta,gamma)=%12.10e %12.10e\n", beta, gamma);
    printf("--------------------------------------------------------\n");
    ierr = compearth_u2beta(1, maxit, 2, &u, tol, &beta);
    compearth_v2gamma(1, &v, &gamma);
    compearth_h2theta(1, &h, &theta);
    printf("--------------------tape appendix 2015------------------\n");
    printf("(beta,gamma)=%12.10e %12.10e\n", beta, gamma);
    printf("theta=%12.10e\n",theta);
    printf("--------------------------------------------------------\n");
    // Now compute the correponding moment tensor
    delta = M_PI/2.0 - beta;
    double gammaDeg = gamma*pi180;
    double deltaDeg = delta*pi180;
    double kappaDeg = kappa*pi180;
    double thetaDeg = theta*pi180;
    double sigmaDeg = sigma*pi180;
    ierr = compearth_TT2CMT(1, &gammaDeg, &deltaDeg, &M0,
                            &kappaDeg, &thetaDeg, &sigmaDeg,
                            M_use, lam, U);
    if (ierr != 0){
        printf("\n");
    }
    printf("M_use:\n");
    printf("%e\n%e\n%e\n%e\n%e\n%e\n", M_use[0], M_use[1], M_use[2],
                                       M_use[3], M_use[4], M_use[5]);
    for (i=0; i<6; i++){M_use_ref[i] = M_use[i];}
    printf("Lambda:\n");
    printf("%e\n%e\n%e\n", lam[0], lam[1], lam[2]);
    printf("U:\n");
    printf("%e %e %e\n", U[0], U[3], U[6]);
    printf("%e %e %e\n", U[1], U[4], U[7]);
    printf("%e %e %e\n", U[2], U[5], U[8]);
    // Finally convert from USE TO NED {xx, yy, zz, xy, xz, yz} for Herrmann's code
    ierr = compearth_convertMT(1, CE_USE, CE_NED, M_use, M_ned);
    printf("M_ned\n");
    printf("%e\n%e\n%e\n%e\n%e\n%e\n", M_ned[0], M_ned[1], M_ned[2],
                                       M_ned[3], M_ned[4], M_ned[5]);
    ierr = compearth_convertMT(1, CE_NED, CE_USE, M_ned, M_use);
    printf("M_use\n");
    for (i=0; i<6; i++){
        printf("%e %e\n", M_use[i], M_use_ref[i] - M_use[i]);
        M_temp[i] = M_ned[i];
    }

    M_temp[0] = M_temp[0]*0.5;
    M_temp[1] = M_temp[1]*0.25;
    M_temp[2] = M_temp[2]*0.2;
    M_temp[3] = M_temp[3]*0.1;
    ierr = compearth_normMT(1, M_temp, CE_TWO_NORM, 2.0, &xnorm);
    printf("xnorm: %f\n", xnorm);
    ierr = compearth_angleMT(1, M_temp, M_use, &angle);
    printf("angle: %f\n", angle);
    return 0;
}
