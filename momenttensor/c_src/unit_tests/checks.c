#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "compearth.h"

int check_Udetcheck(void);
int check_lam2lune(void);

int main( )
{
    int ierr;
    ierr = check_Udetcheck();
    if (ierr != 0){printf("failed Udetcheck\n"); return EXIT_FAILURE;}
    printf("udetcheck was successful\n");

    ierr = check_lam2lune(); 
    if (ierr != 0){printf("failed lam2lune\n"); return EXIT_FAILURE;}
    printf("lam2lune was succesful\n");
    return EXIT_SUCCESS;
}

int check_Udetcheck(void)
{
    const char *fcnm = "check_Udetcheck\0";
    const double U[9] = {0.6948,   0.3171,   0.9502,  // col 1
                         0.0344,   0.4387,   0.3816,  // col 2
                         0.7655,   0.7952,   0.1869}; // col 3
    const double Ur[9] = {0.6948,  0.3171,   0.9502,  // col 1
                         -0.0344, -0.4387,  -0.3816,  // col 2 negated
                          0.7655,  0.7952,   0.1869};
    double Un[9], Un2[9];
    int i, ierr;
    ierr = compearth_Udetcheck(1, U, Un);
    if (ierr != 0)
    {
        printf("%s: Error in Udetcheck1\n", fcnm);
        return EXIT_FAILURE;
    }
    ierr = compearth_Udetcheck(1, Un, Un2);
    if (ierr != 0)
    {   
        printf("%s: Error in Udetcheck2\n", fcnm);
        return EXIT_FAILURE;
    }
    for (i=0; i<9; i++)
    {
        if (fabs(Ur[i] - Un[i]) > 1.e-12)
        {
            printf("%s: Udetcheck failed; %f %f\n", fcnm, Ur[i], Un[i]);
            return EXIT_FAILURE;
        }
        if (fabs(Un[i] - Un2[i]) > 1.e-12)
        {
            printf("%s: Udetcheck2 failed; %f %f\n", fcnm, Un[i], Un2[i]);
            return EXIT_FAILURE;
        }
    } 
    return EXIT_SUCCESS;
}
//============================================================================//
int check_CMT2TT(void)
{
    const double M[6] = {3.1083, 3.0444, 3.3823, -4.8550, -1.9493, 1.1105};

    return EXIT_SUCCESS;
}
//============================================================================//
int check_lam2lune(void)
{
    const char *fcnm = "check_lam2lune\0";
    const double lamRef1[6] = {
        -7.9621992143e+15,-8.2644837610e+15,-8.2644837610e+15,
        -7.5915701171e+15,-8.4370629711e+15,-8.4370629711e+15};
    double *gamma0, *gamma, *delta, *delta0, *lam, *M0, *M00;
    double *thetadc, *lamdev, *lamiso;
    int i, ib, ierr, ig, indx;
    const int ng = 100;
    const int nb = 100;
    int nmt = ng*nb;
    const double dg = (30.0 - -30.0)/(double) (ng - 1);
    const double db = (89.0 - -89.0)/(double) (nb - 1);
    gamma0 = (double *) calloc((size_t) nmt, sizeof(double));
    delta0 = (double *) calloc((size_t) nmt, sizeof(double));
    M00    = (double *) calloc((size_t) nmt, sizeof(double));
    gamma  = (double *) calloc((size_t) nmt, sizeof(double));
    delta  = (double *) calloc((size_t) nmt, sizeof(double));
    M0     = (double *) calloc((size_t) nmt, sizeof(double)); 
    lam  = (double *) calloc((size_t) (3*nmt), sizeof(double)); 
    for (ig=0; ig<ng; ig++)
    {
        for (ib=0; ib<nb; ib++)
        {
            indx = ig*nb + ib;
            gamma0[indx] =-30.0 + (double) ig*dg;
            delta0[indx] =-89.0 + (double) ib*db; // use latitude instead of colat
            M00[indx] = 1.e16;
        }
    }
    compearth_lune2lam(nmt, gamma0, delta0, M00, lam);
    // Verify the first two computations are right
    for (i=0; i<6; i++)
    {
        if (fabs(lam[i] - lamRef1[i])/1.e16 > 1.e-10)
        {
            printf("%s: error in lune2lam: %.10e %.10e\n",
                   fcnm, lam[i], lamRef1[i]);
            return EXIT_FAILURE;
        }
    }
    // Now perform the inverse operation
    thetadc = NULL;
    lamdev = NULL;
    lamiso = NULL;
    ierr = compearth_lam2lune(nmt, lam, 
                              gamma, delta, M0,
                              thetadc, lamdev, lamiso);
    if (ierr != 0)
    {
        printf("%s: error calling lam2lune\n", fcnm);
        return EXIT_FAILURE;
    }
    // See how close i was
    for (i=0; i<nmt; i++)
    {
        if (fabs(gamma[i] - gamma0[i]) > 1.e-12)
        {
            printf("%s: gamma different %e %e\n", fcnm, gamma[i], gamma0[i]);
        }
        if (fabs(delta[i] - delta0[i]) > 1.e-12)
        {
            printf("%s: delta different %e %e\n", fcnm, delta[i], delta[i]);
        }
        if (fabs(M0[i] - M00[i])/1.e15 > 1.e-12)
        {
           printf("%s: M0 different %e %e\n", fcnm, M0[i], M00[i]);
        }
    }
    free(gamma0);
    free(delta0);
    free(M00);
    free(gamma);
    free(delta);
    free(M0);
    free(lam);
    return EXIT_SUCCESS;
}
