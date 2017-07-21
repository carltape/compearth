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

int check_CMT2TT(void)
{
    const double M[6] = {3.1083, 3.0444, 3.3823, -4.8550, -1.9493, 1.1105};

    return EXIT_SUCCESS;
}

int check_lam2lune(void)
{
    double *gvec, *dvec, *lam, *M00;
    int ib, ig, indx;
    const int ng = 100;
    const int nb = 100;
    const int nmt = ng*nb;
    const double dg = (30.0 - -30.0)/(double) (ng - 1);
    const double db = (89.0 - -89.0)/(double) (nb - 1);
    gvec = (double *) calloc((size_t) nmt, sizeof(double));
    dvec = (double *) calloc((size_t) nmt, sizeof(double));
    M00  = (double *) calloc((size_t) nmt, sizeof(double));
    lam  = (double *) calloc((size_t) (3*nmt), sizeof(double)); 
    for (ig=0; ig<ng; ig++)
    {
        for (ib=0; ib<nb; ib++)
        {
            indx = ig*nb + ib;
            gvec[indx] =-30.0 + (double) ig*dg;
            dvec[indx] =-89.0 + (double) ib*db; // use latitude instead of colat
            M00[indx] = 1.e16;
        }
    }
    compearth_lune2lam(nmt, gvec, dvec, M00, lam); 
    
    free(gvec);
    free(dvec);
    free(M00);
    free(lam);
    return EXIT_SUCCESS;
}
