#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

int compearth_lam2phizeta(const int n, const double *__restrict__ lam,
                          double *__restrict__ phi,
                          double *__restrict__ zeta)
{
    const char *fcnm = "compearth_lam2phizeta\0";
    double *lamSort, lam1, lam2, lam3, rho;
    int i, ierr;
    const double sqrt2 = sqrt(2.0);
    const double deg = 180.0/M_PI;
    lamSort = (double *) calloc((size_t) (3*n), sizeof(double));
    ierr = compearth_lamsort(n, lam, lamSort); 
    if (ierr != 0)
    {
        printf("%s: Error sorting eigentriples\n", fcnm);
        free(lamSort);
        return -1;
    }
    for (i=0; i<n; i++)
    {
        lam1 = lamSort[3*i];
        lam2 = lamSort[3*i+1];
        lam3 = lamSort[3*i+2];
        rho = sqrt(lam1*lam1 + lam2*lam2 + lam3*lam3);
        // TT 2013 Eqns 24
        phi[i] = atan2(lam1 - 2.0*lam2 + lam3,
                       sqrt2*(lam1 + lam2 + lam3))*deg;
        zeta[i] = acos( sqrt(2.0*(lam1 - lam2)*(lam2 - lam3))/rho )*deg;
    }
    free(lamSort);
    return 0;
}
