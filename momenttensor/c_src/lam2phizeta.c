#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

int compearth_lam2phizeta(const int n, const double *__restrict__ lam,
                          double *__restrict__ phi,
                          double *__restrict__ zeta)
{
    double lamSort[3*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double lam1, lam2, lam3, rho;
    int i, ierr, j, nloc;
    const double sqrt2 = sqrt(2.0);
    const double deg = 180.0/M_PI;
    for (i=0; i<n; i=i+CE_CHUNKSIZE)
    {
        nloc = MIN(CE_CHUNKSIZE, n - i); 
        ierr = compearth_lamsort(nloc, lam, lamSort);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error sorting eigentriples\n", __func__);
            return -1;
        }
        for (j=0; j<nloc; j++)
        {
            lam1 = lamSort[3*j];
            lam2 = lamSort[3*j+1];
            lam3 = lamSort[3*j+2];
            rho = sqrt(lam1*lam1 + lam2*lam2 + lam3*lam3);
            // TT 2013 Eqns 24
            phi[i+j] = atan2(lam1 - 2.0*lam2 + lam3,
                             sqrt2*(lam1 + lam2 + lam3))*deg;
            zeta[i+j] = acos( sqrt(2.0*(lam1 - lam2)*(lam2 - lam3))/rho )*deg;
        }
    }
    return 0;
}
