#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#define COMPEARTH_PRIVATE_WRAP360 1
#include "compearth.h"

/*!
 * @brief Converts a fault normal vector to strike/dip.
 *
 * @param[in] n      Number of vectors.
 * @param[in] Xin    [3 x n] array of fault normal vectors where the first
 *                   (x,y,z) triplets are given by indices
 *                   0, 1, 2, and the second vector 3, 4, 5...
 *
 * @param[out] Xout  [2 x n] array of strikes and dips (degrees) where the
 *                   first fault's strikes and dips are given by
 *                   indices 0, 1, and the second by 2, 3...
 *
 * @author Carl Tape and translated to C by Ben Baker.
 *
 * @copyright MIT
 *
 */
void compearth_normal2strdip(const int n,
                             const double *__restrict__ Xin,
                             double *__restrict__ Xout)
{
    double n1th[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double n1ph[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double rho[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double x[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double y[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double z[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double kap1, theta1;
    int i, j, nmtLoc;
    const double deg = 180.0/M_PI;
    if (n < 1 || Xin == NULL || Xout == NULL)
    {
        if (n < 1){fprintf(stderr, "%s: No moment tensors\n", __func__);}
        if (Xin == NULL){fprintf(stderr, "%s: Xin is NULL\n", __func__);}
        if (Xout == NULL){fprintf(stderr, "%s: Xout is NULL\n", __func__);}
    }
    // Loop on chunks so I don't have to allocate scratch memory
    for (j=0; j<n; j=j+CE_CHUNKSIZE)
    {
        nmtLoc = MIN(CE_CHUNKSIZE, n - j);
        // Unpack the points
        for (i=0; i<nmtLoc; i++)
        {
            x[i] = Xin[3*(j+i)+0];
            y[i] = Xin[3*(j+i)+1];
            z[i] = Xin[3*(j+i)+2];
        }
        // dip and strike
        // dip: theta is [0,90] because we `force' the normal to point up
        // (It would be prudent to ensure that the third entry of the normal
        // vectors are always positive in this south-east-up basis.)
        // strike: kappa is unrestricted but we force it to [0,360]
        compearth_xyz2tp(nmtLoc, x, y, z, n1th, n1ph, rho);
        for (i=0; i<n; i++)
        {
            theta1 = n1th[i]*deg;
            kap1 = 90.0 - n1ph[i]*deg;
            kap1 = wrap360(kap1);
            Xout[2*(j+i)+0] = kap1;
            Xout[2*(j+i)+1] = theta1;
        }
    } // Loop on chunks
    return;
}

