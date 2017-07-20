#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "compearth.h"

static double wrap360(const double lon);
/*!
 * @brief Converts a fault normal vector to strike/dip
 *
 * @param[in] n      number of vectors
 * @param[in] Xin    [3 x n] fault normal vectors where the first
 *                   vectors (x,y,z) triplets are given by indices
 *                   0, 1, 2, and the second vector 3, 4, 5...
 *
 * @param[out] Xout  [2 x n] strikes and dips (degrees) where the
 *                   first fault's strikes and dips are given by
 *                   indices 0, 1, and the second by 2, 3...
 *
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @copyright MIT
 *
 */
void compearth_normal2strdip(const int n,
                             const double *__restrict__ Xin,
                             double *__restrict__ Xout)
{
    double *x, *y, *z, *n1th, *n1ph, *rho;
    double n1th64[64], n2ph64[64], rho64[64], x64[64], y64[64], z64[64],
           kap1, theta1;
    int i;
    const double deg = 180.0/M_PI;
    // Set space for (x,y,z) triplets 
    if (n <= 64)
    {
        x = x64;
        y = y64;
        z = z64;
        n1th = n1th64;
        n1th = n1th64;
        rho = rho64;
    }
    else
    {
        x = (double *) calloc((size_t) n, sizeof(double));
        y = (double *) calloc((size_t) n, sizeof(double));
        z = (double *) calloc((size_t) n, sizeof(double));
        n1th = (double *) calloc((size_t) n, sizeof(double));
        n1ph = (double *) calloc((size_t) n, sizeof(double));
        rho = (double *) calloc((size_t) n, sizeof(double));
    }
    // Unpack the points 
    for (i=0; i<n; i++)
    {
        x[i] = Xin[3*i+0];
        y[i] = Xin[3*i+1]; 
        z[i] = Xin[3*i+2];
    }
    // dip and strike
    // dip: theta is [0,90] because we `force' the normal to point up
    // (It would be prudent to ensure that the third entry of the normal
    // vectors are always positive in this south-east-up basis.)
    // strike: kappa is unrestricted but we force it to [0,360]
    compearth_xyz2tp(n, x, y, z, n1th, n1ph, rho);
    for (i=0; i<n; i++)
    {
        theta1 = n1th[i]*deg;
        kap1 = 90.0 - n1ph[i]*deg;
        kap1 = wrap360(kap1);
        Xout[2*i+0] = kap1;
        Xout[2*i+1] = theta1;
    }
    if (n > 64)
    {
        free(x);
        free(y);
        free(z);
        free(n1th);
        free(n1ph);
        free(rho);
    }
    x = NULL;
    y = NULL;
    z = NULL;
    n1th = NULL;
    n1ph = NULL;
    rho = NULL;
    return;
}

/*!
 * @brief Wraps angle in degrees to [0,360]
 *
 * @param[in] lon   angle to wrap (degrees)
 *
 * @result wrapped angle in range [0,360]
 *
 * @author Carl Tape translated to C by Ben Baker
 *
 * @copyright MIT
 *
 */
static double wrap360(const double lon)
{
    double lonw;
    bool lpos;
    lpos = false; 
    if (lon > 0.0){lpos = true;}
    lonw = fmod(lon, 360.0); 
    if (lonw == 0.0 && lpos){lonw = 360.0;}
    return lonw;
}
