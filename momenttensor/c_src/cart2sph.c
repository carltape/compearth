#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "compearth.h"
/*!
 * @brief Converts (x,y,z) to spherical coordiantes (r,theta,phi)
 *
 * @param[in] n      number of points in arrays
 * @param[in] x      x cartesian ordinates [n]
 * @param[in] y      y cartesian ordinates [n]
 * @param[in] z      z cartesian ordinates [n]
 *
 * @param[in] theta  colatitude (radians) [n]
 * @param[in] phi    longitude (radians) [n] 
 * @param[in] r      radius [n]
 *
 * @author Ben Baker
 *
 * @copyright MIT
 *
 */
void compearth_matlab_cart2sph(const int n,
                               const double *__restrict__ x,
                               const double *__restrict__ y,
                               const double *__restrict__ z,
                               double *__restrict__ theta,
                               double *__restrict__ phi,
                               double *__restrict__ r)
{
    double x2, y2, z2;
    int i;
    #pragma omp simd
    for (i=0; i<n; i++)
    {
        x2 = x[i]*x[i];
        y2 = y[i]*y[i];
        z2 = z[i]*z[i];
        theta[i] = atan2(y[i], x[i]);
        phi[i] = atan2(z[i], sqrt(x2 + y2));
        r[i] = sqrt(x2 + y2 + z2);
    }
    return;
}
