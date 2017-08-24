#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "compearth.h"
/*!
 * @brief Converts (x,y,z) to spherical coordiantes (r,theta,phi)
 *
 * @param[in] n      Number of points in arrays.
 * @param[in] x      x cartesian ordinates. This is an array of dimension [n].
 * @param[in] y      y cartesian ordinates. This is an array of dimension [n].
 * @param[in] z      z cartesian ordinates. This is an array of dimension [n].
 *
 * @param[in] theta  Colatitude (radians). This is an array of dimension [n].
 * @param[in] phi    Longitude (radians). This is an array of dimension [n].
 * @param[in] r      Radius. This is an array of dimension [n].
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
