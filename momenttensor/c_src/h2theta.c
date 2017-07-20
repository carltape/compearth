#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"
/*!
 * @brief Computes dip angle from point h
 *
 * @param[in] n       number of points in arrays
 * @param[in] h       h in rectilinear space s.t. \f$ h \in [0,1] \f$ [n]
 *
 * @param[out] theta  dip angle \f$ \theta \in [0, \pi/2] \f$ [n]
 *
 * @copyright MIT
 *
 */
void compearth_h2theta(const int n, const double *__restrict__ h,
                       double *__restrict__ theta)
{
    int i;
    #pragma omp simd
    for (i=0; i<n; i++)
    {
        theta[i] = acos(h[i]);
        theta[i] = fmax(0.0, fmin(theta[i], M_PI_2));
    }
    return;
}
