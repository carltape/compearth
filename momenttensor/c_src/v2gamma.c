#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief Commputes lune longitude from point v
 *
 * @param[in] n       number of points in array
 * @param[in] v       v in rectilinear space s.t. \f$ v \in [-1/3, 1/3] \f$ [n]
 *
 * @param[out] gamma  lune longitude \f$ \gamma \in [-\pi/6, \pi/6 \f$  [n]
 *
 * @date 2016 - Ben Baker converted Carl Tape's v2gamma.m to C
 *
 * @copyright MIT
 *
 */
void compearth_v2gamma(const int n, const double *__restrict__ v,
                       double *__restrict__ gamma)
{
    int i;
    const double third = 1.0/3.0;
    #pragma omp simd
    for (i=0; i<n; i++)
    {
        gamma[i] = third*asin(3.0*v[i]);
    }
    return;
}
