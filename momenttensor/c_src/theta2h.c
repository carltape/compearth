#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief Computes h from dip angle theta as in Equation 24c of 
 *        Tape and Tape 2015
 *
 * @param[in] n         number of points in array
 * @param[in] theta     dip angle \f$ \theta \in [0, \pi/2 \f$ [n]
 *
 * @param[out] h        h of Equation 24c \f$ h \in [-1,1] \f$ [n]
 *
 * @copyright MIT
 *
 */
void compearth_theta2h(const int n, const double *__restrict__ theta,
                        double *__restrict__ h)
{
    int i;
    #pragma omp simd
    for (i=0; i<n; i++)
    {
        h[i] = cos(theta[i]);
    }
    return;
}
