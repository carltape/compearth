#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"
/*!
 * @brief Computes dip angle from h.
 *
 * @param[in] n       Number of points in arrays.
 * @param[in] h       h in rectilinear space such that
 *                       \f$ h \in [0,1] \f$.
 *                    This is an array of dimension [n].
 *
 * @param[out] theta  Dip angle (radians) such that
 *                      \f$ \theta \in [0, \pi/2] \f$.
 *                    This is an array of dimension [n].
 *
 * @author Ben Baker, ISTI
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
