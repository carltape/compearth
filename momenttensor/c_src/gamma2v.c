#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief Computes v from lune longitude gamma as in Equation 24b of 
 *        Tape and Tape 2015.
 *
 * @param[in] n         number of points in arrays
 * @param[in] gamma     lune longitude (radians) [n].
 *                      note that \f$ \gamma \in [-\pi/6, \pi/6] \f$. 
 *
 * @param[out] v        v of equation 24b [n].  note that
 *                      \f$ v \in [-1/3,1/3] \f$.
 *
 * @date 2016 - Ben Baker converted Carl Tape's gamma2v.m to C
 *
 * @copyright MIT
 *
 */
void compearth_gamma2v(const int n, const double *__restrict__ gamma,
                       double *__restrict__ v)
{
    int i;
    const double third = 1.0/3.0;
    #pragma omp simd
    for (i=0; i<n; i++) 
    {    
        v[i] = third*sin(3.0*gamma[i]);
    }    
    return;
}
