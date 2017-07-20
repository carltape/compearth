#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief Computes u from lune colatitude beta as in Equation 24a of
 *        Tape and Tape 2015
 *
 * @param[in] n         number of points in arrays
 * @param[in] beta      lune colatitude \f$ \beta \in [0, \pi] \f$ [n]
 *
 * @param[out] u        u of Equation 24a \f$ u \in [0, 3\pi/4 ] \f$ [n]
 *
 * @date 2016 - Ben Baker converted Carl Tape's beta2u.m to C
 *
 * @copyright MIT
 *
 */
void compearth_beta2u(const int n, const double *__restrict__ beta,
                      double *__restrict__ u)
{
    double sin2beta, sin4beta;
    int i;
    #pragma omp simd
    for (i=0; i<n; i++) 
    {
        sin2beta = sin(2.0*beta[i]);
        sin4beta = sin(4.0*beta[i]);
        u[i] = 0.75*beta[i] - 0.5*sin2beta + 0.0625*sin4beta;
    }
    return;
}
