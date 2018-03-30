#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"
#ifdef COMPEARTH_USE_MKL
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wstrict-prototypes"
#endif
#include <mkl_cblas.h>
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#else
#include <cblas.h>
#endif

/*!
 * @brief Converts v-w coordinates to lune coordinates.
 *
 * @param[in] nv      Number of v coordinates.
 * @param[in] v       v coordinates s.t. \f$ v \in [-1/3, 1/3] \f$.  This
 *                    an array of dimension [nv].
 * @param[in] nw      Number of w coordinates.
 * @param[in] w       w coordinates (similar to lune latitudes) s.t. 
 *                    \f$ w \in [-3\pi/8, 3/pi/8] \f$.  This is an
 *                    array of dimension [nw].
 *
 * @param[out] gamma  Longitude (degrees) \f$ \gamma \in [-30,30] \f$.
 *                    This is an array of dimension [nv].
 * @param[out] delta  lune latitude (degrees) \f$ \delta \in [-90,90] \f$.
 *                    This is an array of dimension [nw].
 *
 * @result 0 indicates success.
 *
 * @author Carl Tape and converted to C by Ben Baker
 *
 * @date 2016
 *
 * @copyright MIT
 *
 */
int compearth_rect2lune(const int nv, const double *__restrict__ v,
                        const int nw, const double *__restrict__ w,
                        double *__restrict__ gamma,
                        double *__restrict__ delta) 
                         
{
    double u[CE_CHUNKSIZE] __attribute__((aligned(64)));
    int i, j, ierr, nuLoc;
/*
    const int maxit = 20;
    const int useHalley = 2;
    const double tol = 1.e-12;
*/
    const double pi38 = 3.0*M_PI/8.0;
    const double pi180i = 180.0/M_PI;
    ierr = 0;
    // convert v -> gamma 
    compearth_v2gamma(nv, v, gamma);
    // convert latitude to colatitude
    for (i=0; i<nw; i=i+CE_CHUNKSIZE)
    {
        nuLoc = MIN(CE_CHUNKSIZE, nw - i);
        for (j=0; j<nuLoc; j++)
        {
            u[j] = pi38 - w[i+j];
        }
        // convert u -> beta
        ierr = compearth_u2beta(nuLoc, u, &delta[i]);
        //ierr = compearth_u2beta(nuLoc, maxit, useHalley, u, tol, &delta[i]);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Computing u2beta\n", __func__);
            return -1;
        }
    }
    // convert to gamma to degrees
    cblas_dscal(nv, pi180i, gamma, 1);
    // convert delta from colatitude to latitude and radians to degrees
    #pragma omp simd
    for (i=0; i<nw; i++)
    {
        delta[i] = 90.0 - pi180i*delta[i];
    }
    return ierr;
}
