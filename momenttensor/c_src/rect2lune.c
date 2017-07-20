#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"
#ifdef COMPEARTH_USE_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

/*!
 * @brief Converts v-w coordinates to lune coordinates
 *
 * @param[in] v       v coordinates s.t. \f$ v \in [-1,3, 1/3] \f$ [nv]
 * @param[in] w       w coordinates (similar to lune latitudes) s.t. 
 *                    \f$ w \in [-3\pi/8, 3/pi/8] \f$ [nw]
 *
 * @param[out] gamma  longitude (degrees) \f$ \gamma \in [-30,30] \f$ [nv]
 * @param[out] delta  lune latitude (degrees) \f$ \delta \in [-90,90] \f$ [nw]
 *
 * @result 0 indicates success
 *
 * @date 2016 - Ben Baker converted Carl Tape's lune2lat.m to C
 *
 * @copyright MIT
 *
 */
int compearth_rect2lune(const int nv, const double *__restrict__ v,
                        const int nw, const double *__restrict__ w,
                        double *__restrict__ gamma,
                        double *__restrict__ delta) 
                         
{
    const char *fcnm = "compearth_rect2lune\0";
    double *u;
    int i, ierr, nu;
    const int maxit = 20;
    const int useHalley = 2;
    const double tol = 1.e-12;
    const double pi38 = 3.0*M_PI/8.0;
    const double pi180i = 180.0/M_PI;
    // convert v -> gamma 
    compearth_v2gamma(nv, v, gamma);
    // convert latitude to colatitude
    nu = nw;
#if __STDC_VERSION__ >= 201112L 
    size_t nbytes;
    nbytes = (size_t) nu*sizeof(double);
#ifdef USE_POSIX
    posix_memalign((void **) &u, 64, nbytes);
#else
    u = (double *) aligned_alloc(64, nbytes);
#endif
#else
    u = (double *) calloc((size_t) nu, sizeof(double));
#endif
    for (i=0; i<nw; i++)
    {
        u[i] = pi38 - w[i];
    }
    // convert u -> beta
    ierr = compearth_u2beta(nu, maxit, useHalley, u, tol, delta);
    if (ierr != 0)
    {
        printf("%s: Computing u2beta\n", fcnm);
    }
    free(u); 
    // convert to gamma to degrees
    cblas_dscal(nv, pi180i, gamma, 1);
    // convert delta from colatitude to latitude and radians to degrees
    #pragma omp simd
    for (i=0; i<nu; i++)
    {
        delta[i] = 90.0 - pi180i*delta[i];
    }
    return ierr;
}
