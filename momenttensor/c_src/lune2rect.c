#include <stdio.h>
#include <stdlib.h>
#include "compearth.h"

/*!
 * @brief Converts lune coordinates (gamma, delta) to (v,w)
 *
 * @param[in] ng     number of longitudes.
 * @param[in] gamma  longitude angles (degrees) 
 *                   s.t. \f$ \gamma \in [-30,30] \f$ [ng].
 * @param[in] nd     number of latitudes 
 * @param[in] delta  latitude angles (degrees) s.t. 
 *                   s.t. \f$ \delta \in [-90,90] \f$ [nd].
 *
 * @param[out] v     v coordinates s.t. \f$ v \in [-1/3, 1/3] \f$ [ng].
 * @param[out] w     w coordinates s.t. \f$ v \in [-3\pi/8, 3\pi/8] \f$ [nd].
 * 
 * @date 2016 - Ben Baker converted Carl Tape's lune2rect.m to C
 *
 * @copyright MIT
 *
 */
void compearth_lune2rect(const int ng, double *__restrict__ gamma,
                         const int nd, double *__restrict__ delta,
                         double *__restrict__ v,
                         double *__restrict__ w)
{
    double *betaRad, *gammaRad;
    int i;
    const double pi38 = 3.0*M_PI/8.0;
    const double pi180 = M_PI/180.0;
#if __STDC_VERSION__ >= 201112L 
    size_t nbytes;
#ifdef USE_POSIX
    nbytes = (size_t) ng*sizeof(double);
    posix_memalign((void **) &gammaRad, 64, nbytes);
    nbytes = (size_t) nd*sizeof(double);
    posix_memalign((void **) &betaRad, 64, nbytes);
#else
    nbytes = (size_t) ng*sizeof(double);
    gammaRad = (double *) aligned_alloc(64, nbytes);
    nbytes = (size_t) nd*sizeof(double);
    betaRad  = (double *) aligned_alloc(64, nbytes);
#endif
#else
    gammaRad = (double *) calloc((size_t) ng, sizeof(double));
    betaRad  = (double *) calloc((size_t) nd, sizeof(double));
#endif
    for (i=0; i<ng; i++)
    {
        gammaRad[i] = pi180*gamma[i];
    }
    for (i=0; i<nd; i++)
    {
        betaRad[i] = pi180*delta[i];
    }
    compearth_gamma2v(ng, gammaRad, v);
    compearth_beta2u(nd, betaRad, w);
    for (i=0; i<nd; i++)
    {
        w[i] = pi38 - w[i];
    }
    free(gammaRad);
    free(betaRad);
    return;
}
