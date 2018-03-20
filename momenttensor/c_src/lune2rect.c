#include <stdio.h>
#include <stdlib.h>
#include "compearth.h"

/*!
 * @brief Converts lune coordinates (gamma, delta) to (v, w) coordinates.
 *
 * @param[in] ng     Number of longitudes.
 * @param[in] gamma  Longitude angles (degrees) such that
 *                   \f$ \gamma \in [-30,30] \f$.  This is an array
 *                   of dimension [ng].
 * @param[in] nd     Number of latitudes.
 * @param[in] delta  latitude angles (degrees) such that
 *                   \f$ \delta \in [-90,90] \f$.  This is an array
 *                   of dimension [nd].
 *
 * @param[out] v     v coordinates such that \f$ v \in [-1/3, 1/3] \f$.
 *                   This is an array of dimension [ng].
 * @param[out] w     w coordinates s.t. \f$ w \in [-3\pi/8, 3\pi/8] \f$.
 *                   This is an array of dimension [nd].
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
    double betaRad[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double gammaRad[CE_CHUNKSIZE] __attribute__((aligned(64)));
    int i, j, ndLoc, ngLoc;
    const double pi38 = 3.0*M_PI/8.0;
    const double pi180 = M_PI/180.0;
    // Compute v from gamma
    for (i=0; i<ng; i=i+CE_CHUNKSIZE)
    {
        ngLoc = MIN(CE_CHUNKSIZE, ng - i); 
        for (j=0; j<ngLoc; j++)
        {
            gammaRad[j] = pi180*gamma[i+j];
        }
        compearth_gamma2v(ngLoc, gammaRad, &v[i]);
    }
    // Compute w from beta
    for (i=0; i<nd; i=i+CE_CHUNKSIZE)
    {
        ndLoc = MIN(CE_CHUNKSIZE, nd - i); 
        for (j=0; j<ndLoc; j++)
        {
            betaRad[j] = pi180*(90.0 - delta[i+j]);
        }
        compearth_beta2u(ndLoc, betaRad, &w[i]);
        for (j=0; j<ndLoc; j++)
        {
            w[i+j] = pi38 - w[i+j];
        }
    }
    return;
}
