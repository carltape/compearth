#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief C translation of Carl Tape's utility for converting lune coordinates
 *        (gamma, delta, M0) to eigenvalues.
 *
 * @param[in] n        Number of moment tensors.
 * @param[in] gamma    gamma (longitude) angles (degrees) on lune
 *                     such that \f$ \gamma \in [-30,30] \f$.
 *                     This is an array of dimension [n].
 * @param[in] delta    delta (latitude) angles (degrees) on lune 
 *                     such that \f$ \delta \in [-90,90] \f$.
 *                     This is an array of dimension [n].
 * @param[in] M0in     Seismic moment.  This is an array of dimension [n].
 *                     If NULL M0 will be set to 1 for all gamma and delta.
 *
 * @param[out] lam     Set of diagonalized moment tensors in GCMT convention
 *                     note: normalized s.t. each MT has moment M0.  This
 *                     is an array of dimension [3*n].
 *
 * @author Carl Tape and converted to C by Ben Baker.
 *
 * @date 2016 - Ben Baker converted Carl Tape's lune2lam.m to C.
 *
 * @copyright MIT
 *
 */
void compearth_lune2lam(const int n,
                        const double *__restrict__ gamma,
                        const double *__restrict__ delta,
                        const double *__restrict__ M0in,
                        double *__restrict__ lam)
{
    double lamhat[3], beta_rad, gamma_rad;
    int i;
    const double pi180 = M_PI/180.0;
    double beta, rho;
    bool lhaveM0;
    const double sqrt2 = M_SQRT2;
    memset(lam, 0, (size_t) (3*n)*sizeof(double));
    lhaveM0 = true;
    if (M0in == NULL)
    {
        lhaveM0 = false;
    }
    for (i=0; i<n; i++)
    {
        // magnitude of lambda vectors (TT2012 p. 490 text)
        rho = 1.0;
        if (lhaveM0){rho = M0in[i]*sqrt2;}
        beta = 90.0 - delta[i]; // colatitude
        // compute points on lune 
        beta_rad = beta*pi180;
        gamma_rad = gamma[i]*pi180;
        compearth_tape2015Eqn7(beta_rad, gamma_rad, lamhat);
        // apply magnitude
        lam[3*i+0] = rho*lamhat[0];
        lam[3*i+1] = rho*lamhat[1];
        lam[3*i+2] = rho*lamhat[2];
    }
    return;
}
