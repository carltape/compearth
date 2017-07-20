#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief Computes eigenvalues lambda corresponding to the matrix vector
 *        multiply in Equation 7 of Tape and Tape 2015 for a given beta,
 *        gamma coordinate on the lune
 *
 * @param[in] beta      colatitude on lune [0,pi]
 * @param[in] gamma     longitude on lune [-pi/6,pi/6]
 *
 * @param[out] lam      location Lambda(beta,gamma) on lune 
 *
 * @author Ben Baker
 *
 * @copyright MIT
 *
 */
void compearth_tape2015Eqn7(const double beta, const double gamma,
                            double lam[3])
{
    double v1, v2, v3, sbeta;
    const double sqrt6i = 1.0/sqrt(6.0);
    const double sqrt2 = sqrt(2.0);
    const double sqrt3 = sqrt(3.0);
    // rotation matrix s.t. delta = 90 is (1,1,1) and delta =-90 is (-1,-1,-1)
    const double m11 = sqrt3, m12 =-1.0, m13 = sqrt2;
    const double m21 = 0.0,   m22 = 2.0, m23 = sqrt2;
    const double m31 =-sqrt3, m32 =-1.0, m33 = sqrt2;
    // cartesian points 
    sbeta = sin(beta);
    v1 = sbeta*cos(gamma);
    v2 = sbeta*sin(gamma);
    v3 = cos(beta);
    // rotate points
    lam[0] = (m11*v1 + m12*v2 + m13*v3)*sqrt6i;
    lam[1] = (m21*v1 + m22*v2 + m23*v3)*sqrt6i;
    lam[2] = (m31*v1 + m32*v2 + m33*v3)*sqrt6i;
    return;
}
