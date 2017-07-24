#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "compearth.h"
#ifdef COMPEARTH_USE_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#define deg 180.0/M_PI;

/*!
 * @brief Returns the angle between two vectors, in degrees.
 *
 * @param[in] n    Length of vectors.
 * @param[in] va   Initial vector.  This is an array of dimension [n].
 * @param[in] vb   Rotate vector.  This is an array of dimension [n].
 *
 * @result Rotation angle in degrees between va and vb.
 *
 * @author Ben Baker
 *
 * @copyright MIT
 *
 */
double compearth_matlab_fangle(const int n,
                               const double *__restrict__ va,
                               const double *__restrict__ vb)
{
    double angle, magVa, magVb, xden, xnum;
    xnum = cblas_ddot(n, va, 1, vb, 1); 
    magVa = cblas_dnrm2(n, va, 1);
    magVb = cblas_dnrm2(n, vb, 1);
    xden = magVa*magVb;
    angle = NAN;
    if (xden > 0.0){xden = acos(xnum/xden)*deg;}
    return angle;
}
