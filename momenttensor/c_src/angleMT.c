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
 * @brief Computes the angle between moment tensors.  This is from 
 *        Equation 15 of Tape and Tape 2016 - A confidence parameter for 
 *        seismic moment tensors
 *
 * @param[in] n       Number of moment tensors.
 * @param[in] M1in    First moment tensor packed {m11, m22, m33, m12, m13, m23}.
 *                    This is an array of dimension [6 x n] with leading
 *                    dimension 6.
 * @param[in] M2in    Second moment tensor packed {m11, m22, m33, m12, m13, m23}.
 *                    This is an array of dimension [6 x n] with leading
 *                    dimension 6.
 *
 * @param[out] theta  Angle between moment tensors M1 and M2 in radians. 
 *                    This is an array of dimension [n].
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 * @copyright MIT
 *
 */
int compearth_angleMT(const int n,
                      const double *__restrict__ M1in,
                      const double *__restrict__ M2in,
                      double *__restrict__ theta)
{
    double M1[9*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double M2[9*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double M1_mag[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double M2_mag[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double xden[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double arg[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double xnum[CE_CHUNKSIZE] __attribute__((aligned(64)));
    int i, ierr, imt, nmtLoc;
    ierr = 0;
    for (imt=0; imt<n; imt=imt+CE_CHUNKSIZE)
    {
        nmtLoc = MIN(CE_CHUNKSIZE, n - imt);
        // Convert to 3x3 matrices
        compearth_Mvec2Mmat(nmtLoc, &M1in[6*imt], 1, M1);
        compearth_Mvec2Mmat(nmtLoc, &M2in[6*imt], 1, M2);
        // Compute norms
        ierr = compearth_normMat(nmtLoc, M1, CE_TWO_NORM, 2.0, M1_mag);
        ierr = compearth_normMat(nmtLoc, M2, CE_TWO_NORM, 2.0, M2_mag);
        // Compute the denomimators
        for (i=0; i<nmtLoc; i++)
        {
            xden[i] = M1_mag[i]*M2_mag[i];
        }
        // Compute the numerators
        for (i=0; i<nmtLoc; i++)
        {
            xnum[i] = cblas_ddot(9, &M1[9*i], 1, &M2[9*i], 1);
        }
        // Compute argument and do all checks
        for (i=0; i<nmtLoc; i++)
        {
            // Can't divide by zero - force this to work b/c M1 or M2 is
            // equivalently 0.
            if (xden[i] == 0.0)
            {
                xnum[i] = 0.0;
                xden[i] = 1.0;
                fprintf(stderr, "%s: Denominator is 0\n", __func__);
            }
            arg[i] = xnum[i]/xden[i];
            // Correct for numerical errors
            if (arg[i] <-1.0)
            {
                fprintf(stdout, "%s: Warning arg is %f setting to -1\n",
                         __func__, arg[i]); 
                arg[i] =-1.0;
            }
            if (arg[i] > 1.0)
            {
                fprintf(stdout, "%s: Warning arg is %f setting to +1\n",
                        __func__, arg[i]);
                arg[i] = 1.0;
            }
        }
        // Finally compute the angle
        for (i=0; i<nmtLoc; i++)
        {
            theta[imt+i] = acos(arg[i]);
        }
    }
    return ierr;
}
