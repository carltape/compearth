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
 * @author Carl Tape and translated to C by Ben Baker.
 *
 * @copyright MIT
 *
 */
int compearth_angleMT(const int n,
                      const double *__restrict__ M1in,
                      const double *__restrict__ M2in,
                      double *__restrict__ theta)
{
    double M1[9] __attribute__((aligned(64)));
    double M2[9] __attribute__((aligned(64)));
    double arg, M1_mag, M2_mag, xnum, xden;
    int i, ierr;
    ierr = 0;
printf("yes\n");
    for (i=0; i<n; i++)
    {
        // Convert to 3x3 matrices
        compearth_Mvec2Mmat(1, &M1in[6*i], 1, M1);
        compearth_Mvec2Mmat(1, &M2in[6*i], 1, M2);
        // Compute norms
        ierr = compearth_normMat(1, M1, CE_TWO_NORM, 2.0, &M1_mag);
        ierr = compearth_normMat(1, M2, CE_TWO_NORM, 2.0, &M2_mag);
        xden = M1_mag*M2_mag;
        if (fabs(xden) < 1.e-15) //xden == 0.0)
        {
            fprintf(stderr, "%s: Division by zero!\n", __func__);
            ierr = 1;
            break;
        }
        xnum = cblas_ddot(9, M1, 1, M2, 1);
        arg = xnum/xden;
        if (arg <-1.0 || arg > 1.0)
        {
            fprintf(stderr, "%s: Invalid argument to acos %f %f %f!\n",
                    __func__, arg, xnum, xden);
            ierr = 1;
            break;
        }
        theta[i] = acos(xnum/xden);
    }
    return ierr;
}
