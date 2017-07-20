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
 * @param[in] n       number of moment tensors
 * @param[in] M1in    first moment tensor packed {m11, m22, m33, m12, m13, m23}
 *                    [6*n]
 * @param[in] M2in    second moment tensor packed {m11, m22, m33, m12, m13, m23}
 *                    [6*n]
 *
 * @param[out] theta  angle between moment tensors M1 and M2 (radians) [n]
 *
 * @result 0 indicates success
 *
 */
int compearth_angleMT(const int n,
                      const double *__restrict__ M1in,
                      const double *__restrict__ M2in,
                      double *__restrict__ theta)
{
    const char *fcnm = "compearth_angleMT\0";
    double M1[9], M2[9], arg, M1_mag, M2_mag, xnum, xden;
    int i, ierr;
    ierr = 0;
    for (i=0; i<n; i++)
    {
        // Convert to 3x3 matrices
        compearth_Mvec2Mmat(1, &M1in[6*i], 1, M1);
        compearth_Mvec2Mmat(1, &M2in[6*i], 1, M2);
        // Compute norms
        ierr = compearth_normMat(1, M1, TWO_NORM, 2.0, &M1_mag);
        ierr = compearth_normMat(1, M2, TWO_NORM, 2.0, &M2_mag);
        xden = M1_mag*M2_mag;
        if (fabs(xden) < 1.e-15) //xden == 0.0)
        {
            printf("%s: Division by zero!\n", fcnm);
            ierr = 1;
            break;
        }
        xnum = cblas_ddot(9, M1, 1, M2, 1);
        arg = xnum/xden;
        if (arg <-1.0 || arg > 1.0)
        {
            printf("%s: Invalid argument to acos %f %f %f!\n",
                   fcnm, arg, xnum, xden);
            ierr = 1;
            break;
        }
        theta[i] = acos(xnum/xden);
    }
    return ierr;
}
