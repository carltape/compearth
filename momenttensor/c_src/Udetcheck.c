#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define COMPEARTH_PRIVATE_DET3X3 1
#include "compearth.h"
/*
#ifdef COMPEARTH_USE_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
*/

//static double det3x3ColumnMajor(const double *__restrict__ A);

/*!
 * @brief Ensure that U is a rotation matrix with det(U) = 1.
 *
 * @param[in] n      Number of rotation matrices.
 * @param[in] Uin    Rotation matrices to check.  This is an array of dimension
 *                   [3 x 3 x n] where each [3 x 3] matrix is in column major
 *                   format.
 * @param[in] Uout   Rotation matrices with determinants all positive.  This
 *                   is an array of dimension [3 x 3 x n] where each [3 x 3]
 *                   matrix is in column major format.
 * 
 * @result 0 indicates success.
 *
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @copyright MIT
 *
 */
int compearth_Udetcheck(const int n,
                        const double *__restrict__ Uin,
                        double *__restrict__ Uout)
{
    double det;
    int i, imt;
    // Copy Uin to Uout
    memcpy(Uout, Uin, 9*(size_t) n*sizeof(double));
    //cblas_dcopy(9*n, Uin, 1, Uout, 1);
    // Fix
    for (imt=0; imt<n; imt++)
    {
        det = det3x3ColumnMajor(&Uout[9*imt]); 
        // Negate the second column
        if (det < 0.0)
        {
            for (i=0; i<3; i++)
            {
                Uout[9*imt+3+i] =-Uout[9*imt+3+i];
            }
        }
    }
    // Verify
    for (imt=0; imt<n; imt++)
    {
        det = det3x3ColumnMajor(&Uout[9*imt]);
        if (det < 0.0)
        {
            fprintf(stderr, "%s: Error det(u) < 0\n", __func__);
            return -1;
        }
    }
    return 0;
}

/*
static double det3x3ColumnMajor(const double *__restrict__ A)
{
    double det;
    det = A[0]*( A[4]*A[8] - A[5]*A[7]) 
        - A[3]*( A[1]*A[8] - A[2]*A[7])
        + A[6]*( A[1]*A[5] - A[2]*A[4]);
    return det;
}
*/
