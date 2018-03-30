#include <stdio.h>
#include <stdlib.h>
#define COMPEARTH_PRIVATE_GEM3 1
#define COMPEARTH_PRIVATE_GEMT3 1
#include "compearth.h"
#ifdef DEBUG
#ifdef COMPEARTH_USE_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#endif

/*!
 * @brief Transforms a set of 3 x 3 tensors using transformation matrix T.
 *
 * @param[in] nmt    Number of moment tensors.
 * @param[in] T      3 x 3 transformation matrix in column major order.
 * @param[in] Min    [3 x 3 x n] input matrices.  Each [3 x 3] matrix is in
 *                   column major order.
 * @param[out] Mout  [3 x 3 x n] output matrices.  Each [3 x 3] matrix is
 *                   in column major order.
 * @result 0 indicates success.
 *
 * @author Carl Tape and translated to C by Ben Baker.
 *
 */
int compearth_transform_mat(const int nmt,
                            const double *__restrict__ T,
                            const double *__restrict__ Min,
                            double *__restrict__ Mout)
{
    double TM[9];
    int imt; //i, j
    // Compute T*Min*T'
    for (imt=0; imt<nmt; imt++)
    {
#ifdef DEBUG
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    3, 3, 3, 1.0, T, 3, &Min[9*imt], 3, 0.0, TM, 3);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    3, 3, 3, 1.0, TM, 3, T, 3, 0.0, &Mout[9*imt], 3);
#else
        gemm3_colMajorNoTransTrans(T, &Min[9*imt], TM);
        gemm3_colMajorNoTransTrans(TM, T, &Mout[9*imt]);
#endif
    }
    return 0;
}

