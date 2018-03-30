#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define COMPEARTH_PRIVATE_GEM3 1
#include "compearth.h"
#ifdef DEBUG
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
#endif

/*!
 * @brief Compute a rotation matrix given an axis and angle
 *
 * @param[in] n     Number of matrices to generate.
 * @param[in] v     Rotation axis.  This is an array of dimension [3 x n]
 *                  with leading dimension 3.
 * @param[in] xi    Rotation angles (degrees).  This is an array of
 *                  dimension [n].
 *
 * @param[in] U     Rotation matrices.  This is an array of dimension
 *                  [3 x 3 x n] where each leading [3 x 3] matrix is
 *                  in column major format.
 *
 * @result 0 indicates success.
 *
 * @date 2016 - Ben Baker translated Carl Tape's rotmat_gen.m to C
 *
 * @copyright MIT
 *
 */
int compearth_eulerUtil_rotmatGen(const int n,
                                  const double *__restrict__ v,
                                  const double *__restrict__ xi,
                                  double *__restrict__ U)
{
    double R1[9], R2[9], R3[9], R4[9], R5[9], R54[9], R543[9], R5432[9];
    double ele, nvph_deg, nvth_deg, rho, vph, vth, vph_deg, vth_deg;
    int i, ierr;
    const double pi180i = 180.0/M_PI;
    ierr = 0;
    for (i=0; i<n; i++)
    {
        compearth_matlab_cart2sph(1, &v[3*i+0], &v[3*i+1], &v[3*i+2],
                                  &vph, &ele, &rho);
        vth = M_PI_2 - ele;
        vth_deg = vth*pi180i;
        vph_deg = vph*pi180i;
        nvph_deg =-vph_deg;
        nvth_deg =-vth_deg;
        ierr = ierr + compearth_eulerUtil_rotmat(1, &nvph_deg, 3, R1);
        ierr = ierr + compearth_eulerUtil_rotmat(1, &nvth_deg, 2, R2);
        ierr = ierr + compearth_eulerUtil_rotmat(1, &xi[i],    3, R3);
        ierr = ierr + compearth_eulerUtil_rotmat(1, &vth_deg,  2, R4);
        ierr = ierr + compearth_eulerUtil_rotmat(1, &vph_deg,  3, R5);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error computing rotation matrices %d\n",
                    __func__, i);
        }
        gemm3_colMajorNoTransNoTrans(R5, R4, R54);
        gemm3_colMajorNoTransNoTrans(R54, R3, R543);
        gemm3_colMajorNoTransNoTrans(R543, R2, R5432);
        gemm3_colMajorNoTransNoTrans(R5432, R1, &U[9*i]);
#ifdef DEBUG
        double U9[9];
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    3, 3, 3, 1.0, R5, 3, R4, 3, 0.0, R54, 3);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    3, 3, 3, 1.0, R54, 3, R3, 3, 0.0, R543, 3);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    3, 3, 3, 1.0, R543, 3, R2, 3, 0.0, R5432, 3);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    3, 3, 3, 1.0, R5432, 3, R1, 3, 0.0, U9, 3);
        for (int k=0; k<9; k++)
        {
            if (fabs(U9[k] - U[9*i+k]) > 1.e-14)
            {
                printf("Computation failed: %f %f\n", U9[k], U[9*i+k]);
                ierr = ierr + 1;
            }
        }
#endif
    }
    return ierr;
}

