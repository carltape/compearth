#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"
#ifdef COMPEARTH_USE_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

static void gemm3_colMajorNoTransNoTrans(const double *__restrict__ A,
                                         const double *__restrict__ B,
                                         double *__restrict__ C);

/*!
 * @brief Compute a rotation matrix given an axis and angle
 *
 * @param[in] n     number of matrices to generate
 * @param[in] v     rotation axis [3*n]
 * @param[in] xi    rotation angles (degrees) [n]
 *
 * @param[in] U     rotation matrices which each [3 x 3] packed in 
 *                  column major format [9*n]
 *
 * @result 0 indicates success
 *
 * @date 2016 - Ben Baker translated Carl Tape's rotmat_gen.m to C
 *
 * @copyright MIT
 *
 */
int compearth_eulerUtil_rotmatGen(const int n,
                                  const double *__restrict__ v,
                                  const double *xi,
                                  double *__restrict__ U)
{
    const char *fcnm = "compearth_eulerUtil_rotmatGen\0";
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
            printf("%s: Error computing rotation matrices %d\n", fcnm, i);
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

static void gemm3_colMajorNoTransNoTrans(const double *__restrict__ A,
                                         const double *__restrict__ B,
                                         double *__restrict__ C)
{
    // column 1
    C[0] = A[0]*B[0] + A[3]*B[1] + A[6]*B[2];
    C[1] = A[1]*B[0] + A[4]*B[1] + A[7]*B[2];
    C[2] = A[2]*B[0] + A[5]*B[1] + A[8]*B[2];
    // column 2
    C[3] = A[0]*B[3] + A[3]*B[4] + A[6]*B[5];
    C[4] = A[1]*B[3] + A[4]*B[4] + A[7]*B[5];
    C[5] = A[2]*B[3] + A[5]*B[4] + A[8]*B[5];
    // column 3
    C[6] = A[0]*B[6] + A[3]*B[7] + A[6]*B[8];
    C[7] = A[1]*B[6] + A[4]*B[7] + A[7]*B[8];
    C[8] = A[2]*B[6] + A[5]*B[7] + A[8]*B[8];
    return; 
}
