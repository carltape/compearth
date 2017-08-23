#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define COMPEARTH_PRIVATE_GEMV3 1
#include "compearth.h"

/*!
 * @brief Converts (full) moment tensors to fault parameters.  Note that
 *        N1 and N2 can point in any direction and you can flip the sign
 *        of both N1 and N2 without changing the system.
 *
 * @param[in] nmt   Number of moment tensors.
 * @param[in] M     6 x nmt set of moment tensors with leading dimension 6.
 *
 * @param[out] nu   Poisson parameter.  This is an array of dimension nmt.
 * @param[out] N1   Unit normal vector for fault plane 1 in same
 *                  basis as M.  This is an array of dimension [3 x nmt] with
 *                  leading dimension 3.
 * @param[out] N2   Unit normal vector for fault plane 2 in same
 *                  basis as M.  This is an array of dimension [3 x nmt] with
 *                  leading dimension 3.
 * @param[out] lam  Eigenvalues of M.  This is an array of dimension [3 x nmt]
 *                  with leading dimension 3.
 *
 * @result 0 indicates success.
 *
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @copyright MIT
 *
 */
int compearth_CMT2faultpar(const int nmt,
                          const double *__restrict__ M,
                          double *__restrict__ nu,
                          double *__restrict__ alpha,
                          double *__restrict__ N1,
                          double *__restrict__ N2,
                          double *__restrict__ lam)
{
    double *U;
    double Ux[9] __attribute__((aligned(64)));
    double Y1[9] __attribute__((aligned(64)));
    double Y2[9] __attribute__((aligned(64)));
    double temp[3] __attribute__((aligned(64)));
    double a2, na2;
    int i, ierr, j;
    const double xvec[3] = {1, 0, 0};
    const int isort = 1;
    // Decomose CMTs for eigenvalues and eigenbasis
    U = (double *) calloc(9*(size_t) nmt, sizeof(double));
    ierr = compearth_CMTdecom(nmt, M, isort, lam, U);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error decomposing CMTs\n", __func__);
        free(U);
        return -1;
    }
    // Get nu and alpha from eigenvalues 
    ierr = compearth_lam2nualpha(nmt, lam, nu, alpha);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Error computing nu and alpha\n", __func__);
        free(U);
        return -1;
    }
    // Eqn 36
    for (i=0; i<nmt; i++)
    {
        for (j=0; j<9; j++){Ux[j] = U[9*i+j];}
        a2  = alpha[i]/2.0;
        na2 =-a2;
        compearth_eulerUtil_rotmat(1, &na2, 2, Y1);
        compearth_eulerUtil_rotmat(1, &a2,  2, Y2);
        // N1 = Ux*Y1*[1, 0, 0]
        gemv3_colMajorNoTrans(Y1, xvec, temp);
        gemv3_colMajorNoTrans(Ux, temp, &N1[3*i]);
        // N2 = Ux*Y2*[1, 0, 0]
        gemv3_colMajorNoTrans(Y2, xvec, temp);
        gemv3_colMajorNoTrans(Ux, temp, &N2[3*i]);
    }
    free(U);
    return 0;
}
