#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define COMPEARTH_PRIVATE_CROSS3 1
#define COMPEARTH_PRIVATE_GEMV3 1
#define COMPEARTH_PRIVATE_GEM3 1
#define COMPEARTH_PRIVATE_GEMT3 1
#include "compearth.h"
#ifdef COMPEARTH_USE_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

/*!
 * @brief C translation of Carl Tape's utility for converting geometrical
 * parameters into moment tensors.
 * https://github.com/carltape/compearth/blob/master/momenttensor/matlab/TT2CMT.m
 *
 * @param[in] gamma      angle (degrees) from DC meridian to MT point [-30,30]
 * @param[in] delta      angle (degrees) from deviatoric plane to MT 
 *                       point [-90,90]
 * @param[in] M0         seismic moment; note rho = sqrt(2)*M0 so set
 *                       M0=1/sqrt(2) for rho=1
 * @param[in] kappa      strike angle (degrees)
 * @param[in] theta      dip angle (degrees)
 * @param[in] sigmaIn    slip (rake) angle (degrees)
 *
 * @param[out] M         moment tensor in Up-South-East coordinates packed
 *                       {rr, tt, pp, rt, rp, tp} [6]
 * @param[out] lam       eigenvalues
 * @param[out] U         3 x 3 bases in USE convention stored in column major
 *                       order [9]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 */
int compearth_tt2cmt(const double gamma,
                     const double delta,
                     const double M0,
                     const double kappa,
                     const double theta,
                     const double sigmaIn, 
                     double M[6], double lam[3], double U[9])
{
    double M9[9], R[9], Yrot[9], V[9], Ux[9], Uxd[9], M6[6],
           north[3], zenith[3], K[3], N[3], S[3], NxS[3], phi, sigma;
    const double neg45 =-45.0;
    int ierr;
    //------------------------------------------------------------------------//
    //
    // basic error handling
    ierr = 1;
    if (gamma <-30.0 || gamma > 30.0)
    {
        fprintf(stderr, "%s: gamma=%f is out of range [-30,30]\n",
                 __func__, gamma);
        goto ERROR;
    }
    if (delta <-90.0 || delta > 90.0)
    {
        fprintf(stderr, "%s: delta=%f is out of range [-90,90]\n",
                __func__, delta);
        goto ERROR;
    }
    // warnings
    ierr = 0;
    sigma = sigmaIn;
    //if (theta == 0.0 && fabs(sigma) < 1.e-14)
    if (fabs(theta) < 1.e-15 && fabs(sigma) < 1.e-14)
    {
        fprintf(stdout,
                "%s: Input fault is horizontal; strike angle %f is undefined\n",
                __func__, kappa);
        fprintf(stdout, "%s: Resetting slip angle to 0\n", __func__);
        sigma = 0.0;
    }
    // moment tensor source type (or pattern)
    compearth_lune2lam(1, &gamma, &delta, &M0, lam);
    // PART 2: moment tensor orientation
    // NOTE: Algorithmically, it would be simpler to compute V directly from the
    // expression in Proposition 2, since this requires fewer calculations.
    // (In the case here, we do not need the fault vectors at all.)
    // The implementaion below is more conceptual.
    // The basis is specified from the components of the north and zenith 
    // vectors.  The output for M and U can be changed by using convert_MT.m
    // or convertv.m.

    // for north-west-up basis (TapeTape2012)
    //north = [1 0 0]'; zenith = [0 0 1]';

    // for south-east-up basis (TapeTape2013)
    north[0] = -1.0; north[1]  = 0.0; north[2]  = 0;
    zenith[0] = 0.0; zenith[1] = 0.0; zenith[2] = 1.0;

    // TT2012, p. 485
    phi = -kappa;
    // TT2012, Eq 27abc
    ierr = ierr + compearth_eulerUtil_rotmat(1, &phi, 3, R);
    gemv3_colMajorNoTrans(R, north, K);
#ifdef DEBUG
    double K3[3]; 
    cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, R, 3,
                north, 1, 0.0, K3, 1);
#endif
    ierr = ierr + compearth_eulerUtil_rotmatGen(1, K, &theta, R);
    gemv3_colMajorNoTrans(R, zenith, N);
#ifdef DEBUG
    double N3[3];
    cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, R, 3,
                zenith, 1, 0.0, N3, 1);
#endif
    ierr = ierr + compearth_eulerUtil_rotmatGen(1, N, &sigma, R);
    gemv3_colMajorNoTrans(R, K, S);
#ifdef DEBUG
    double S3[3];
    cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, R, 3,
                K, 1, 0.0, S3, 1);
    for (int k=0; k<3; k++)
    {
        if (fabs(K3[k] - K[k]) > 1.e-14 ||
            fabs(N3[k] - N[k]) > 1.e-14 ||
            fabs(S3[k] - S[k]) > 1.e-14)
        {
            fprintf(stderr, "%s: Failed gemv\n", __func__);
            ierr = ierr + 1;
        }
    }
#endif
    // TT2012, Eq 28 (or Proposition 2)
    compearth_eulerUtil_rotmat(1, &neg45, 2, Yrot);
    cross3(N, S, NxS); //ierr = cross(3, N, S, NxS);
    V[0] = S[0]; V[3] = NxS[0]; V[6] = N[0];
    V[1] = S[1]; V[4] = NxS[1]; V[7] = N[1];
    V[2] = S[2]; V[5] = NxS[2]; V[8] = N[2];
    gemm3_colMajorNoTransNoTrans(V, Yrot, Ux);
#ifdef DEBUG
    double Uxref[9];
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1.0, V, 3, Yrot, 3, 0.0, Uxref, 3);
#endif
    // Ux*diag(lam)
    Uxd[0] = Ux[0]*lam[0]; Uxd[3] = Ux[3]*lam[1]; Uxd[6] = Ux[6]*lam[2];
    Uxd[1] = Ux[1]*lam[0]; Uxd[4] = Ux[4]*lam[1]; Uxd[7] = Ux[7]*lam[2];
    Uxd[2] = Ux[2]*lam[0]; Uxd[5] = Ux[5]*lam[1]; Uxd[8] = Ux[8]*lam[2];
    // Ux*diag(lam)*ux
    gemm3_colMajorNoTransTrans(Uxd, Ux, M9);
#ifdef DEBUG
    double Mref[9];
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                3, 3, 3, 1.0, Uxd, 3, Ux, 3, 0.0, Mref, 3);
    for (int k=0; k<9; k++)
    {
        if (fabs(Ux[k] - Uxref[k]) > 1.e-14 ||
            fabs(M9[k] - Mref[k])/M0 > 1.e-14)
        {
            fprintf(stderr, "%s: Failed gemm %f %f %f %f\n", __func__,
                    Ux[k], Uxref[k], M9[k], Mref[k]);
            ierr = ierr  + 1;
        }
    }
#endif
    cblas_dcopy(9, Ux, 1, U, 1);
    compearth_Mvec2Mmat(1, M9, 2, M6);
    // convert from north-west-up to up-south-east
    // i1 = 3; i2 = 1;
    // M = convert_MT(i1,i2,M);
    // U = convertv(i1,i2,U);

    // convert moment tensor from south-east-up to up-south-east
    // (note: U is still in south-east-up)
    ierr = compearth_convertMT(1, CE_SEU, CE_USE, M6, M);
    if (ierr != 0)
    {
        fprintf(stderr,
                "%s: Error converting moment tensor coordinate systems\n",
                __func__);
    }
//printf("%e\n%e\n%e\n%e\n%e\n%e\n",M[0],M[1],M[2],M[3],M[4],M[5]);
//printf("U:\n");
//printf("%e %e %e\n", U[0], U[3], U[6]); 
//printf("%e %e %e\n", U[1], U[4], U[7]);
//printf("%e %e %e\n", U[2], U[5], U[8]);
ERROR:;  
    return ierr;
}
