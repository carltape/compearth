#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define COMPEARTH_PRIVATE_CROSS3 1
#define COMPEARTH_PRIVATE_GEMV3 1
#define COMPEARTH_PRIVATE_GEM3 1
#define COMPEARTH_PRIVATE_GEMT3 1
#include "compearth.h"

/*!
 * @brief C translation of Carl Tape's utility for converting geometrical
 * parameters into moment tensors.
 * https://github.com/carltape/compearth/blob/master/momenttensor/matlab/TT2CMT.m
 *
 * @param[in] nmt        Number of moment tensors.
 * @param[in] gamma      Angle (degrees) from DC meridian to MT point [-30,30].
 *                       This is an array with dimension [nmt].
 * @param[in] delta      Angle (degrees) from deviatoric plane to MT 
 *                       point [-90,90].  This is an array with dimension [nmt].
 * @param[in] M0         Seismic moment.  Note that rho = sqrt(2)*M0 so set
 *                       M0=1/sqrt(2) for rho=1.  This is an array of dimension
 *                       [nmt].
 * @param[in] kappa      Strike angles (degrees) in range [0,360].  This is 
 *                       an array of  dimension [nmt].
 * @param[in] theta      Dip angles (degrees) in range [0,90].  This is an
 *                       array of dimension [nmt].
 * @param[in] sigmaIn    Slip (rake) angles (degrees) in range [-180,180].
 *                       This is an array of dimension [nmt].
 *
 * @param[out] M         Moment tensor in Up-South-East coordinates packed
 *                       {rr, tt, pp, rt, rp, tp}.  This is an array of
 *                       dimension [6 x nmt] with leading dimension 6.
 * @param[out] lam       [3 x nmt] set of eigenvalues corresponding to U with
 *                       leading dimension 3.
 * @param[out] U         3 x 3 bases in SEU convention.  This is an array of
 *                       dimension [3 x 3 x nmt] where each leading 3 x 3
 *                       matrix is in column major order.
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 */
int compearth_TT2CMT(const int nmt,
                     const double *__restrict__ gamma,
                     const double *__restrict__ delta,
                     const double *__restrict__ M0,
                     const double *__restrict__ kappa,
                     const double *__restrict__ theta,
                     const double *__restrict__ sigmaIn, 
                     double *__restrict__ M,
                     double *__restrict__ lam,
                     double *__restrict__ U)
{
    double M9[9*CE_CHUNKSIZE], R[9*CE_CHUNKSIZE],
           M6[6*CE_CHUNKSIZE], K[3*CE_CHUNKSIZE],
           N[3*CE_CHUNKSIZE],
           sigma[CE_CHUNKSIZE], phi[CE_CHUNKSIZE];
    double Yrot[9], V[9], Uxd[9], NxS[3], S[3], *lam3, *Ux;
    const double neg45 =-45.0;
    int i, ierr, ierr1, ierr2, imt, iwarn1, nmtLoc;
    // for north-west-up basis (TapeTape2012)
    //north = [1 0 0]'; zenith = [0 0 1]';
    // for south-east-up basis (TapeTape2013)
    const double north[3]  = {-1.0, 0.0, 0.0};
    const double zenith[3] = { 0.0, 0.0, 1.0};
    //------------------------------------------------------------------------//
    //
    // basic error handling
    ierr = 1;
    if (nmt < 1 || gamma == NULL || delta == NULL || M0 == NULL ||
        kappa == NULL || theta == NULL || sigmaIn == NULL)
    {
        if (nmt < 1)
        {
            fprintf(stderr, "%s: nmt=%d must be positive\n", __func__, nmt);
        }
        if (gamma == NULL){fprintf(stderr, "%s: gamma is NULL\n", __func__);}
        if (kappa == NULL){fprintf(stderr, "%s: kappa is NULL\n", __func__);}
        if (M0 == NULL){fprintf(stderr, "%s: M0 is NULL\n", __func__);}
        if (kappa == NULL){fprintf(stderr, "%s: kappa is NULL\n", __func__);}
        if (theta == NULL){fprintf(stderr, "%s: theta is NULL\n", __func__);}
        if (sigmaIn == NULL)
        {
            fprintf(stderr, "%s: sigmaIn is NULL\n", __func__);
        }
        return -1;
    }
    // Bounds checks and warnings on undefined strike angles
    ierr1 = 0;
    ierr2 = 0;
    iwarn1 = 0;
    for (i=0; i<nmt; i++)
    {
        if (fabs(theta[i]) < 1.e-15 && fabs(sigmaIn[i]) < 1.e-14)
        {
            iwarn1 = iwarn1 + 1;
        }
    }
    if (ierr1 != 0)
    {
        fprintf(stderr, "%s: All gammas must be in range [-30,30]\n", __func__);
        return -1; 
    }
    if (ierr2 != 0)
    {
        fprintf(stderr, "%s: All betas must be in range [-90,90]\n", __func__);
        return -1;
    }
    if (iwarn1 != 0)
    {
        fprintf(stdout,
                "%s: %d input faults are horizontal; strike angles undefined\n",
                __func__, iwarn1);
        fprintf(stdout, "%s: Slip angles will be set to 0\n", __func__);
    }
    // Loop on moment tensor chunks
    for (imt=0; imt<nmt; imt=imt+CE_CHUNKSIZE)
    {
        nmtLoc = MIN(CE_CHUNKSIZE, nmt - imt);
        // Override slip angle
        for (i=0; i<nmtLoc; i++)
        {
            sigma[i] = sigmaIn[imt+i];
            if (fabs(theta[imt+i]) < 1.e-15 && fabs(sigma[imt+i]) < 1.e-14)
            {
                sigma[i] = 0.0;
            } 
        }
        // Moment tensor source type (or pattern)
        compearth_lune2lam(nmtLoc, &gamma[imt], &delta[imt],
                           &M0[imt], &lam[3*imt]);
        // PART 2: moment tensor orientation
        // NOTE: Algorithmically, it would be simpler to compute V directly
        // from the expression in Proposition 2, since this requires fewer
        // calculations.   (In the case here, we do not need the fault 
        // vectors at all.) The implementaion below is more conceptual.
        // The basis is specified from the components of the north and zenith
        // vectors.  The output for M and U can be changed by using convert_MT.m
        // or convertv.m.
        
        // TT2012, p. 485
        for (i=0; i<nmtLoc; i++)
        {
            phi[i] =-kappa[imt+i];
        }
        // TT2012, Eq 27abc
        ierr = compearth_eulerUtil_rotmat(nmtLoc, phi, 3, R);
        if (ierr != 0){goto ERROR;}
        for (i=0; i<nmtLoc; i++)
        {
            gemv3_colMajorNoTrans(&R[9*i], north, &K[3*i]);
        }
        ierr = compearth_eulerUtil_rotmatGen(nmtLoc, K,
                                             &theta[imt], R);
        if (ierr != 0){goto ERROR;}
        for (i=0; i<nmtLoc; i++)
        {
            gemv3_colMajorNoTrans(&R[9*i], zenith, &N[3*i]);
        }
        ierr = compearth_eulerUtil_rotmatGen(nmtLoc, N, sigma, R);
        if (ierr != 0){goto ERROR;}
        // TT2012, Eq 28 (or Proposition 2)
        compearth_eulerUtil_rotmat(1, &neg45, 2, Yrot);
        lam3 = (double *) &lam[3*imt];
        Ux = (double *) &U[9*imt];
        for (i=0; i<nmtLoc; i++)
        {
            gemv3_colMajorNoTrans(&R[9*i], &K[3*i], S); // moved to here
            cross3(&N[3*i], S, NxS); //ierr = cross(3, N, S, NxS);
            V[0] = S[0]; V[3] = NxS[0]; V[6] = N[3*i+0];
            V[1] = S[1]; V[4] = NxS[1]; V[7] = N[3*i+1];
            V[2] = S[2]; V[5] = NxS[2]; V[8] = N[3*i+2];
            gemm3_colMajorNoTransNoTrans(V, Yrot, Ux);
            // Ux*diag(lam)
            Uxd[0] = Ux[0]*lam3[3*i+0];
            Uxd[1] = Ux[1]*lam3[3*i+0];
            Uxd[2] = Ux[2]*lam3[3*i+0];
            Uxd[3] = Ux[3]*lam3[3*i+1];
            Uxd[4] = Ux[4]*lam3[3*i+1];
            Uxd[5] = Ux[5]*lam3[3*i+1];
            Uxd[6] = Ux[6]*lam3[3*i+2];
            Uxd[7] = Ux[7]*lam3[3*i+2];
            Uxd[8] = Ux[8]*lam3[3*i+2];
            // Ux*diag(lam)*ux
            gemm3_colMajorNoTransTrans(Uxd, Ux, &M9[9*i]);
//printf("%e %e %e %e %e %e %e %e %e\n", M9[9*i], M9[9*i+1], M9[9*i+2],
// M9[9*i+3], M9[9*i+4], M9[9*i+5], M9[9*i+6], M9[9*i+7], M9[9*i+8]);
//printf("%e %e %e %e %e %e\n", M9[9*i],M9[9*i+1], M9[9*i+2],M9[9*i+3],M9[9*i+4],M9[9*i+5]);
        }
        lam3 = NULL;
        Ux = NULL;
        // cblas_dcopy(9*nmtLoc, Ux, 1, U[9*i], 1); -> now happens w/ ptrs
        compearth_Mvec2Mmat(nmtLoc, M9, 2, M6);
        // convert from north-west-up to up-south-east
        // i1 = 3; i2 = 1;
        // M = convert_MT(i1,i2,M);
        // U = convertv(i1,i2,U);

        // convert moment tensor from south-east-up to up-south-east
        // (note: U is still in south-east-up)
        ierr = compearth_convertMT(nmtLoc, CE_SEU, CE_USE, M6, &M[6*imt]);
/*
for (i=0; i<nmtLoc; i++)
{
//printf("%e %e %e\n", lam[3*i], lam[3*i+1], lam[9*i+2]);
//printf("%e %e %e %e %e %e %e\n", M[6*(imt+i)+0], M[6*(imt+i)+1], M[6*(imt+i)+2], M[6*(imt+i)+3], M[6*(imt+i)+3], M[6*(imt+i)+4], M[6*(imt+i)+5]);
}
//getchar();
*/
        if (ierr != 0)
        {
            fprintf(stderr,
                "%s: Error converting moment tensor coordinate systems\n",
                __func__);
        }
    } // Loop on moment tensors
//printf("%e\n%e\n%e\n%e\n%e\n%e\n",M[0],M[1],M[2],M[3],M[4],M[5]);
//printf("U:\n");
//printf("%e %e %e\n", U[0], U[3], U[6]); 
//printf("%e %e %e\n", U[1], U[4], U[7]);
//printf("%e %e %e\n", U[2], U[5], U[8]);
ERROR:;  
    return ierr;
}
