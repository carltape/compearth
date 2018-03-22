#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#define COMPEARTH_PRIVATE_CROSS3 1
#define COMPEARTH_PRIVATE_NORM3 1
#define COMPEARTH_PRIVATE_DOT3 1
#define COMPEARTH_PRIVATE_WRAP360 1
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

//#define MAXMT 64

static void setZero(const int nmt, const double tol, double *__restrict__ X);
static void faultVec2Ang(const double *__restrict__ S,
                         const double *__restrict__ N,
                         double *theta, double *sigma,
                         double *kappa, double *__restrict__ K, int *ierr);
static int pickP1(const double thetaA, const double sigmaA, const double kappaA,
                  const double thetaB, const double sigmaB, const double kappaB,
                  const double tol, const int p1, const int p2);

/*!
 * @brief Converts a moment tensor to six parameters of Tape and Tape 2012.
 *        The reverse program is TT2CMT.c
 *
 * @param[in] nmt           Number of moment tensors.
 * @param[in] Min           [6 x nmt] array of moment tensors in CMT
 *                          convention (i.e. UP-SOUTH-EAST).  M is packed:
 *                          \f$
 *                           M = \{ M_{rr}, M_{\theta \theta}, M_{\phi \phi}
 *                                  M_{r \theta}, M_{r \phi}, M_{\theta \phi} \}
 *                          \f$.
 *                          The leading dimension is 6.
 * @param[in] ldisplay      If true then display the results. \n
 *                          Otherwise, this routine will be quiet unless
 *                          an error is encountered.
 * 
 * @param[out] gamma        Angle from DC meridian to MT point.
 *                          Note that \f$ \gamma \in [-30, 30] \f$. 
 *                          This is an array of dimension [nmt].
 * @param[out] delta        Angle from deviatoric plane to MT point.
 *                          Note that \f$ \delta \in [-90, 90] \f$.
 *                          This is an array of dimension [nmt].
 * @param[out] M0           Seismic moment in N-m.  This is an array of
 *                          dimension [nmt].
 * @param[out] kappa        Strike angle \f$ \kappa \in [0,360] \f$.
 *                          This is an array of dimension [nmt].
 * @param[out] theta        Dip angle \f$ \theta \in [0,90] \f$.
 *                          This is an array of dimension [nmt].
 * @param[out] sigma        Slip (or rake) angle \f$ \sigma \in [-90,90] \f$.
 *                          This is an array of dimension [nmt].
 * @param[out] K            If K is not NULL then it is the strike vector
 *                          (SOUTH-EAST-UP).  In this case K is an array
 *                          of dimension [3 x nmt] with leading dimension 3.
 * @param[out] N            If N is not NULL then it is the normal vector
 *                          (SOUTH-EAST-UP).  In this case N is an array
 *                          of dimension [3 x nmt] with leading dimension 3.
 * @param[out] S            If S Is not NULL then it is the slip vector
 *                          (SOUTH-EAST-UP).  In this case S is an array
 *                          of dimension [3 x nmt] with leading dimension 3.
 * @param[out] thetadc      If thetadc is not NULL then it is angle between 
 *                          the moment tensor and the double couple where
 *                          \f$ \theta_{DC} \in [0,90] \f$.
 *                          In this case thetadc is an array of dimension [nmt].
 * @param[out] lam          If lam is not NULL then these are the eigenvalues
 *                          corresponding to the moment tensor.  In this
 *                          case lam is an array of dimension [3 x nmt] 
 *                          with leading dimension 3.
 * @param[out] U            If U is not NULL then this is the basis in 
 *                          SOUTH-EAST-UP for the moment tensor.  In this
 *                          case U an array of dimension [3 x 3 x nmt] with
 *                          leading dimension 9 and the understanding that
 *                          each 3 x 3 matrix is in column major format.
 *
 * @result 0 indicates success.
 *
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @date July 2017
 * 
 * @copyright MIT
 * 
 */
int compearth_CMT2TT(const int nmt, const double *__restrict__ Min,
                     const bool ldisplay,
                     double *__restrict__ gamma,
                     double *__restrict__ delta,
                     double *__restrict__ M0,
                     double *__restrict__ kappa,
                     double *__restrict__ theta,
                     double *__restrict__ sigma,
                     double *__restrict__ K, double *__restrict__ N,
                     double *__restrict__ S, double *__restrict__ thetadc,
                     double *__restrict__ lam, double *__restrict__ U)
{
    double *lamdev, *lamiso, Vwork[9], Yrot[9], M[6],
           Sloc[12], Kloc[12], Nloc[12], kappaL[4], sigmaL[4], thetaL[4];
    double *thetadcWork;
    double Uwork[9*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double lamWork[3*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double Nwork[3*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double Swork[3*CE_CHUNKSIZE] __attribute__((aligned(64)));
    //double lamSpace[3*MAXMT];
    //double uSpace[9*MAXMT];
    bool lwantLam, lwantK, lwantN, lwantS, lwantU;
    int itemp[4], ierr, ierrAll, match, imt, j, jmt, nmtLoc, nMatch;
    const int isort = 1;
    const double rotAngle = 45.0;
    const double tol = 1.e-6;
    ierr = 0;
    if (nmt < 1 || Min == NULL || gamma == NULL || delta == NULL ||
        M0 == NULL || kappa == NULL || theta == NULL || sigma == NULL)
    {
        if (nmt < 1){fprintf(stderr, "%s: No moment tensors\n", __func__);}
        if (Min == NULL){fprintf(stderr, "%s: Min is NULL\n", __func__);}
        if (gamma == NULL){fprintf(stderr, "%s: gamma is NULL\n", __func__);}
        if (delta == NULL){fprintf(stderr, "%s: delta is NULL\n", __func__);}
        if (M0 == NULL){fprintf(stderr, "%s: M0 is NULL\n", __func__);}
        if (kappa == NULL){fprintf(stderr, "%s: kappa is NULL\n", __func__);}
        if (theta == NULL){fprintf(stderr, "%s: theta is NULL\n", __func__);}
        if (sigma == NULL){fprintf(stderr, "%s: sigma is NULL\n", __func__);}
        return -1; 
    }
    // Determine if eigenvalue/eigenvector decomposition is requested
    lwantLam = false;
    if (lam != NULL){lwantLam = true;}
    lwantU = false;
    if (U != NULL){lwantU = true;}
    // Determine the desired fault vectors
    lwantK = false;
    lwantN = false;
    lwantS = false;
    if (K != NULL){lwantK = true;}
    if (N != NULL){lwantN = true;}
    if (S != NULL){lwantS = true;} 
    // Loop on moment tensor chunks
    for (jmt=0; jmt<nmt; jmt=jmt+CE_CHUNKSIZE)
    {
        nmtLoc = MIN(CE_CHUNKSIZE, nmt - jmt);
        // This is the numerically expensive part b/c of eigendecomposition
        ierrAll = 0;
        for (imt=0; imt<nmtLoc; imt++)
        {
            // KEY: Convert M into another basis.
            // YOU MUST ALSO CHANGE north AND zenith IN fault2vecang BELOW
            // --> U will be with respect to this basis (from CMTdecom.m)
            // ierr = compearth_convertMT(1, CE_USE, CE_NWU, &Min[6*imt], M);
            ierr = compearth_convertMT(1, CE_USE, CE_SEU,
                                       &Min[6*(jmt+imt)], M);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Error switching basis\n", __func__);
                ierrAll = ierrAll + 1; //break;
            }
            // PART 1: moment tensor source type (or pattern)
            // Decompose moment tensor into eigenvalues + basis (M = U*lam*U')
            // NOTE: ordering of eigenvalues is important.
            ierr = compearth_CMTdecom(1, M, isort, &lamWork[3*imt],
                                      &Uwork[9*imt]);
            if (ierr != 0)
            {
                fprintf(stderr, "%s: Error decomposing CMT\n", __func__);
                ierrAll = ierrAll + 1; //break;
            }
        }
        if (ierrAll != 0)
        {
            fprintf(stderr, "%s: Error during eigendecomposition\n", __func__);
            ierr = 1;
            goto ERROR;
        }
        // Compute the lune coordinates and magnitude from eigenvalues
        lamdev = NULL;
        lamiso = NULL;
        thetadcWork = NULL;
        if (thetadc != NULL){thetadcWork = &thetadc[jmt];}
        ierr = compearth_lam2lune(nmtLoc, lamWork, &gamma[jmt], &delta[jmt],
                                  &M0[jmt], thetadcWork,
                                  lamdev, lamiso);
        // Part 2: moment tensor orientation; TT2012, Section 6.3
        ierr = compearth_eulerUtil_rotmat(1, &rotAngle, 2, Yrot);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error computing rotation matrix\n", __func__);
            return -1;
        }
        // Compute candidate fault vectors
        for (imt=0; imt<nmtLoc; imt++)
        {
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        3, 3, 3, 1.0, &Uwork[9*imt], 3, Yrot, 3,
                        0.0, Vwork, 3); // V = U*Yrot (TT2012, p. 487)
            Swork[3*imt+0] = Vwork[0];
            Swork[3*imt+1] = Vwork[1];
            Swork[3*imt+2] = Vwork[2];
            Nwork[3*imt+0] = Vwork[6];
            Nwork[3*imt+1] = Vwork[7];
            Nwork[3*imt+2] = Vwork[8];
        }
        // Save lambda and U
        if (lwantLam){cblas_dcopy(3*nmtLoc, lamWork, 1, &lam[3*jmt], 1);}
        if (lwantU){cblas_dcopy(9*nmtLoc, Uwork, 1, &U[9*jmt], 1);}
        // Reassign ~0 elements to 0; ~1 elements to 1, and ~-1 elements to -1.
        setZero(nmtLoc, tol, Swork);
        setZero(nmtLoc, tol, Nwork);
        // Compute fault angles for four possible combinations (TT2012 Fig 15)
        for (imt=0; imt<nmtLoc; imt++)
        {
            for (j=0; j<3; j++)
            {
                Sloc[3*0+j] = Swork[3*imt+j]; Nloc[3*0+j] = Nwork[3*imt+j];
                Sloc[3*1+j] =-Swork[3*imt+j]; Nloc[3*1+j] =-Nwork[3*imt+j];
                Sloc[3*2+j] = Nwork[3*imt+j]; Nloc[3*2+j] = Swork[3*imt+j];
                Sloc[3*3+j] =-Nwork[3*imt+j]; Nloc[3*3+j] =-Swork[3*imt+j];
            }
            faultVec2Ang(&Sloc[0], &Nloc[0], &thetaL[0], &sigmaL[0],
                         &kappaL[0], &Kloc[0], &ierr);
            faultVec2Ang(&Sloc[3], &Nloc[3], &thetaL[1], &sigmaL[1],
                         &kappaL[1], &Kloc[3], &ierr);
            faultVec2Ang(&Sloc[6], &Nloc[6], &thetaL[2], &sigmaL[2],
                         &kappaL[2], &Kloc[6], &ierr);
            faultVec2Ang(&Sloc[9], &Nloc[9], &thetaL[3], &sigmaL[3],
                         &kappaL[3], &Kloc[9], &ierr);
            // There are four combinations of N and S that represent a double
            // copule moment tensor, as shown in Figure 15 of TT2012.
            // From these four combinations, there are two possible fault
            // planes.  We want to isoate the combination that is within the
            // boundingregion shown in Figures 16 and B1.
            memset(itemp, 0, 4*sizeof(int));
            nMatch = 0;
            for (j=0; j<4; j++)
            {
                if (thetaL[j] <= 90.0 + tol && fabs(sigmaL[j]) <= 90.0 + tol)
                { 
                    itemp[nMatch] = j; //bmatch[j] = true;
                    nMatch = nMatch + 1;
                }
            }
            if (nMatch == 1)
            {
                match = itemp[0];
            }
            else if (nMatch == 2)
            {
                match = pickP1(thetaL[itemp[0]], sigmaL[itemp[0]], kappaL[itemp[0]],
                               thetaL[itemp[1]], sigmaL[itemp[1]], kappaL[itemp[1]],
                               tol, itemp[0], itemp[1]);
                if (match < 0)
                {
                    fprintf(stderr, "%s: Failed to pick a fault plane\n",
                             __func__);
                    ierr = 1;
                    goto ERROR;
                }
            }
            else if (nMatch == 3)
            {
                fprintf(stdout,
                  "%s: Warning mt on bdry of orientation domain 3 candidates\n",
                  __func__);
                fprintf(stdout, "%s: thetas: %e %e %e\n", __func__,
                        thetaL[0], thetaL[1], thetaL[2]);
                fprintf(stdout, "%s: sigmas: %e %e %e\n", __func__,
                        sigmaL[0], sigmaL[1], sigmaL[2]);
                fprintf(stdout, "%s: kappas: %e %e %e\n", __func__,
                        kappaL[0], kappaL[1], kappaL[2]);
                // Just take the first one
                match = itemp[0]; 
            }
            else if (nMatch == 4)
            {
                fprintf(stderr, "%s: Error not yet programmed\n", __func__);
                ierr = 1;
                goto ERROR;
            }
            else
            {
                fprintf(stderr, "%s: Error no match\n", __func__);
                ierr = 1;
                goto ERROR;
            }
            // Select the angle
            kappa[jmt+imt] = kappaL[match];
            sigma[jmt+imt] = sigmaL[match];
            theta[jmt+imt] = thetaL[match];
            // Fault vectors
            if (lwantK)
            {
                for (j=0; j<3; j++){K[3*(jmt+imt)+j] = Kloc[3*match+j];}
            }
            if (lwantN)
            {
                for (j=0; j<3; j++){N[3*(jmt+imt)+j] = Nloc[3*match+j];}
            }
            if (lwantS)
            {
                for (j=0; j<3; j++){S[3*(jmt+imt)+j] = Sloc[3*match+j];}
            }
        }
    } // Loop on moment tensor chunks
    if (ldisplay)
    {
        fprintf(stderr, "%s: ldisply not yet supported\n", __func__);
    }
ERROR:;
/*
    Uwork = NULL;
    lamWork = NULL;
*/
    return ierr; 
}
/*!
 * @brief Returns fault angles in degrees.  Assumes input vectors are in 
 *        the South-East-Up basis.
 *
 * @author Carl Tape and converted to C by Ben Baker
 *
 * @copyright MIT
 *
 */
static void faultVec2Ang(const double *__restrict__ S,
                         const double *__restrict__ N,
                         double *theta, double *sigma,
                         double *kappa, double *__restrict__ K, int *ierr)
{
    double v[3], costh, vnorm;
    int ierr1;
    const double deg = 180.0/M_PI;
    // South-East-Up (as In TT2012)
    const double zenith[3] = {0, 0, 1};
    const double negZenith[3] = {0, 0, -1};
    const double north[3] = {-1, 0, 0}; 
    *ierr = 0;
    *kappa = (double) NAN;
    *theta = (double) NAN;
    *sigma = (double) NAN;
    // Strike vector from TT2012, Eqn 29
    cross3(zenith, N, v);
    vnorm = norm3(v);
    if (vnorm < DBL_EPSILON) //== 0.0)
    {
        fprintf(stderr,
              "%s: Horizontal fault -- strike vector is same as slip vector\n",
                __func__);
        K[0] = S[0];
        K[1] = S[1];
        K[2] = S[2];
    }
    else
    {
        K[0] = v[0]/vnorm;
        K[1] = v[1]/vnorm;
        K[2] = v[2]/vnorm;
    }
    // Figure 14
    *kappa = compearth_eulerUtil_fangleSigned(3, north, K, negZenith, &ierr1);
    if (ierr1 != 0){*ierr = *ierr + 1;}
    *kappa = wrap360(*kappa);
    // Figure 14
    costh = dot3(N, zenith);
    *theta = acos(costh)*deg;
    // Figure 14
    *sigma = compearth_eulerUtil_fangleSigned(3, K, S, N, &ierr1);
    if (ierr1 != 0){*ierr = *ierr + 1;}
    return;
}

static void setZero(const int nmt, const double tol, double *__restrict__ X)
{
   double dmax;
   int i, idmax;
   // Compute the largest element of abs(X)
   idmax = (int) cblas_idamax(3*nmt, X, 1); 
   dmax = X[idmax];
   // Elements near zero whilst trying to eliminate round-off errors
   #pragma omp simd
   for (i=0; i<3*nmt; i++)
   {
       if (fabs(X[i]/dmax) < tol){X[i] = 0.0;}
   }
   #pragma omp simd 
   for (i=0; i<3*nmt; i++)
   { 
       if (fabs(X[i] - 1.0) < tol){X[i] =-1.0;}
   }
   #pragma omp simd
   for (i=0; i<3*nmt; i++)
   {
       if (fabs(X[i] + 1.0) < tol){X[i] = 1.0;}
   }
   return;
}

static int pickP1(const double thetaA, const double sigmaA, const double kappaA,
                  const double thetaB, const double sigmaB, const double kappaB,
                  const double tol, const int p1, const int p2)
{
    int ipick;
    ipick =-1;
    if (fabs(thetaA - 90.0) < tol)
    {
        if (kappaA < 180.0){ipick = p1;}
        if (kappaB < 180.0){ipick = p2;}
        return ipick;
    }
    if (fabs(sigmaA - 90.0) < tol)
    {
        if (kappaA < 180.0){ipick = p1;}
        if (kappaB < 180.0){ipick = p2;}
        return ipick;
    }
    if (fabs(sigmaA + 90.0) < tol)
    {
        if (kappaA < 180.0){ipick = p1;}
        if (kappaB < 180.0){ipick = p2;}
        return ipick;
    }
    fprintf(stderr, "%s: Error no selection criterion was met\n", __func__);
    fprintf(stderr, "thetaA,sigmaA,kappaA=%f,%f,%f\n", thetaA, sigmaA, kappaA);
    fprintf(stderr, "thetaB,sigmaB,kappaB=%f,%f,%f\n", thetaB, sigmaB, kappaB);
    return ipick;
}

