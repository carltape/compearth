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
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#define MAXMT 64

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
    const char *fcnm = "compearth_CMT2TT\0";
    double *lamdev, *lamiso, *lamWork, *Uwork, Vwork[9], Yrot[9], M[6],
           Sloc[12], Kloc[12], Nloc[12], kappaL[4], sigmaL[4], thetaL[4];
    double lamSpace[3*MAXMT];
    double uSpace[9*MAXMT];
    bool lwantLam, lwantK, lwantN, lwantS, lwantU;
    int itemp[4], ierr, ierrAll, match, imt, j, nMatch;
    const int isort = 1;
    const double rotAngle = 45.0;
    const double tol = 1.e-6;
    ierr = 0;
    if (nmt < 1 || Min == NULL || gamma == NULL || delta == NULL ||
        M0 == NULL || kappa == NULL || theta == NULL || sigma == NULL)
    {
        if (nmt < 1){printf("%s: No moment tensors\n", fcnm);}
        if (Min == NULL){printf("%s: Min is NULL\n", fcnm);}
        if (gamma == NULL){printf("%s: gamma is NULL\n", fcnm);}
        if (delta == NULL){printf("%s: delta is NULL\n", fcnm);}
        if (M0 == NULL){printf("%s: M0 is NULL\n", fcnm);}
        if (kappa == NULL){printf("%s: kappa is NULL\n", fcnm);}
        if (theta == NULL){printf("%s: theta is NULL\n", fcnm);}
        if (sigma == NULL){printf("%s: sigma is NULL\n", fcnm);}
        return -1; 
    }
    // Pair up pointers
    lwantLam = false;
    if (lam != NULL){lwantLam = true;}
    if (lwantLam)
    {
        lamWork = lam;
    }
    else
    {
        if (nmt <= MAXMT)
        {
            lamWork = lamSpace;
        }
        else
        {
            lamWork = (double *) calloc((size_t) (3*nmt), sizeof(double)); 
        } 
    }
    lwantU = false;
    if (U != NULL){lwantU = true;}
    if (lwantU)
    {
        Uwork = U;
    }
    else
    {
        if (nmt <= MAXMT)
        {
            Uwork = uSpace;
        }
        else
        {
            Uwork = (double *) calloc((size_t) (9*nmt), sizeof(double));
        }
    } 
    // Determine the desired fault vectors
    lwantK = false;
    lwantN = false;
    lwantS = false;
    if (K != NULL){lwantK = true;}
    if (N != NULL){lwantN = true;}
    if (S != NULL){lwantS = true;} 
    // This is the numerically expensive part b/c of eigendecomposition
    ierrAll = 0;
    for (imt=0; imt<nmt; imt++)
    {
        // KEY: Convert M into another basis.
        // YOU MUST ALSO CHANGE north AND zenith IN fault2vecang BELOW
        // --> U will be with respect to this basis (from CMTdecom.m)
        // ierr = compearth_convertMT(1, USE, NWU, &Min[6*imt], M);
        ierr = compearth_convertMT(1, USE, SEU, &Min[6*imt], M);
        if (ierr != 0)
        {
            printf("%s: Error switching basis\n", fcnm);
            ierrAll = ierrAll + 1; //break;
        }
        // PART 1: moment tensor source type (or pattern)
        // Decompose moment tensor into eigenvalues + basis (M = U*lam*U')
        // NOTE: ordering of eigenvalues is important.
        ierr = compearth_CMTdecom(1, M, isort, &lamWork[3*imt], &Uwork[9*imt]);
        if (ierr != 0)
        {
            printf("%s: Error decomposing CMT\n", fcnm);
            ierrAll = ierrAll + 1; //break;
        }
    }
    if (ierrAll != 0)
    {
        printf("%s: Error during eigendecomposition\n", fcnm);
        ierr = 1;
        goto ERROR;
    }
    // Compute the lune coordinates and magnitude from eigenvalues
    lamdev = NULL;
    lamiso = NULL;
    ierr = compearth_lam2lune(nmt, lamWork, gamma, delta, M0, thetadc,
                              lamdev, lamiso);
    // Part 2: moment tensor orientation; TT2012, Section 6.3
    ierr = compearth_eulerUtil_rotmat(1, &rotAngle, 2, Yrot);
    if (ierr != 0)
    {
        printf("%s: Error computing rotation matrix\n", fcnm);
        return -1;
    }
    // Compute candidate fault vectors
    for (imt=0; imt<nmt; imt++)
    {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    3, 3, 3, 1.0, &Uwork[9*imt], 3, Yrot, 3,
                    0.0, Vwork, 3); // V = U*Yrot (TT2012, p. 487)
        S[3*imt+0] = Vwork[0];
        S[3*imt+1] = Vwork[1];
        S[3*imt+2] = Vwork[2];
        N[3*imt+0] = Vwork[6];
        N[3*imt+1] = Vwork[7];
        N[3*imt+2] = Vwork[8];
    }
    // Reassign ~0 elements to 0; ~1 elements to 1, and ~-1 elements to -1.
    setZero(nmt, tol, S);
    setZero(nmt, tol, N);
    // Compute fault angles for four possible combinations (TT2012 Figure 15)
    for (imt=0; imt<nmt; imt++)
    {
        for (j=0; j<3; j++)
        {
            Sloc[3*0+j] = S[3*imt+j]; Nloc[3*0+j] = N[3*imt+j];
            Sloc[3*1+j] =-S[3*imt+j]; Nloc[3*1+j] =-N[3*imt+j];
            Sloc[3*2+j] = N[3*imt+j]; Nloc[3*2+j] = S[3*imt+j];
            Sloc[3*3+j] =-N[3*imt+j]; Nloc[3*3+j] =-S[3*imt+j];
        }
        faultVec2Ang(&Sloc[0], &Nloc[0], &thetaL[0], &sigmaL[0],
                     &kappaL[0], &Kloc[0], &ierr);
        faultVec2Ang(&Sloc[3], &Nloc[3], &thetaL[1], &sigmaL[1],
                     &kappaL[1], &Kloc[3], &ierr);
        faultVec2Ang(&Sloc[6], &Nloc[6], &thetaL[2], &sigmaL[2],
                     &kappaL[2], &Kloc[6], &ierr);
        faultVec2Ang(&Sloc[9], &Nloc[9], &thetaL[3], &sigmaL[3],
                     &kappaL[3], &Kloc[9], &ierr);
        // There are four combinations of N and S that represent a double couple
        // moment tensor, as shown in Figure 15 of TT2012.
        // From these four combinations, there are two possible fault planes.
        // We want to isoate the combination that is within the bounding
        // region shown in Figures 16 and B1.
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
                printf("%s: Failed to pick a fault plane\n", fcnm);
                ierr = 1;
                goto ERROR;
            }
        }
        else if (nMatch == 3)
        {
            printf("%s: Warning mt on bdry of orientation domain 3 candiates\n",
                   fcnm);
            printf("%s: thetas: %e %e %e\n", fcnm, thetaL[0], thetaL[1], thetaL[2]);
            printf("%s: sigmas: %e %e %e\n", fcnm, sigmaL[0], sigmaL[1], sigmaL[2]);
            printf("%s: kappas: %e %e %e\n", fcnm, kappaL[0], kappaL[1], kappaL[2]);
            // Just take the first one
            match = itemp[0]; 
        }
        else if (nMatch == 4)
        {
            printf("%s: Error not yet programmed\n", fcnm);
            ierr = 1;
            goto ERROR;
        }
        else
        {
            printf("%s: Error no match\n", fcnm);
            ierr = 1;
            goto ERROR;
        }
        // Select the angle
        kappa[imt] = kappaL[match];
        sigma[imt] = sigmaL[match];
        theta[imt] = thetaL[match];
        // Fault vectors
        if (lwantK)
        {
            for (j=0; j<3; j++){K[3*imt+j] = Kloc[3*match+j];}
        }
        if (lwantN)
        {
            for (j=0; j<3; j++){N[3*imt+j] = Nloc[3*match+j];}
        }
        if (lwantS)
        {
            for (j=0; j<3; j++){S[3*imt+j] = Sloc[3*match+j];}
        }
 
    }
    // Clean up
ERROR:;
    Uwork = NULL;
    lamWork = NULL;
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
    const char *fcnm = "faultVec2Ang\0";
    double v[3], costh, vnorm;
    int ierr1;
    const double deg = 180.0/M_PI;
    // South-East-Up (as In TT2012)
    const double zenith[3] = {0, 0, 1};
    const double negZenith[3] = {0, 0, -1};
    const double north[3] = {-1, 0, 0}; 
    *ierr = 0;
    *kappa = NAN;
    *theta = NAN;
    *sigma = NAN;
    // Strike vector from TT2012, Eqn 29
    cross3(zenith, N, v);
    vnorm = norm3(v);
    if (vnorm < DBL_EPSILON) //== 0.0)
    {
        printf("%s: Horizontal fault -- strike vector is same as slip vector\n",
               fcnm);
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
   double absX[3], dmax;
   int i, idmax;
   // Compute the largest element of abs(X)
   idmax = cblas_idamax(3*nmt, X, 1); 
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
    const char *fcnm = "pickP1";
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
    printf("%s: Error no selection criterion was met\n", fcnm);
    return ipick;
}

