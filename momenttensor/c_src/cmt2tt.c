#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#define COMPEARTH_PRIVATE_CROSS3 1
#define COMPEARTH_PRIVATE_NORM3 1
#define COMPEARTH_PRIVATE_DOT3 1
#include "compearth.h"
#ifdef COMPEARTH_USE_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#define MAXMT 64

static void setZero(const int nmt, const double tol, double *__restrict__ X);
static double wrap360(const double lon);

/*!
 * @brief Converts a moment tensor to six parameters of Tape and Tape 2012.
 *        The reverse program is TT2CMT.c
 *
 * @param[in] nmt           Number of moment tensors.
 * @param[in] M             [6 x nmt] array of moment tensors in CMT
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
int compearth_cmt2tt(const int nmt, const double *__restrict__ Min,
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
           S1[3], S2[3], S3[3], S4[3], N1[3], N2[3], N3[3], N4[3];
    double lamSpace[3*MAXMT];
    double uSpace[9*MAXMT];
    bool lwantLam, lwantU;
    int ierr, ierrAll, imt, j;
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
            S1[j] = S[3*imt+j]; N1[j] = N[3*imt+j];
            S2[j] =-S[3*imt+j]; N2[j] =-N[3*imt+j];
            S3[j] = N[3*imt+j]; N3[j] = S[3*imt+j];
            S4[j] =-N[3*imt+j]; N4[j] =-S[3*imt+j];
        }
    }
    // Clean up
ERROR:;
    Uwork = NULL;
    lamWork = NULL;
    return ierr; 
}

static void faultVec2Ang(const double *__restrict__ S,
                         const double *__restrict__ N,
                         double *theta, double *sigma, double *kappa, double *K)
{
    const char *fcnm = "faultVec2Ang\0";
    double v[3], costh, vnorm;
    const double deg = 180.0/M_PI;
    // South-East-Up (as In TT2012)
    const double zenith[3] = {0, 0, 1};
    const double north[3] = {-1, 0, 0}; 
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

    // Figure 14
    costh = dot3(N, zenith);
    *theta = acos(costh)*deg;

    // Figure 14

    *kappa = wrap360(*kappa);
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

/*
static void cross3(const double *__restrict__ a, 
                   const double *__restrict__ b,
                   double *__restrict__ c)
{
    // Compute cross product
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return;
}

static double norm3(const double *__restrict__ a)
{
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

static double dot3(const double *__restrict__ a, const double *__restrict__ b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
*/

static double wrap360(const double lon)
{
    double lonw;
    bool lpos;
    lpos = false; 
    if (lon > 0.0){lpos = true;}
    lonw = fmod(lon, 360.0); 
    if (lonw == 0.0 && lpos){lonw = 360.0;}
    return lonw;
}
