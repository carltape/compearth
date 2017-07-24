#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "compearth.h"

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
    double lamWork[3], M[6], Uwork[9];
    int ierr, imt;
    const int isort = 1;
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
    for (imt=0; imt<nmt; imt++)
    {
        // KEY: Convert M into another basis.
        // YOU MUST ALSO CHANGE north AND zneith IN fault2vecang BELOW
        // --> U will be with respect to this basis (from CMTdecom.m)
        // ierr = compearth_convertMT(USE, NWU, &Min[6*imt], M);
        ierr = compearth_convertMT(USE, SEU, &Min[6*imt], M);
        if (ierr != 0)
        {
            printf("%s: Error switching basis\n", fcnm);
            break;
        }
        // PART 1: moment tensor source type (or pattern)
        // Decompose moment tensor into eigenvalues + basis (M = U*lam*U')
        // NOTE: ordering of eigenvalues is important.
        ierr = compearth_CMTdecom(1, M, isort, lamWork, Uwork);
        if (ierr != 0)
        {
            printf("%s: Error decomposing CMT\n", fcnm);
            break;
        }
        // Compute the lune coordinates and magnitude from eigenvalues

    }
    return ierr; 
}
