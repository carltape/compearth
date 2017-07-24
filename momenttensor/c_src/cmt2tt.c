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
 * @param[in] M             [6 x n] array of moment tensors in CMT
 *                          convention (i.e. UP-SOUTH-EAST).  M is packed:
 *                          \f$
 *                           M = \{M_{rr}, M_{\theta \theta}, M_{\phi \phi}
 *                                 M_{r \theta}, M_{r \phi}, M_{\theta \phi}
 *                          \f$.
 *                          The leading dimension is 6.
 * @param[in] ldisplay      If true then display the results. \n
 *                          Otherwise, this routine will be quiet unless
 *                          an error is encountered.
 * 
 * @param[out] gamma        Angle from DC meridian to MT point.
 *                          Note that $\f \gamma \in [-30, 30] \f$. 
 *                          This is an array of dimension [nmt].
 * @param[out] delta        Angle from deviatoric plane to MT point.
 *                          Note that \f$ \delta \in [-90, 90] \f$.
 *                          This is an array of dimension [nmt].
 * @param[out] M0           Seismic moment in N-m.  This is an array of
 *                          dimension [nmt].
 * @param[in,out] K
 * @param[in,out] N
 * @param[in,out] S
 * @param[in,out] thetadc
 * @param[in,out] lam
 * @param[in,out] U 
 *
 * @result 0 indicates success.
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
                     double *__restrict__ sigma )
{
    const char *fcnm = "compearth_CMT2TT\0";
    double lam[3], M[6], U[9];
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
        ierr = compearth_CMTdecom(1, M, isort, lam, U);
        if (ierr != 0)
        {
            printf("%s: Error decomposing CMT\n", fcnm);
            break;
        }
        // Compute the lune coordinates and magnitude from eigenvalues

    }
    return ierr; 
}
