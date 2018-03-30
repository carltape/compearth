#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

//#define MAXMT 64

/*!
 * @brief Convert eigenvalues to lune coordiantes (gamma, delta, M0)
 *
 * @param[in] nmt       Number of moment tensors.
 * @param[in] lam       [3 x nmt] set of eigenvalues moment tensors.
 *                      The leading dimension is 3.
 *
 * @param[out] gamma    Angle from the double couple meridian to lune point.
 *                      \f$ \gamma \in [-30,30] \f$.  This is an array
 *                      of dimension [nmt].
 * @param[out] delta    Angle from DC meridian to lune point.
 *                      \f$ \delta \in [-90,90] \f$.  This is an array
 *                      of dimension [nmt].
 * @param[out] M0       If M0 is not NULL then this is the seismic moment,
 *                      \f$ M_0 = \frac{|\lam|_2}{\sqrt{2} \f$.
 *                      This is an array or dimension [nmt].
 * @param[out] thetadc  If theta0 is not NULL then this is the angle
 *                      from the double couple meridian to lune point.
 *                      \f$ \theta_{DC} \in [0,90] \f$.  This is an array
 *                      of dimension [nmt].
 * @param[out] lamdev   If lamdev is not NULL then these are the eigenvalues
 *                      of the deviatoric component.  This is an array of
 *                      dimension [3 x nmt] with leading dimension 3.
 * @param[out] lamiso   If lamiso is not NULL then these are the eigenvalues
 *                      of the isotropic component.  This is an array of
 *                      dimension [3 x nmt] with leading dimension 3.
 *
 * @result 0 indicates success.
 *
 * @author Carl Tape and converted to C by Ben Baker.
 *
 * @copyright MIT
 *
 */
int compearth_lam2lune(const int nmt, const double *__restrict__ lam,
                       double *__restrict__ gamma,
                       double *__restrict__ delta,
                       double *__restrict__ M0,
                       double *__restrict__ thetadc,
                       double *__restrict__ lamdev,
                       double *__restrict__ lamiso)
{
    double lamW[3*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double bdot, lam1, lam2, lam3, lammag, lamsum, third_lamsum, xnum, xden;
    int i, ierr, imt, nmtLoc;
    bool lwantLamDev, lwantLamIso, lwantM0, lwantThetaDC;
    const double sqrt3 = sqrt(3.0);
    const double deg = 180.0/M_PI;
    const double third = 1.0/3.0;
    // Error checks
    ierr = 0;
    if (nmt < 1 || lam == NULL || gamma == NULL || delta == NULL)
    {
        if (nmt < 1){fprintf(stderr, "%s: No moment tensors\n", __func__);}
        if (lam == NULL){fprintf(stderr, "%s: lam is NULL\n", __func__);}
        if (gamma == NULL){fprintf(stderr, "%s: gamma is NULL\n", __func__);}
        if (delta == NULL){fprintf(stderr, "%s: delta is NULL\n", __func__);}
        return -1;
    }
    // Figure out the things that I want
    lwantM0 = false;
    if (M0 != NULL){lwantM0 = true;}
    lwantThetaDC = false;
    if (thetadc != NULL){lwantThetaDC = true;}
    lwantLamIso = false;
    if (lamiso != NULL){lwantLamIso = true;}
    lwantLamDev = false;
    if (lamdev != NULL){lwantLamDev = true;}
    for (imt=0; imt<nmt; imt=imt+CE_CHUNKSIZE)
    {
        nmtLoc = MIN(CE_CHUNKSIZE, nmt - imt); 
        // Sort the eigenvalues 
        ierr = compearth_lamsort(nmtLoc, &lam[3*imt], lamW);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error sorting eigenvalues\n", __func__);
            goto ERROR;
        }
        for (i=0; i<nmtLoc; i++)
        {
            lam1 = lamW[3*i];
            lam2 = lamW[3*i+1];
            lam3 = lamW[3*i+2];
            lamsum = lam1 + lam2 + lam3;
            lammag = sqrt(lam1*lam1 + lam2*lam2 + lam3*lam3);
            // Compute Tape and Tape 2012a, Eq 21 (and 23)
            bdot = (lamW[3*i] + lamW[3*i+1] + lamW[3*i+2])/(sqrt3*lammag);
            // Numerical safety 2: is abs(bdot) > 1 -> adjust bdoet to +1 or -1
            bdot = fmax(-1.0, fmin(bdot, 1.0));
            delta[imt+i] = 0.0;
            // Numerical safety 1: if trace(M) == 0 then delta = 0
            if (lamsum != 0.0){delta[imt+i] = 90.0 - acos(bdot)*deg;}
            // Tape and Tape 2012a, Eqn 21a
            xnum =-lamW[3*i] + 2.0*lamW[3*i+1] - lamW[3*i+2]; 
            xden = sqrt3*(lamW[3*i] - lamW[3*i+2]);
            gamma[imt+i] = atan2(xnum, xden)*deg; // N.B. Carl uses atan
            // Extra output; compute seismic moment
            if (lwantM0){M0[imt+i] = lammag*M_SQRT1_2;} 
            // Compute thetadc -- the angle between the DC and the lune point
            if (lwantThetaDC)
            {
                thetadc[imt+i] = acos((lam1 - lam3)/(M_SQRT2*lammag))*deg;
            }
            // Isotropic component of moment tensor
            third_lamsum = third*lamsum;
            if (lwantLamIso)
            {
                lamiso[3*(imt+i)+0] = third_lamsum;
                lamiso[3*(imt+i)+1] = third_lamsum;
                lamiso[3*(imt+i)+2] = third_lamsum;
            }
            if (lwantLamDev)
            {
                lamdev[3*(imt+i)+0] = lamW[3*i+0] - third_lamsum;
                lamdev[3*(imt+i)+1] = lamW[3*i+1] - third_lamsum;
                lamdev[3*(imt+i)+2] = lamW[3*i+2] - third_lamsum;
            }
        }
    } // Loop on MT chunks
ERROR:;
    return ierr;
}
                       
