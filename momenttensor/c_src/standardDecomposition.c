#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#define COMPEARTH_PRIVATE_UPDOWN_ABS_ARGSORT3 1
#include "compearth.h"

/*!
 * @brief A utility for calculating the `standard' isotropic/double couple/CLVD
 *        decomposition of a moment tensor.  This is based on cmopad's standard
 *        decomposition but uses compearth's functionality.
 *
 * @param[in] nmt       Number of moment tensors
 * @param[in] M         Moment tensors to decompose.  This is an array
 *                      of dimension [6 x nmt] with leading dimension 6.
 *                      The bases are given by basis.  Each moment tensor
 *                      is packed \f$ \{m_{11}, m_{22}, m_{33},
 *                                      m_{12}, m_{13}, m_{23} \} \f$.
 *                      The momen tensor terms have units of (N-m).
 * @param[in] basis     Basis of input moment tensors. 
 *
 * @param[out] M0       Scalar moment (N-m) of moment tensors using the
 *                      definition of Hudson and Bowers: M0 = Miso + Mdev.
 *                      This is an array of dimension [nmt].
 * @param[out] Mw       Moment magnitudes computed with the Harvard CMT
 *                      convention:
 *                      \f$ 
 *                         M_w
 *                       = \frac{2}{3} \left ( \log_10(M0) - 16.1 \right )
 *                      \f$.
 *                      This is an array of dimension [nmt].
 * @param[out] fp1      Strike, dip, and rake (respectively) in degrees 
 *                      of fault plane 1.  This is an array of dimension
 *                      [3 x nmt] with leading dimension 3.
 * @param[out] fp2      Strike, dip, and rake (respectively) in degrees
 *                      of fault plane 2.  This is an array of dimension
 *                      [3 x nmt] with leading dimension 3.
 * @param[out] pAxis    This is the azimuth (degrees), plunge (degrees), and
 *                      length (N-m) of the pressure axis.  This is an array
 *                      of dimension [3 x nmt] with leading dimension 3.
 * @param[out] bAxis    This is the azimuth (degrees), plunge (degrees), and
 *                      length (N-m) of the null axis.  This is an array 
 *                      of dimension [3 x nmt] with leading dimension 3.
 * @param[out] tAxis    This is the azimuth (degrees), plunge (degrees), and
 *                      length (N-m) of the tension axis.  This is an array
 *                      of dimension [3 x nmt] with leading dimension 3. 
 * @param[out] isoPct   This is the isotropic percent in range [0,100]
 *                      of the input moment tensors.  This is an array
 *                      of dimension [nmt].
 * @param[out] devPct   This is the deviatoric percent in range [0,100]
 *                      of the input moment tensors. Note, devPct + isoPct
 *                      should equal 100.  This is an array of dimension
 *                      [nmt]. 
 * @param[out] dcPct    This is the double couple percent in range [0,100]
 *                      of the input moment tensors.  This is an array
 *                      of dimension [nmt].
 * @param[out] clvdPct  This is the compensated linear vector dipole [0,100]
 *                      of the input moment tensors.  This is an array
 *                      of dimension [nmt].
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 * @copyright MIT
 *
 */
int compearth_standardDecomposition(const int nmt,
                                    const double *__restrict__ M,
                                    enum compearthCoordSystem_enum basis,
                                    double *__restrict__ M0,
                                    double *__restrict__ Mw,
                                    double *__restrict__ fp1,
                                    double *__restrict__ fp2,
                                    double *__restrict__ pAxis,
                                    double *__restrict__ bAxis,
                                    double *__restrict__ tAxis,
                                    double *__restrict__ isoPct,
                                    double *__restrict__ devPct,
                                    double *__restrict__ dcPct,
                                    double *__restrict__ clvdPct)
{
    double Muse[6*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double MisoW[6*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double MdevW[6*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double M0iso[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double M0devi[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double U[9*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double lam[3*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double lamDevI[3*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double gamma[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double delta[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double M0work[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double kappa1[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double theta1[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double sigma1[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double kappa2[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double theta2[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double sigma2[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double pl1[CE_CHUNKSIZE], pl2[CE_CHUNKSIZE], pl3[CE_CHUNKSIZE];
    double az1[CE_CHUNKSIZE], az2[CE_CHUNKSIZE], az3[CE_CHUNKSIZE];
    double K[3*CE_CHUNKSIZE];
    double S[3*CE_CHUNKSIZE];
    double N[3*CE_CHUNKSIZE];
    double *thetadc = NULL; //[CE_CHUNKSIZE];
    double lamT[3], F;
    int perm[3], i, ierr, imt, nmtLoc;
    const enum magType_enum imag = CE_HARVARD_CMT;
    //const int isort = 1; // sort in descending order
    const double epsilon = 100.0*DBL_EPSILON;
    if (nmt < 1 || M == NULL)
    {
        if (nmt < 1){fprintf(stderr, "%s: No moment tensors\n", __func__);}
        if (M == NULL){fprintf(stderr, "%s: M is NULL\n", __func__);}
    }
    // Chunk the computations to reduce memory demand
    for (imt=0; imt<nmt; imt=imt+CE_CHUNKSIZE)
    {
        nmtLoc = MIN(CE_CHUNKSIZE, nmt - imt);
        // Decompose the moment tensor into a double couple and isotropic part
        ierr = compearth_CMTdecomIso(nmtLoc, &M[6*imt], MisoW, MdevW, NULL);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error in iso/dev decomposition\n", __func__);
            return -1;
        }
        //--------------------------------------------------------------------//
        // There's always this wrinkle of whether or not I want to decompose  //
        // the full moment tensor or the deviatoric + isotropic mts.          //
        // Effectively, what I want to say is that it doesn't matter and      //
        // computing one eigendecomposition is sufficient. That means I want: //
        //                                                                    //
        //    M = X*(m0_iso*I + m_dev)*X' (where X is orthonormal in LAPACK)  //
        //      = X*m0_iso*I*X' + X*m_dev*X'                                  //
        //      = m0_iso X X' + X m_dev X'                                    //
        //      = m0_iso I + X m_dev X'.                                      //
        //                                                                    //
        // Thus, if I have the eigendecomposition of the full tensor          //
        // then I have the eigenvectors of the deviatoric moment tensor.      //
        // Now, this probably needs more analysis if I have a degeneracy      //
        // like a linear vector dipole model but I really should be           //
        // reevaluating my life choices when I want to compute a moment tensor//
        // decomposition of anything that's not substantially double couple.  //
        // Anyway, with the eigenvectors from the full MT or the deviatoric   //
        // MT I can then compute these tension, null, and plunge vectors.     //
        //--------------------------------------------------------------------//
        // Convert M into the working basis for CMT2TT
        ierr = compearth_convertMT(nmtLoc, basis, CE_USE, &M[6*imt], Muse);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: convertMT failed\n", __func__);
            return -1;
        }
        // Get the eigendecomposition and lots of other goodies
        ierr = compearth_CMT2TT(nmtLoc, Muse,
                                false, // ldisply
                                gamma, delta, M0work,
                                kappa1, theta1, sigma1,
                                K, N, S, thetadc, lam, U);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: CMT2TT failed\n", __func__);
            return -1;
        }
        // Compute the deviatoric eigenvalues from the full eignevalues and
        // get the isotropic scalar moment.
        for (i=0; i<nmtLoc; i++)
        {
            // Gather isotropic scalar moment
            M0iso[i] = fabs(MisoW[6*i]);
            // Compute deviatoric eigenvalues by removing isotropic part 
            lamT[0] = lam[3*i+0] - MisoW[6*i];
            lamT[1] = lam[3*i+1] - MisoW[6*i];
            lamT[2] = lam[3*i+2] - MisoW[6*i];
            // Sort deviatoric in ascending order based on absolute value
            argsort3_absUpDown(lamT, true, perm);
            // Apply the argsort
            lamDevI[3*i]   = lamT[perm[0]];
            lamDevI[3*i+1] = lamT[perm[1]];
            lamDevI[3*i+2] = lamT[perm[2]];
        }
        // Represent eigenbasis as plunge/azimuth.
        ierr = compearth_U2pa(nmtLoc, U,
                              pl1, az1, pl2, az2, pl3, az3);
        if (ierr != 0)
        {   
            fprintf(stderr, "%s: Error converting u to plunge/azimuth\n", 
                    __func__);
            return -1;
        }
        // Extract the pressure, null, and tension vectors.
        for (i=0; i<nmtLoc; i++)
        {
            // Eigenvalues are sorted highest to lowest
            lamT[0] = lam[3*i+0] - MisoW[6*i];
            lamT[1] = lam[3*i+1] - MisoW[6*i];
            lamT[2] = lam[3*i+2] - MisoW[6*i];
            // The first eigenvector will go to positive eigenvalue (tension)
            tAxis[3*(imt+i)+0] = az1[i];
            tAxis[3*(imt+i)+1] = pl1[i];
            tAxis[3*(imt+i)+2] = lamT[0];
            // The second eigenvector will be null/neutral (intermediate) 
            bAxis[3*(imt+i)+0] = az2[i];
            bAxis[3*(imt+i)+1] = pl2[i];
            bAxis[3*(imt+i)+2] = lamT[1];
            // The third eigenvector will go to negative eigenvalue (pressure)
            pAxis[3*(imt+i)+0] = az3[i];
            pAxis[3*(imt+i)+1] = pl3[i];
            pAxis[3*(imt+i)+2] = lamT[2];
        }
        // Compute the auxiliary fault plane
        ierr = compearth_auxiliaryPlane(nmtLoc,
                                        kappa1, theta1, sigma1,
                                        kappa2, theta2, sigma2);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error computing auxiliary fault plane\n", 
                    __func__);
            return -1; 
        }
        // I'm going to add some arbitrary priority.  The smaller strike angle
        // gets to be plane 1.
        for (i=0; i<nmtLoc; i++)
        {
            if (kappa1[i] < kappa2[i])
            {
                fp1[3*(imt+i)+0] = kappa1[i];
                fp1[3*(imt+i)+1] = theta1[i];
                fp1[3*(imt+i)+2] = sigma1[i];
                fp2[3*(imt+i)+0] = kappa2[i];
                fp2[3*(imt+i)+1] = theta2[i];
                fp2[3*(imt+i)+2] = sigma2[i];
            }
            else
            {
                fp1[3*(imt+i)+0] = kappa2[i];
                fp1[3*(imt+i)+1] = theta2[i];
                fp1[3*(imt+i)+2] = sigma2[i];
                fp2[3*(imt+i)+0] = kappa1[i];
                fp2[3*(imt+i)+1] = theta1[i];
                fp2[3*(imt+i)+2] = sigma1[i];
            } 
        }
        // Compute the magnitude (Jost and Herrmann Eqn 19)
        for (i=0; i<nmtLoc; i++)
        {
            M0devi[i] = 0.5*(fabs(lamDevI[3*i+1]) + fabs(lamDevI[3*i+2]));
        }
        for (i=0; i<nmtLoc; i++)
        {
            M0[imt+i] = M0iso[i] + M0devi[i]; // scalar moment Bowers and Hudson
            isoPct[imt+i] = (M0iso[i]/M0[imt+i])*100.0;
            F = 0.5;
            if (M0devi[i] > epsilon){F =-lamDevI[3*i+0]/lamDevI[3*i+2];}
            dcPct[imt+i] = (1.0 - 2.0*fabs(F))*(1.0 - 0.01*isoPct[imt+i])*100.0;
            devPct[imt+i] = 100.0 - isoPct[imt+i];
            clvdPct[imt+i] = 100.0 - isoPct[imt+i] - dcPct[imt+i];
        }
        // Compute the moment magnitude
        ierr = compearth_m02mw(nmtLoc, imag, &M0[imt], &Mw[imt]);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: m02mw failed\n", __func__);
            return -1;
        }
    } // Loop on MT chunks
    return 0;
}
