#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#define COMPEARTH_PRIVATE_UPDOWN_ABS_ARGSORT3 1
#include "compearth.h"
#ifdef COMPEARTH_USE_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

/*!
 * @brief A utility for calculating the `standard' decomposition of a moment 
 *        tensor.  This is based on cmopad's standard decomposition but
 *        uses compearth's functionality.
 *
 * @param[in] nmt     Number of moment tensors
 * @param[in] M       Moment tensors to decompose.  This is an array
 *                    of dimension [6 x nmt] with leading dimension 6.
 *                    The bases are given by basis.  Each moment tensor
 *                    is packed \f$ {m_{11}, m_{22}, m_{33},
 *                                   m_{12}, m_{13}, m_{23} \} \f$.
 * @param[in] basis   Basis of input moment tensors. 
 *
 */
int compearth_standardDecomposition(const int nmt,
                                    const double *__restrict__ M,
                                    enum compearthCoordSystem_enum basis,
                                    double *__restrict__ M0,
                                    double *__restrict__ pAxis,
                                    double *__restrict__ bAxis,
                                    double *__restrict__ tAxis,
                                    double *__restrict__ isoPct,
                                    double *__restrict__ dcPct,
                                    double *__restrict__ devPct,
                                    double *__restrict__ clvdPct)
{
    double Muse[6*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double MisoW[6*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double MdevW[6*CE_CHUNKSIZE] __attribute__((aligned(64)));
    //double MdevUSE[6*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double M0iso[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double M0devi[CE_CHUNKSIZE] __attribute__((aligned(64)));
    double U[9*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double lam[3*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double lamDevI[3*CE_CHUNKSIZE] __attribute__((aligned(64)));
    double gamma[CE_CHUNKSIZE];
    double delta[CE_CHUNKSIZE];
    double M0work[CE_CHUNKSIZE];
    double kappa1[CE_CHUNKSIZE];
    double theta1[CE_CHUNKSIZE];
    double sigma1[CE_CHUNKSIZE]; 
    double kappa2[CE_CHUNKSIZE];
    double theta2[CE_CHUNKSIZE];
    double sigma2[CE_CHUNKSIZE];
    double pl1[CE_CHUNKSIZE], pl2[CE_CHUNKSIZE], pl3[CE_CHUNKSIZE];
    double az1[CE_CHUNKSIZE], az2[CE_CHUNKSIZE], az3[CE_CHUNKSIZE];
    double K[3*CE_CHUNKSIZE];
    double S[3*CE_CHUNKSIZE];
    double N[3*CE_CHUNKSIZE];
    double *thetadc = NULL; //[CE_CHUNKSIZE];
    double lamT[3], F;
    int perm[3], i, ierr, imt, nmtLoc;
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
        // TODO remove this line
        //ierr = compearth_CMTdecom(nmtLoc, Muse, isort, lam, U);
        // Compute the deviatoric eigenvalues from the full eignevalues and
        // get the isotropic scalar moment.
        for (i=0; i<nmtLoc; i++)
        {
            // Gather isotropic scalar moment
            M0iso[i] = MisoW[6*i];
            // Compute deviatoric eigenvalues by removing isotropic part 
            lamT[0] = lam[3*i+0] - M0iso[i];
            lamT[1] = lam[3*i+1] - M0iso[i];
            lamT[2] = lam[3*i+2] - M0iso[i];
            // Sort deviatoric in ascending order based on absolute value
            argsort3_absUpDown(lamT, true, perm);
            // Apply the argsort
            lamDevI[3*i]   = lamT[perm[i]];
            lamDevI[3*i+1] = lamT[perm[i+1]];
            lamDevI[3*i+2] = lamT[perm[i+2]];
        }
/*
        // Compute the isotropic magnitude: 1/3*trace(M_{iso}).  However,
        // M_{iso} has constant diagonal so it is sufficient to simply grab
        // an element instead of sum the digonal and divide by 3. 
        //cblas_dcopy(nmtLoc, MisoW, 3, M0iso, 1); 
        // The computation proceeds with the deviatoric moment tensor in SEU 
        // coordinates
        ierr = compearth_convertMT(nmtLoc, basis, CE_SEU, MdevW, MdevUSE);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to convert MT basis\n", __func__);
            return -1;
        }
        // Decompose the deviatoric moment tensor
        ierr = compearth_CMTdecom(nmtLoc, MdevUSE, isort, lam, U);
        // Sort the eigenvalues in descending order of absolute value
        for (i=0; i<nmtLoc; i++)
        {
            argsort3_absUpDown(&lam[3*i], true, perm);
            lamDevI[3*i]   = lam[3*i + perm[i]];
            lamDevI[3*i+1] = lam[3*i + perm[i+1]];
            lamDevI[3*i+2] = lam[3*i + perm[i+2]];
        }
*/
double p[3], b[3], t[3];
        // Represent eigenbasis as plunge/azimuth.
        // TODO: An SVD lurks in here - I should probably make the
        // orthogonalization style optional b/c it' probably not worth it.
        ierr = compearth_U2pa(nmtLoc, U,
                              pl1, az1, pl2, az2, pl3, az3); //&p[2], &p[1], &b[2], &b[1], &t[2], &t[1]);
        if (ierr != 0)
        {   
            fprintf(stderr, "%s: Error converting u to plunge/azimuth\n", 
                    __func__);
            return -1;
        }
        // Extract the plunge, null, and tension vectors
        for (i=0; i<nmtLoc; i++)
        {
            pAxis[3*i+0] = pl1[i];
            pAxis[3*i+1] = az1[i]; 
            pAxis[3*i+1] = lam[3*i];

            bAxis[3*i+0] = pl2[i];
            bAxis[3*i+1] = az2[i]; 
            bAxis[3*i+1] = lam[3*i+1];
 
            tAxis[3*i+0] = pl3[i];
            tAxis[3*i+1] = az3[i]; 
            tAxis[3*i+1] = lam[3*i+2];
        }
        // Fill in the eigenvalues.  They were sorted in descending order. 
        cblas_dcopy(nmtLoc, &lam[0], 3, p, 3); // p[0] = lam[0];
        cblas_dcopy(nmtLoc, &lam[1], 3, b, 3); // b[0] = lam[1];
        cblas_dcopy(nmtLoc, &lam[2], 3, t, 3); // t[0] = lam[2];
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

        }
    printf("sdr 1: %f %f %f\n", kappa1[0], theta1[0], sigma1[0]);
    printf("sdr 2: %f %f %f\n", kappa2[0], theta2[0], sigma2[0]);

        // Compute the magnitude (Jost and Herrmann Eqn 19)
        for (i=0; i<nmtLoc; i++)
        {
            M0devi[i] = 0.5*(fabs(lamDevI[3*i+1]) + fabs(lamDevI[3*i+2]));
        }
        for (i=0; i<nmtLoc; i++)
        {
            M0[imt+i] = M0iso[i] + M0devi[i]; // scalar moment Bowers and Hudson
            isoPct[imt+i] = M0iso[i]/M0[imt+i]*100.0;
            F = 0.5;
            if (M0devi[i] > epsilon){F =-lamDevI[3*i+0]/lamDevI[3*i+2];}
            dcPct[imt+i] = (1.0 - 2.0*fabs(F))*(1.0 - 0.01*isoPct[i])*100.0;
            devPct[imt+i] = 100.0 - isoPct[imt+i];
            clvdPct[imt+i] = 100.0 - isoPct[imt+i] - dcPct[imt+i];
        }

 
    }
    return 0;
}
