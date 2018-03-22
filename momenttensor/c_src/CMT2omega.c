#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

/*!
 * @brief Computes the omega angles between two moment tensors.
 *
 * @param[in] nmt1    Number of moment tensors in M1.  This can be 1 otherwise
 *                    it must be equal to nmt2.  If it is 1 then the angle will
 *                    be between M1 and all the moment tensors in M2.
 * @param[in] M1      First array of moment tensors.  M1 is packed
 *                    \f$ M = \{M_{11} M_{22} M_{33} M_{12} M_{13} M_{23}\} \f$.
 *                    This is an array of dimension [6 x nmt1] with leading
 *                    dimension 6.
 * @param[in] nmt2    Number of moment tensors in M2.  This can be 1 otherwise
 *                    it must be equal to nmt1.  If it is 1 then the angle will
 *                    be between M2 and all the moment tensors in M1.
 * @param[in] M2      Second array of moment tensors.  M2 is packed
 *                    \f$ M = \{M_{11} M_{22} M_{33} M_{12} M_{13} M_{23}\} \f$.
 *                    This is an array of dimension [6 x nmt2] with leading
 *                    dimension 6.
 *
 * @param[out] omega  This is a measure between two moment tensors and combines
 *                    differences in eigenvalues (source types) and orientation
 *                    but not magnitude.
 *                    This is an array of dimension [MAX(nmt1, nmt2)].
 *
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @copyright MIT
 *
 */
int compearth_CMT2omega(const int nmt1, const double *__restrict__ M1,
                        const int nmt2, const double *__restrict__ M2,
                        double *__restrict__ omega)
{
    double Mwork[6*CE_CHUNKSIZE] __attribute__((aligned(64)));
    int ierr, imt, nmtLoc;
    const double pi180i = 180.0/M_PI;
    ierr = 0;
    if (nmt1 < 1 || M1 == NULL || nmt2 < 1 || M2 == NULL || omega == NULL)
    {
        if (nmt1 < 1)
        {
            fprintf(stderr, "%s: Error nmt1 must be positive\n", __func__);
        }
        if (nmt2 < 1)
        {
            fprintf(stderr, "%s: Error nmt2 must be positive\n", __func__);
        }
        if (M1 == NULL){fprintf(stderr, "%s: Error M1 is NULL\n", __func__);}
        if (M2 == NULL){fprintf(stderr, "%s: Error M2 is NULL\n", __func__);}
        if (omega == NULL)
        {
            fprintf(stderr, "%s: Error omega is NULL\n", __func__);
        }
        return -1;
    }
    // Straight compuation
    if (nmt1 == nmt2)
    {
        ierr = compearth_angleMT(nmt1, M1, M2, omega);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error computing angleMT 1\n", __func__);
        }
        // Convert to degrees
        cblas_dscal(nmt1, pi180i, omega, 1);
        return ierr; 
    }
    else
    {
        if (nmt1 == 1 || nmt2 == 1)
        {
            // Copy nmt1 and repeat computation with it
            if (nmt1 == 1)
            {
                for (imt=0; imt<CE_CHUNKSIZE; imt++)
                {
                    cblas_dcopy(6, &M1[0], 1, &Mwork[6*imt], 1);
                }
                for (imt=0; imt<nmt2; imt=imt+CE_CHUNKSIZE)
                {
                    nmtLoc = MIN(CE_CHUNKSIZE, nmt2 - imt);
                    ierr = compearth_angleMT(nmtLoc, Mwork, &M2[6*imt],
                                             &omega[imt]); 
                    if (ierr != 0)
                    {   
                        fprintf(stderr, "%s: Error calling angleMT at %d\n",
                                __func__, __LINE__);
                        return ierr;
                    }
                    // Convert to degrees
                    cblas_dscal(nmtLoc, pi180i, &omega[imt], 1);
                }
            }
            // Copy nmt2 and repeat computation with it
            else
            {
                for (imt=0; imt<CE_CHUNKSIZE; imt++)
                {   
                    cblas_dcopy(6, &M2[0], 1, &Mwork[6*imt], 1); 
                }
                for (imt=0; imt<nmt1; imt=imt+CE_CHUNKSIZE)
                {
                    nmtLoc = MIN(CE_CHUNKSIZE, nmt1 - imt);
                    ierr = compearth_angleMT(nmtLoc, &M1[6*imt], Mwork,
                                             &omega[imt]); 
                    if (ierr != 0)
                    {
                        fprintf(stderr, "%s: Error calling angleMT at %d\n",
                                __func__, __LINE__);
                        return ierr;
                    }
                    // Convert to degrees
                    cblas_dscal(nmtLoc, pi180i, &omega[imt], 1);
                }
            }
        }
        else
        {
            fprintf(stderr, "%s: nmt1=%d /= nmt2=%d; this is invalid\n",
                    __func__, nmt1, nmt2);
            return -1;
        }
    }
    return ierr; 
}
