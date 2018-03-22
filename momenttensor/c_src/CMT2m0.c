#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compearth.h"
#ifdef COMPEARTH_USE_MKL
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wstrict-prototypes"
#endif
#include <mkl_lapacke.h>
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#else
#include <lapacke.h>
#endif

/*!
 * @brief Convert from moment tnesor to scalar seismic moment.
 *
 * @param[in] nm     Number of moment tensors.
 * @param[in] im0    =1 for Silver and Jordan (1982) formula. \n
 *                   =2 for GCMT formula. \n
 *                   =3 for old `Caltech Way.'
 * @param[in] M      [nm x 6] array of moment tensors packed
 *                   {M11, M22, M33, M12, M13, M23} with leading
 *                   dimension 6.  The Mij element should be of units
 *                   of Newton-meters (although it should not affect
 *                   the function here).
 *
 * @param[out] M0    Scalar moments with units of Newton-meters.  This
 *                   is an array of dimension [nm].
 *
 * @result 0 indicates success.
 * 
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @copyright MIT
 *
 */
int compearth_CMT2m0(const int nm, const int im0,
                     const double *__restrict__ M,
                     double *__restrict__ M0)
{
    double Mt9[9], lams[3], work[100], M11, M22, M33, M12, M13, M23;
    const double sqrt2i = 1.0/sqrt(2.0);
    const int lwork = 100;
    int i, ierr, info;
    ierr = 0;
    if (im0 == 2 || im0 == 3)
    {
        for (i=0; i<nm; i++)
        {
            // Fill upper-triangle; column 1
            Mt9[0] = M[6*i];
            Mt9[1] = 0.0;
            Mt9[2] = 0.0;
            // column 2 
            Mt9[3] = M[6*i+3];
            Mt9[4] = M[6*i+1];
            Mt9[5] = 0.0;
            // column 3
            Mt9[6] = M[6*i+4];
            Mt9[7] = M[6*i+5];
            Mt9[8] = M[6*i+2];
            // get eigenvalues in ascending order
            info = LAPACKE_dsyev_work(LAPACK_COL_MAJOR, 'N', 'U', 3,
                                      Mt9, 3, lams, work, lwork);
            if (info != 0)
            {
                fprintf(stderr, "%s: Error computing eigenvalues\n", __func__);
                return -1;
            }
            // GCMT double couple formulation
            if (im0 == 2)
            {
                M0[i] = 0.5*(fabs(lams[0]) + fabs(lams[2]));
            }
            else
            {
                M0[i] = fabs(lams[2]);
                if (fabs(lams[0]) > fabs(lams[2])){M0[i] = fabs(lams[0]);}
            }
        }
    }
    // Silver and Jordan
    else if (im0 == 1)
    {
        for (i=0; i<nm; i++)
        {
            M11 = M[6*i];
            M22 = M[6*i+1];
            M33 = M[6*i+2];
            M12 = M[6*i+3];
            M13 = M[6*i+4];
            M23 = M[6*i+5];
            M0[i] = sqrt2i*sqrt((      M11*M11 + M22*M22 + M33*M33
                                + 2.0*(M12*M12 + M13*M13 + M23*M23)));
        }
    }
    else
    {
        fprintf(stderr, "%s: Invalid magnitude type %d\n", __func__, im0); 
        for (i=0; i<nm; i++){M0[i] = 0.0;}
        return -1;
    }
    return ierr; 
}
