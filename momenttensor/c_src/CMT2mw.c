#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief Convert from moment tnesor to moment magnitude 
 *
 * @param[in] nm     Number of moment tensors.
 * @param[in] im0    =1 for Silver and Jordan (1982) formula. \n
 *                   =2 for GCMT formula. \n
 *                   =3 for old `Caltech Way.'
 * @param[in] M      [nm x 6] array of moment tensors packed
 *                   {M11, M22, M33, M12, M13, M23} with leading
 *                   dimension 6.  The Mij element should be of units
 *                   of Newton-meters.
 *
 * @param[out] Mw    Moment magnitudes.  This is an array of dimension [nm].
 *
 * @result 0 indicates success.
 * 
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @copyright MIT
 *
 */
int compearth_CMT2mw(const int nm, const int im0,
                     const double *__restrict__ M,
                     double *__restrict__ Mw)
{
    double M0[CE_CHUNKSIZE] __attribute__((aligned(64)));
    int i, ierr, nmLoc;
    // Loop on moment tensor chunks
    for (i=0; i<nm; i=i+CE_CHUNKSIZE)
    {
        nmLoc = MIN(CE_CHUNKSIZE, nm - i);
        ierr = compearth_CMT2m0(nmLoc, im0, &M[i], M0);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error computing m0\n", __func__);
            return -1;
        }
        ierr = compearth_m02mw(nmLoc, CE_KANAMORI_1978, M0, &Mw[i]);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error computing Mw\n", __func__);
            return -1;
        }
    } // Loop on chunks
    return 0;
}
