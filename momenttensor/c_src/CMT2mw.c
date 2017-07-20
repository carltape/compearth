#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief Convert from moment tnesor to moment magnitude 
 *
 * @param[in] im0    =1 for Silver and Jordan (1982) formula.
 *                   =2 for GCMT formula.
 *                   =3 for old `Caltech Way'
 * @param[in] M      [nm x 6] array of moment tensors packed
 *                   {M11, M22, M33, M12, M13, M23} with leading
 *                   dimension 6.  The Mij element should be of units
 *                   of Newton-meters.
 *
 * @param[out] Mw    moment magnitudes [nw]
 *
 * @result 0 indicates success
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
    const char *fcnm = "compearth_CMT2mw\0";
    int ierr;
    double *M0, M064[64];
    if (nm > 64)
    {
        M0 = (double *) calloc((size_t) nm, sizeof(double));
    }
    else
    {
        M0 = M064;
    }
    ierr = compearth_CMT2m0(nm, im0, M, M0);
    if (ierr != 0)
    {
        printf("%s: Error computing m0\n", fcnm);
        if (nm > 64){free(M0);}
        return -1;
    }
    ierr = compearth_m02mw(nm, KANAMORI_1978, M0, Mw);
    if (ierr != 0)
    {
        printf("%s: Error computing Mw\n", fcnm);
        return -1;
    }
    if (nm > 64){free(M0);}
    M0 = NULL;
    return 0;
}
