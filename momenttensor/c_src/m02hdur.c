#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"
/*!
 * @brief Empirical formula of Dahlen and Tromp (1999) pg. 178 for converting
 *        seismic moment to half duration.
 *
 * @param[in] nm     number of scalar moments
 * @param[in] M0     scalar moment in Newton-meters [nm]
 *
 * @param[out] hdur  half duration in seconds [nm]
 *
 * @result 0 indicates success.
 *        -1 indicates invalid inputs.
 *        >0 indicates that some M0's were negative.
 *
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @copyright MIT
 *
 */
int compearth_m02hdur(const int nm, const double *__restrict__ M0,
                      double *__restrict__ hdur)
{
    const char *fcnm = "compearth_m02hdur\0";
    int i, ierr;
    if (nm < 1 || M0 == NULL || hdur == NULL)
    {
        if (nm < 1){printf("%s: No scalar moments\n", fcnm);}
        if (M0 == NULL){printf("%s: M0 is NULL\n", fcnm);}
        if (hdur == NULL){printf("%s: hdur is NULL\n", fcnm);}
        return -1;
    }
    ierr = 0;
    for (i=0; i<nm; i++)
    {
        if (M0[i] < 0.0){ierr = ierr + 1;}
        hdur[i] = 2.4e-6*cbrt(fmax(M0[i], 0.0));
    }
    if (ierr != 0)
    {
        printf("%s: %dscalar moments were negative - hdurs are zero\n",
               fcnm, ierr);
    }
    return ierr;
}
