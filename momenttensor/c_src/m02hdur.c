#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"
/*!
 * @brief Empirical formula of Dahlen and Tromp (1999) pg. 178 for converting
 *        seismic moment to half duration.
 *
 * @param[in] nm     Number of scalar moments.
 * @param[in] M0     Scalar moment in Newton-meters.  This is an array
 *                   of dimension [nm].
 *
 * @param[out] hdur  Half duration in seconds.  This is an array of
 *                   dimension [nm].
 *
 * @result 0 indicates success. \n
 *        -1 indicates invalid inputs. \n
 *        >0 indicates that some M0's were negative.
 *
 * @author Carl Tape and translated to C by Ben Baker.
 *
 * @copyright MIT
 *
 */
int compearth_m02hdur(const int nm, const double *__restrict__ M0,
                      double *__restrict__ hdur)
{
    int i, ierr;
    if (nm < 1 || M0 == NULL || hdur == NULL)
    {
        if (nm < 1){fprintf(stderr, "%s: No scalar moments\n", __func__);}
        if (M0 == NULL){fprintf(stderr, "%s: M0 is NULL\n", __func__);}
        if (hdur == NULL){fprintf(stderr, "%s: hdur is NULL\n", __func__);}
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
        fprintf(stderr, "%s: %d scalar moments were negative; hdurs are zero\n",
                __func__, ierr);
    }
    return ierr;
}
