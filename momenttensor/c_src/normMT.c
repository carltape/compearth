#include <stdio.h>
#include <stdlib.h>
#include "compearth.h"
/*!
 * @brief Comptues matrix (Frobenius) onrm for symmetric matrix (moment tensor).
 *
 * @param[in] n       Number of moment tensors.
 * @param[in] M       Input symmetric matrix: {M11, M22, M33, M12, M13, M23}.
 *                    This is an array of dimension [n x 6] with leading
 *                    dimension 6.
 * @param[in] Lnorm   Matrix norm: \n
 *                      TWO_NORM (2). \n
 *                      ONE_NORM (1). \n
 *                      P_NORM (in this case must set p). \n
 *                      INFINITY_NORM. \n
 *                      NEGATIVE_INFINITY_NORM. \n
 * @param[in] p       If using a p-norm this is the value for p (> 0).
 *
 * @param[out] Mnorm  Matrix norms for each moment tensor.  This is an array
 *                    of dimension [n].
 *
 * @result 0 indicates success.
 *
 */
int compearth_normMT(const int n,
                     const double *M,
                     const enum normType_enum Lnorm,
                     const double p,
                     double *mnorm)
{
    const char *fcnm = "compearth_normMT\0";
    double Mn[9];
    int i, ierr;
    ierr = 0;
    for (i=0; i<n; i++)
    {
        // Transform to 3 x 3 matrix
        compearth_Mvec2Mmat(1, &M[6*i], 1, Mn);
        // Compute norm
        ierr = compearth_normMat(1, Mn, Lnorm, p, &mnorm[i]);
        if (ierr != 0)
        {
            printf("%s: Error computing matrix norm!\n", fcnm);
        }
    }
    return ierr;
}
