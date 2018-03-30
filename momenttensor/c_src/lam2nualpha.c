#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief Converts eigenvalues to \f$ \nu \f$ and \f$ \alpha \f$.
 *
 * @param[in] n       Number of eigenvalue triples.
 * @param[in] lam     [3 x n] set of eigenvalue triples.  The leading dimension
 *                    is 3.
 * @param[out] nu     Unitless Poisson parameters.  This is an array of
 *                    dimension n.
 * @param[out] alpha  Angles between fault normal and slip vector
 *                    where \f$ \alpha \in [0,180] \f$.  This is an array
 *                    of dimension n.
 *
 * @result 0 indicates success.
 *
 * @author Carl Tape and translated to C by Ben Baker
 *
 * @copyright MIT
 *
 */
int compearth_lam2nualpha(const int n, const double *__restrict__ lam,
                          double *__restrict__ nu, double *__restrict__ alpha)
{
    double lam1, lam2, lam3;
    int i;
    const double deg = 180.0/M_PI;
    if (n < 1 || lam == NULL || alpha == NULL || nu == NULL)
    {
        if (n < 1){fprintf(stderr, "%s: Error no moment tensors\n", __func__);}
        if (lam == NULL){fprintf(stderr, "%s: Error lam is NULL\n", __func__);}
        if (alpha == NULL)
        {
            fprintf(stderr, "%s: Error alpha is NULL\n", __func__);
        }
        if (nu == NULL){fprintf(stderr, "%s: Error nu is NULL\n", __func__);} 
        return -1;
    }
    for (i=0; i<n; i++)
    {
        lam1 = lam[3*i];
        lam2 = lam[3*i+1];
        lam3 = lam[3*i+2];
        alpha[i] = acos((lam1 - 2.0*lam2 + lam3)/(lam1 - lam3))*deg;
        nu[i] = lam2/(lam1 + lam3);
    }
    return 0;
}
