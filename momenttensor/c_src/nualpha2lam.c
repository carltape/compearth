#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief Converts \f$ \nu \f$ and \f$ \alpha \f$ to unit \f$ \Lambda \f$
 *        vector.  This is hte inverse function of lam2nualpha.  See
 *        Tape and Tape (2013), "The classical model for moment tensors".
 *
 * @param[in] n       Number of eigenvalue triples.
 * @param[in] nu      Unitless Poisson parameters.  This is an array of
 *                    dimension n.
 * @param[in] alpha   Angles between fault normal and slip vector
 *                    where \f$ \alpha \in [0,180] \f$.  This is an array
 *                    of dimension n.
 *
 * @param[out] lam    [3 x n] set of eigenvalue triples.  The leading dimension
 *                    is 3.
 *
 * @author Carl Tape and converted to C by Ben Baker
 *
 * @copyright MIT
 *
 */
int compearth_nualpha2lam(const int n,
                          const double *__restrict__ nu,
                          const double *__restrict__ alpha,
                          double *__restrict__ lam)
{
    double cos_1p2n, cosa, lam1, lam2, lam3, mag;
    int i;
    const double rad = M_PI/180.0;
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
        cosa = cos(alpha[i]*rad);
        // TapeTape2013 Eq 30
        cos_1p2n = cosa/(1.0 - 2.0*nu[i]);
        lam1 = cos_1p2n + 1.0;
        lam2 = 2.*nu[i]*cos_1p2n;
        lam3 = cos_1p2n - 1.0;
        // Unit norm
        mag = sqrt(lam1*lam1 + lam2*lam2 + lam3*lam3);
        lam[3*i+0] = lam1/mag;
        lam[3*i+1] = lam2/mag;
        lam[3*i+2] = lam3/mag; 
    }
    return 0;
}
