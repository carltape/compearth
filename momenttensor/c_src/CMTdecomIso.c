#include <stdio.h>
#include <stdlib.h>
#include "compearth.h"
/*!
 * @brief Decomposes a general moment tensor into its isotropic and 
 *        deviatoric components.
 *
 * @param[in] n      Number of moment tensors.
 * @param[in] M      [6 x n] general moment tensors.  The leading dimension
 *                   is 6.
 * @param[out] Miso  [6 x n] isotropic moment tensors.  The leading dimension
 *                   is 6.
 * @param[out] Mdev  [6 x n] deviatoric moment tensors.  The leading dimension
 *                   is 6.
 * @param[out] trM   trace(M).  This is an array of dimension [n].
 *
 * @result 0 indicates success.
 *
 * @copyright MIT
 *
 */
int compearth_CMTdecomIso(const int n, const double *__restrict__ M,
                          double *__restrict__ Miso,
                          double *__restrict__ Mdev,
                          double *__restrict__ trM)
{
    const double third = 1.0/3.0;
    int i;
    if (n < 1 || M == NULL || Miso == NULL || Mdev == NULL || trM == NULL)
    {
        if (n < 1){fprintf(stderr, "%s: No moment tensors\n", __func__);}
        if (M == NULL){fprintf(stderr, "%s: M is NULL\n", __func__);}
        if (Miso == NULL){fprintf(stderr, "%s: Miso is NULL\n", __func__);}
        if (Mdev == NULL){fprintf(stderr, "%s: Mdev is NULL\n", __func__);}
        if (trM == NULL){fprintf(stderr, "%s: trM is NULL\n", __func__);}
        return -1;
    }
    for (i=0; i<n; i++)
    {
        trM[i] = M[6*i] + M[6*i+1] + M[6*i+2];
        Miso[6*i+0] = third*trM[i];
        Miso[6*i+1] = third*trM[i];
        Miso[6*i+2] = third*trM[i];
        Miso[6*i+3] = 0.0;
        Miso[6*i+4] = 0.0;
        Miso[6*i+5] = 0.0;
    }
    for (i=0; i<6*n; i++)
    {
        Mdev[i] = M[i] - Miso[i]; 
    }
    return 0;
}
