#include <stdio.h>
#include <stdlib.h>
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

static double compearth_norm64f(const int n, const double *__restrict__ x,
                                const enum ceNormType_enum norm,
                                const double p, int *ierr);
/*!
 * @brief Computes the matrix (Frobenius) norm for a matrix
 *
 * @param[in] n       number of matrices
 * @param[in] M       3 x 3 input matrix as a length 9 array [9*n]
 * @param[in] Lnorm   matrix norm
 *                      TWO_NORM (2)
 *                      ONE_NORM (1)
 *                      P_NORM (in this case must set p)
 *                      INFINITY_NORM
 *                      NEGATIVE_INFINITY_NORM
 * @param[in] p       if using a p-norm this is the value for p (> 0)
 *
 * @param[out] mnorm  matrix norms for each matrix
 *
 * @result 0 indicates success 
 *
 */
int compearth_normMat(const int n,
                      const double *__restrict__ M,
                      const enum ceNormType_enum Lnorm,
                      const double p,
                      double *__restrict__ mnorm)
{
    int i, ierr;
    ierr = 0;
    for (i=0; i<n; i++)
    {
        mnorm[i] = compearth_norm64f(9, &M[9*i], Lnorm, p, &ierr);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error computing matrix norm!\n", __func__);
            mnorm[i] = 0.0;
            break;
        }
    }
    return ierr;
}

/*!
 * @brief Computes the norm of a vector.  This is from ISTI's ISCL.
 *
 * @param[in] n      Length of array x.
 * @param[in] x      Array of dimension [n] of which to compute norm.
 * @param[in] norm   Type of norm to compute.
 * @param[in] p      If performing a P norm then this must be defined to
 *                   a real number greater than or equal to 1.  Otherwise,
 *                   it will not be accessed.
 *
 * @param[out] ierr  0 indicates success
 *
 * @result P-norm of array x.
 *
 * @author Ben Baker
 *
 * @date August 2017 - redefined variables so that this may exist in 
 *       compearth without collisions with ISCL.
 * 
 */
static double compearth_norm64f(const int n, const double *__restrict__ x,
                                const enum ceNormType_enum norm,
                                const double p, int *ierr)
{
    double xnorm;
    int i;
    *ierr = 0;
    xnorm = 0.0;
    if (n < 1 || x == NULL)
    {
        if (n < 1){fprintf(stderr, "%s: Error no data\n", __func__);}
        if (x == NULL){fprintf(stderr, "%s: Error x is NULL\n", __func__);}
        *ierr = 1;
        return xnorm;
    }
    // 2 norm
    if (norm == CE_TWO_NORM)
    {
        xnorm = cblas_dnrm2(n, x, 1);
    // 1 norm
    }
    else if (norm == CE_ONE_NORM)
    {
        xnorm = cblas_dasum(n, x, 1);
    }
    // p norm
    else if (norm == CE_P_NORM)
    {
        if (p <= 0.0)
        {
            fprintf(stderr, "%s: Invalid p value %f\n", __func__, p);
            *ierr = 1;
            return xnorm;
        }
        #pragma omp simd reduction(+:xnorm)
        for (i=0; i<n; i++)
        {
            xnorm = xnorm + pow(fabs(x[i]), p);
        }
        xnorm = pow(xnorm, 1./p);
    }
    // infinity norm
    else if (norm == CE_INFINITY_NORM)
    {
        xnorm = fabs(x[cblas_idamax(n, x, 1)]);
    }
    // negative infinity norm
    else if (norm == CE_NEGATIVE_INFINITY_NORM)
    {
        xnorm = fabs(x[0]);
        #pragma omp simd reduction(min:xnorm)
        for (i=1; i<n; i++)
        {
            xnorm = fmin(fabs(x[i]), xnorm);
        }
    }
    // not sure - default to 2-norm
    else
    {
        fprintf(stderr, "%s: Defaulting to 2-norm\n", __func__);
        xnorm = cblas_dnrm2(n, x, 1);
    }
    return xnorm;
}
